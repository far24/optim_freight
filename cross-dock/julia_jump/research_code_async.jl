#*****************************************
# Necessary Package
#*****************************************
using Random
using JuMP
using GLPK, Cbc
using LinearAlgebra
using CSV
using Distances

# vehicle routing problem with cross-dock
# Wen et al. (2009)

#*****************************************
# Network Property
#*****************************************
# suppliers are the pick-up nodes
P = [1,2]

# customers are the delivery nodes
D = [3,4]

# cross-docks
cd_pick_start = [D[end]+1]
cd_pick_end = [D[end]+2]
cd_del_start = [D[end]+3]
cd_del_end = [D[end]+4]



N = vcat(P,D, cd_pick_start, cd_pick_end, cd_del_start, cd_del_end)

# vehicles
K = [1]

# ****************************************
# parameters
#*****************************************
# number of pickup/delivery Nodes
n = length(P)

d_i = [10,20,10,20, 1,1,1,1]
veh_cap = 40

# location of nodes on a graph
#c = [9,5,4,8,13,11,8,8]
c = [90,50,150,130,80,80,80,80]
t_ij = pairwise(Euclidean(), c'; dims=2)
t_ij

#=
# [ai , bi] = the time window for node i (i ∈ N);
a_i = zeros(N)
b_i = rand(5:10, N)
# di = the amount of demand of request i (i ∈ P);
d_i_pickup = rand(5:10, P)
d_i_delivery = d_i_pickup
d_i = vcat(d_i_pickup, d_i_delivery)

# Q = the vehicle capacity;
Q = 50
# A = the fixed time for unloading and reloading at the cross-dock;
A = rand(5:10, O)
# B = the time for unloading and reloading a pallet.
A = rand(5:10, O)
=#


# very big number
M= 99999

# ****************************************
# Create a JuMP model
# ****************************************
cd_modl = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************

# BINARY
# ^^^^^^^^^
# x_ijk = 1 if vehicle k travels from node i to node j ((i, j ) ∈ E, k ∈ K)
# u_ik = 1 if vehicle k unloads demand/request i at the cross-dock (i ∈ P, k ∈ K); 0 otherwise;
# r_ik =1 if vehicle k reloads demand/request i at the cross-dock (i ∈ P, k ∈ K); 0 otherwise;
# q_k = 1 if vehicle k has to unload at the cross-dock (k ∈ K); 0 otherwise;
# qq_k = 1 if vehicle k has to reload at the cross-dock (k ∈ K); 0 otherwise;
@variable(cd_modl, x_ijk[i=N,j=N,k=K], Bin)
@variable(cd_modl, u_ik[i=P, k=K], Bin)
@variable(cd_modl, r_ik[i=P, k=K], Bin)
@variable(cd_modl, g_k[k=K], Bin)
@variable(cd_modl, h_k[k=K], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# l_ik = the time at which vehicle k leaves node i (i ∈ N, k ∈ K); 0 otherwise;
# f_k = finish time of unloading operation of truck k at cross-dock (k ∈ K)
# g_i = finish time of unloading operation of request i at cross-dock (i ∈ P)
# w_k = start time of reloading operation of truck k at cross-dock (k ∈ K)
@variable(cd_modl, s_ik[i=N, k=K] >=0 )
@variable(cd_modl, f_k[k=K]>= 0)
@variable(cd_modl, w_k[k=K] >= 0)
@variable(cd_modl, v_i[i=P] >= 0)


# ****************************************
# Constraints
#*****************************************

# --------------------[Pick-up Process]
# 01: only one vehicle can arrive at a pickup node
@constraint(cd_modl, p_nd_exit[i= P],
                    sum(x_ijk[i,j,k] for j=vcat(P, cd_pick_end) for k=K if i != j) == 1)


# 02: start from cross dock // all the trucks must be utilized
@constraint(cd_modl, p_nd_start[i= cd_pick_start, k=K],
                    sum(x_ijk[i,j,k] for j=P if i != j) == 1)

# 03: end at cross dock
@constraint(cd_modl, p_nd_end[j= cd_pick_end, k=K],
                sum(x_ijk[i,j,k] for i=P if i != j) == 1)

# 04: capacity of vehicle
@constraint(cd_modl, p_veh_cap[k=K],
                sum(d_i[i] * x_ijk[i,j,k] for i=P for j=vcat(P,cd_pick_end) if i != j) <= veh_cap)

# 05: continuous routes/ consecutive movement
@constraint(cd_modl, p_move[h=P, k=K],
                    sum(x_ijk[i,h,k] for i=vcat(cd_pick_start,P) if i != h)
                    - sum(x_ijk[h,j,k] for j=vcat(P,cd_pick_end) if j != h) == 0)


# --------------------[delivery Process]
# 01: only one vehicle can arrive at a pickup node
@constraint(cd_modl, d_nd_exit[i=D],
                    sum(x_ijk[i,j,k] for j=vcat(D, cd_del_end) for k=K if i != j) == 1)

# 02: start from cross dock // all the trucks must be utilized
@constraint(cd_modl, d_nd_start[i= cd_del_start, k=K],
                sum(x_ijk[i,j,k] for j=D if i != j) == 1)

# 03: end at cross dock
@constraint(cd_modl, d_nd_end[j= cd_del_end, k=K],
                sum(x_ijk[i,j,k] for i=D if i != j) == 1)

# 04: capacity of vehicle
@constraint(cd_modl, d_veh_cap[k=K],
                sum(d_i[i] * x_ijk[i,j,k] for i=D for j=vcat(D,cd_del_end) if i != j) <= veh_cap)

# 05: continuous routes/ consecutive movement
@constraint(cd_modl, d_move[h=D, k=K],
                    sum(x_ijk[i,h,k] for i=vcat(cd_del_start, D) if i != h)
                    - sum(x_ijk[h,j,k] for j=vcat(D,cd_del_end) if j != h) == 0)


# --------------------[time constraints]
# 01: arrival time of vehicle at pickup nodes
@constraint(cd_modl, arr_time_pick[i=vcat(cd_pick_start, P, cd_pick_end),j=vcat(cd_pick_start, P, cd_pick_end),k=K; i != j],
                    s_ik[j,k]
                    - s_ik[i,k]
                    - t_ij[i,j]
                    + M*(1-x_ijk[i,j,k])
                     >= 0)

# arrival time of vehicle at delivery nodes
@constraint(cd_modl, arr_time_del[i=vcat(cd_del_start, D, cd_del_end),j=vcat(cd_del_start, D, cd_del_end),k=K; i != j],
                 s_ik[j,k]
                 - s_ik[i,k]
                 - t_ij[i,j]
                 + M*(1-x_ijk[i,j,k])
                  >= 0)

# --------------------[consolidation decision]

# consolidation decisions
@constraint(cd_modl, consolidation[i=P,k=K],
                u_ik[i,k]
                - r_ik[i,k]
                - sum(x_ijk[i,j,k] for j = vcat(P,cd_pick_end) if i!=j)
                + sum(x_ijk[i+n,j,k] for j = vcat(D,cd_del_end)if i+n !=j)
                == 0)
# unloading and loading does not happen for same vehicles at cross dock
@constraint(cd_modl, load_unload[i=P,k=K],
                u_ik[i,k]
                + r_ik[i,k]
                <= 1)



# check (1 and 2) if the vehicles has to stop at cross-dock using decision from request
# this variable is used to calculate the unloading time of vehicle next
@constraint(cd_modl, unload_check1[k=K],
                (1/M)*sum(u_ik[i,k] for i=P) <= g_k[k])

@constraint(cd_modl, unload_check2[k=K],
                g_k[k] <= sum(u_ik[i,k] for i=P))


# Calculate Finish time of trucks for unloading
@constraint(cd_modl, unld_f_tm_veh[k=K],
                f_k[k]
                == s_ik[cd_pick_end[1],k]
                + 1*g_k[k]
                + 2*sum(d_i[i]*u_ik[i,k] for i=P)
                )

# unloading time of request i = finish time of veh k (f_k) if veh k unloads
# req i (u_ik) at the cross dock
@constraint(cd_modl, unld_req_tm[i=P, k=K],
            v_i[i]>= f_k[k] - M*(1-u_ik[i,k])
            )

# reloading in trucks (w_k) starts after the finish time (f_k) of unloading at cross docks
# of the same truck
@constraint(cd_modl, unld_reld[k=K],
            w_k[k]>= f_k[k])


@constraint(cd_modl, reld_start_tm[i=P, k=K],
            w_k[k]>= v_i[i] - M*(1-r_ik[i,k]))

# check (1 and 2) if the vehicles has to stop at cross-dock using decision from request
# this variable is used to calculate the unloading time of vehicle next
@constraint(cd_modl, reload_check1[k=K],
                (1/M)*sum(r_ik[i,k] for i=P) <= h_k[k])

@constraint(cd_modl, reload_check2[k=K],
                h_k[k] <= sum(u_ik[i,k] for i=P))


# Calculate start time of departure at cross dock of trucks
@constraint(cd_modl, leave_cd_time[k=K],
                s_ik[cd_del_start[1],k]
                == w_k[k]
                + 1*h_k[k]
                + 2*sum(d_i[i]*r_ik[i,k] for i=P)
                )

# ****************************************
# Objective
#*****************************************
@objective(cd_modl, Min,
sum(t_ij[i,j] * x_ijk[i,j,k] for i=N for j=N for k=K))


optimize!(cd_modl)


if termination_status(cd_modl) == MOI.OPTIMAL
    optimal_solution = value.(x_ijk)
    optimal_objective = objective_value(cd_modl)
elseif termination_status(cd_modl) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x_ijk)
    suboptimal_objective = objective_value(cd_modl)
else
    error("The model was not solved correctly.")
end
