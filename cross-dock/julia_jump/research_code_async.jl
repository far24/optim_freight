#*****************************************
# Necessary Package
#*****************************************
using Random
using JuMP
using LinearAlgebra
using GLPK

using CSV

# vehicle routing problem with cross-dock
# Wen et al. (2009)

#*****************************************
# Network Property
#*****************************************
# suppliers are the pick-up nodes
P = [1,2,3]

# customers are the delivery nodes
D = [4,5,6]

# cross-docks
cross_pickup = [7]
cross_delivery = [8]

N = vcat(P,D, cross_pickup, cross_delivery)

# vehicles
K = [1, 2]

# ****************************************
# parameters
#*****************************************
# number of pickup/delivery Nodes
n = length(P)

d_i = [10,20,30, 30, 20, 10]
veh_cap = 40

# location of nodes on a graph
c = c = [9,5,4,8,13,11,8,8]
t_ij = pairwise(Euclidean(), c'; dims=2)



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
@variable(cd_modl, q_k[k=K], Bin)
@variable(cd_modl, qq_k[k=K], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# l_ik = the time at which vehicle k leaves node i (i ∈ N, k ∈ K); 0 otherwise;
# f_k = finish time of unloading operation of truck k at cross-dock (k ∈ K)
# g_i = finish time of unloading operation of request i at cross-dock (i ∈ P)
# s_k = start time of reloading operation of truck k at cross-dock (k ∈ K)
@variable(cd_modl, l_ik[i=N, k=K] >= 0)
@variable(cd_modl, f_k[k=K] >= 0)
@variable(cd_modl, s_k[k=K] >= 0)
@variable(cd_modl, g_i[i=P] >= 0)


# ****************************************
# Constraints
#*****************************************

# --------------------[Pick-up Process]
# 01: only one vehicle can arrive at a pickup node
@constraint(cd_modl, p_nd_entry[i= P],
                    sum(x_ijk[i,j,k] for j=P for k=K if i != j) == 1)


# 02: start from cross dock // all the trucks must be utilized
@constraint(cd_modl, p_nd_start[i= cross_pickup, k=K],
                    sum(x_ijk[i,j,k] for j=P if i != j) == 1)

# 03: end at cross dock
@constraint(cd_modl, p_nd_end[j= cross_pickup, k=K],
                sum(x_ijk[i,j,k] for i=P if i != j) == 1)

# 04: capacity of vehicle
@constraint(cd_modl, p_veh_cap[k=K],
                sum(d_i[i] * x_ijk[i,j,k] for i=P for j=P if i != j) <= veh_cap)

# 05: continuous routes/ consecutive movement
@constraint(cd_modl, p_move[h=P, k=K],
                    sum(x_ijk[i,h,k] for i=P if i != h)
                    - sum(x_ijk[h,j,k] for j=P if j != h) == 0)


# --------------------[delivery Process]
# 01: only one vehicle can arrive at a pickup node
@constraint(cd_modl, d_nd_entry[i=D],
                    sum(x_ijk[i,j,k] for j=D for k=K if i != j) == 1)

# 02: start from cross dock // all the trucks must be utilized
@constraint(cd_modl, d_nd_start[i= cross_delivery, k=K],
                sum(x_ijk[i,j,k] for j=D if i != j) == 1)

# 03: end at cross dock
@constraint(cd_modl, d_nd_end[j= cross_delivery, k=K],
                sum(x_ijk[i,j,k] for i=D if i != j) == 1)

# 04: capacity of vehicle
@constraint(cd_modl, d_veh_cap[k=K],
                sum(d_i[i] * x_ijk[i,j,k] for i=D for j=D if i != j) <= veh_cap)

# 05: continuous routes/ consecutive movement
@constraint(cd_modl, d_move[h=D, k=K],
                    sum(x_ijk[i,h,k] for i=D if i != h)
                    - sum(x_ijk[h,j,k] for j=D if j != h) == 0)


# --------------------[time constraints]
# 01: arrival time of vehicle at node if
@constraint(cd_modl, arr_time[i=N,j=N,k=K; i != j],
                    l_ik[j,k]
                    - l_ik[i,k]
                    - t_ij[i,j]
                    + M*(1-x_ijk[i,j,k])
                     >= 0)

@constraint(cd_modl, arr_time_init[i=N,j=N,k=K; i != j],
                    l_ik[j,k] <= M*x_ijk[i,j,k])


# --------------------[consolidation decision]

# consolidation decisions
@constraint(cd_modl, consolidation[i=P,k=K],
                u_ik[i,k]
                - r_ik[i,k]
                - sum(x_ijk[i,j,k] for j = vcat(P,cross_pickup) if i!=j)
                + sum(x_ijk[i+n,j,k] for j = vcat(D,cross_delivery)if i+n !=j)
                == 0)
# unloading and loading does not happen for same vehicles
@constraint(cd_modl, load_unload[i=P,k=K],
                u_ik[i,k]
                + r_ik[i,k]
                <= 1)

# check for unloading process of truck at cross-docks
@constraint(cd_modl, unload_check1[k=K],
                (1/M)* sum(u_ik[i,k] for i=P) - q_k[k] <= 0)

@constraint(cd_modl, unload_check2[k=K],
                q_k[k] - sum(u_ik[i,k] for i=P) <= 0)

# Calculate Finish time of trucks for unloading
@constraint(cd_modl, unld_f_tm_veh[k=K],
                f_k[k]
                == l_ik[cross_pickup[1],k]
                + 2*q_k[k]
                + 0.05*sum(d_i[i]*u_ik[i,k] for i=P)
                )

# calculate finish time of unloading demand/request "i"
@constraint(cd_modl, unld_f_tm_req[i=P,k=K],
            g_i[i] >= f_k[k] - M*(1-u_ik[i,k]))

# calculate start time of reloading in the trucks
@constraint(cd_modl, strt_reld_tm_veh[i=P,k=K],
            s_k[k] >= g_i[i] - M*(1-r_ik[i,k]))

# reloading starts after unloading
@constraint(cd_modl, unld_reld[k=K],
            s_k[k]>= f_k[k])


# check 01 for reloading process of truck at cross-docks
@constraint(cd_modl, reload_check1[k=K],
                (1/M)* sum(r_ik[i,k] for i=P) - qq_k[k] <= 0)


# check 02 for reloading process of truck at cross-docks
@constraint(cd_modl, reload_check2[k=K],
                qq_k[k] - sum(r_ik[i,k] for i=P) <= 0)

# Calculate Finish time of trucks for reloading
@constraint(cd_modl, reld_s_tm_veh[k=K],
                l_ik[cross_delivery[1],k]
                ==  s_k[k]
                + 1*qq_k[k]
                + 0.05*sum(d_i[i]*r_ik[i,k] for i=P)
                )



# ****************************************
# Objective
#*****************************************
