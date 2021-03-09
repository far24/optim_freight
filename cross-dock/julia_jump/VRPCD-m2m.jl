# vehicle routing problem with cross-dock
# Wen et al. (2009)

#*****************************************
# Necessary Package
#*****************************************
using Random
using JuMP
using GLPK, Cbc
using LinearAlgebra
using CSV
using Distances

#*****************************************
# Network Property
#*****************************************
# suppliers are the pick-up nodes
P = [1,2]

# customers are the delivery nodes
C = [3,4]

# cross-docks
cd_pick_start = [C[end]+1]
cd_pick_end = [C[end]+2]
cd_del_start = [C[end]+3]
cd_del_end = [C[end]+4]

N = vcat(P,C, cd_pick_start, cd_pick_end, cd_del_start, cd_del_end)

# nodes in pickup Process
S = vcat(P, cd_pick_start, cd_pick_end)
D = vcat(C, cd_del_start, cd_del_end)

# vehicles
#K = [1,2]
K_S = [1,2]
K_D = [3,4]



# ****************************************
# parameters
#*****************************************
# number of pickup/delivery Nodes
n = length(P)

# Q = the vehicle capacity;
Q = 80

# location of nodes on a graph
#c = [9,5,4,8,13,11,8,8]

c = [10, 20, 65, 80, 55,55,55,55]
t_ij = pairwise(Euclidean(), c'; dims=2)
t_ij

# many-to-many relationship between supplier and customer
CV = [
1 1 0 0;
1 0 0 1;
0 1 1 0;
]

# quantity of products from suppliers (rows) to customers (cols)
prod = [
20 30 0 0;
15 0 0 25;
0 25 25 0;
]

# total product supplied by supplier i
p_i = sum(prod, dims = 2)

# total demand of customer i
d_i = sum(prod, dims = 1)

#  TIME RELATED
# service time pickup location
sp_i = 10
# service time delivery location
sd_i = 10
# duration of unloading Process at CD
u = 5
# duration of loading Process at CD
l = 5
# maximum route duration
T = 99999

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
# see Nikolopoulou et al. (2016)
# x_ijk: 1 if veh k travels from node i to node j
@variable(pick_mod, x_ijk[i= vcat(cd_pick_start,P,cd_pick_end) ,j= vcat(cd_pick_start,P,cd_pick_end),k=K_S], Bin)
@variable(del_mod, z_ijk[i= vcat(cd_del_start, D, cd_del_end) ,j= vcat(cd_del_start, D, cd_del_end),k=K_D], Bin)

# y_ik: 1 if node i is serviced by veh k
@variable(pick_mod, y_ik[i= vcat(cd_pick_start,P,cd_pick_end), k=K_S], Bin)
@variable(del_mod, v_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# tp_ik :: Time at which vehicle k leaves pickup node i (i [ S, k [ KS)
# td_hl :: Time at which vehicle l leaves delivery node h (h [ D, l [ KD)
@variable(pick_mod, atp_ik[i=vcat(cd_pick_start,P,cd_pick_end), k=K_S] >=0)
@variable(del_mod, atd_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D] >=0)

# RT_k :: Release time of pickup vehicle k at the cross-dock (k [ KS)
# DT_l :: Starting time of delivery vehicle l at the cross dock (l [ KD)
@variable(cd_modl, rt_k[k=K_S] >=0)
@variable(cd_modl, dt_k[k=K_D] >=0)


# ****************************************
# Constraints
#*****************************************

# --------------------[Pick-up Process]-----------------------------------------
# Each node is visited once by one vehicle
@constraint(pick_mod, p_one_veh[i= P],
                    sum(y_ik[i,k] for k=K_S) == 1)

@constraint(pick_mod, p_cd_one_veh[i= vcat(cd_pick_start,cd_pick_end)],
                    sum(y_ik[i,k] for k= K_S ) <= length(K_S))

@constraint(pick_mod, p_nd_entry[j= P, k= K_S],
                    sum(x_ijk[i,j,k] for i=vcat(cd_pick_start,P) if i != j)
                        == y_ik[j,k]
            )

@constraint(pick_mod, p_nd_exit[i= P, k=K_S],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P) if i != j)
                        == y_ik[i,k]
            )

# flow conservation
@constraint(pick_mod, p_mvmt[l= P, k=K_S],
                    sum(x_ijk[i,l,k] for i=vcat(cd_pick_start,P) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P) if l != j)
                    == 0
            )

# if vehicle is used its route starts from CD
@constraint(pick_mod, p_veh_start[i= cd_pick_start, k = K_S],
                    sum(x_ijk[i,j,k] for j=P) <= y_ik[i,k])


# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(pick_mod, p_veh_end[j= cd_pick_end, k = K_S],
                sum(x_ijk[i,j,k] for i=P) <= y_ik[j,k])


# Capacity Constraints
@constraint(pick_mod, p_veh_cap[k=K_S],
                    sum(p_i[i] * y_ik[i,k] for i=P) <= Q)

# -----------[time constraints]
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(pick_mod, p_veh_arr[i=vcat(cd_pick_start,P), j=vcat(P,cd_pick_end), k=K_S; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )

# Compute the departure time to zero if not visited
@constraint(pick_mod, p_init_dep_time[i=vcat(cd_pick_start, P, cd_pick_end), k=K_S],
                atp_ik[i,k] <= M*y_ik[i,k])


# --------------------[delivery Process]-----------------------------------------
# Each node is visited once by one vehicle
@constraint(del_mod, d_one_veh[i= D],
                    sum(v_ik[i,k] for k=K_D) == 1)

@constraint(del_mod, d_cd_one_veh[i= vcat(cd_del_start,cd_del_end), k=K_D],
                    sum(v_ik[i,k] ) <= length(K_D))

# only one vehicle enters each node
@constraint(del_mod, d_nd_entry[j= D, k=K_D],
                    sum(z_ijk[i,j,k] for i=vcat(cd_del_start,D) if i != j)
                        == v_ik[j,k]
            )

# only one vehicle leaves each node
@constraint(del_mod, d_nd_exit[i= D,k=K_D],
                    sum(z_ijk[i,j,k] for j=vcat(cd_del_end,D) if i != j)
                        == v_ik[i,k]
            )

# flow conservation
@constraint(del_mod, d_mvmt[l= D, k=K_D],
                    sum(z_ijk[i,l,k] for i=vcat(cd_del_start,D) if i != l)
                    - sum(z_ijk[l,j,k] for j=vcat(cd_del_end,D) if l != j)
                    == 0
            )

# if vehicle is used its route starts from CD
@constraint(del_mod, d_veh_start[i= cd_del_start, k = K_D],
                    sum(z_ijk[i,j,k] for j=D) <= v_ik[i,k])


# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(del_mod, d_veh_end[j= cd_del_end, k = K_D],
                sum(z_ijk[i,j,k] for i=D) <= v_ik[j,k])


# Capacity Constraints
@constraint(del_mod, d_veh_cap[k=K_D],
                    sum(d_i[i] * v_ik[i,k] for i=D) <= Q)

# -----------[time constraints]
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(del_mod, d_veh_arr[i=vcat(cd_del_start,D), j=vcat(D,cd_del_end), k=K_D; i !=j],
                    atd_ik[j,k]
                    >= atd_ik[i,k] + sp_i + t_ij[i,j]  - M*(1-z_ijk[i,j,k])
            )

# Compute the departure time to zero if not visited
@constraint(del_mod, d_init_dep_time[i=D, k=K_D],
                atd_ik[i,k] <= M*v_ik[i,k])


# -----------[At Cross Dock:: time constraints]
# release time of veh k (in pickup process)




# ****************************************
# Objective
#*****************************************
@objective(cd_modl, Min,
sum(t_ij[i,j] * x_ijk[i,j,k] for i=vcat(cd_pick_start,P) for j=vcat(P,cd_pick_end) for k=K_S)
+ sum(t_ij[i,j] * z_ijk[i,j,k] for i=vcat(cd_del_start,D) for j=vcat(D,cd_del_end) for k=K_D)
)


print("------------------------starting optimization--------------------------")
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




for k=K
    print("\nveh: ", k, "\t")
    for i= N
        for j=N
            if value.(x_ijk[i,j,k]) == 1
                print(i, "-->", j, "\t")
            end
        end
    end
end



for k=K
    print("\nveh: ", k)
    for i= N
        print("\tnode: ", i, "\t\t time: ", value.(s_ik[i,k]), "\n")
    end
end
