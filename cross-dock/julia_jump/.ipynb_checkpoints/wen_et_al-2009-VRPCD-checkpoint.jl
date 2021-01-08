# vehicle routing problem with cross-dock
# Wen et al. (2009)

#*****************************************
# Necessary Package
#*****************************************
using Random
using JuMP
using LinearAlgebra
using GLPK

#*****************************************
# Network Property
#*****************************************
# suppliers are the pick-up nodes
supplier = ["p1", "p2"]
# customers are the delivery nodes
customer = ["d1", "d2"]
# request
# cross-docks
cross_pickup = ["o1", "o2"]
cross_delivery = ["o3", "o4"]
# vehicles
veh = ["v1", "v2"]

# ****************************************
# parameters
#*****************************************
P = length(supplier)
D = length(customer)
OP = length(cross_pickup)
OD = length(cross_delivery)
O = OP + OD
N = P + D + O
K = length(veh)

n = P

Random.seed!(314)
# ci j = the travel time between node i and node j ((i, j ) ∈ E);
c_ij = rand(10:60, N,N)
# assign 0 to the diagonals
c_ij[diagind(c_ij)] .= 0

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
M= 100000
# ****************************************
# Create a JuMP model
# ****************************************
cd_modl = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************
# x_ijk = 1 if vehicle k travels from node i to node j ((i, j ) ∈ E, k ∈ K)
#        0 otherwise;
# u_ik = 1 if vehicle k unloads request i at the cross-dock (i ∈ P, k ∈ K)
# 0 otherwise;
# r_ik =1 if vehicle k reloads request i at the cross-dock (i ∈ P, k ∈ K)
# 0 otherwise;
# gk = 1 if vehicle k has to unload at the cross-dock (k ∈ K)
# 0 otherwise;
# h_k = 1 if vehicle k has to reload at the cross-dock (k ∈ K)
# 0 otherwise;
# s_ik = the time at which vehicle k leaves node i (i ∈ N, k ∈ K);
@variable(cd_modl, x_ijk[i=1:N,j=1:N,1:K], Bin)
@variable(cd_modl, r_ik[1:P, 1:K], Bin)
@variable(cd_modl, u_ik[1:P, 1:K], Bin)
@variable(cd_modl, h_k[1:K], Bin)
@variable(cd_modl, g_k[1:K], Bin)
@variable(cd_modl, s_ik[1:N, 1:K], Bin)

#*****************************************
# Constraints
#*****************************************
#01: Each node visited once by the vehicle
@constraint(cd_modl, node_visit[i=1:P+D],
                    sum(x_ijk[i,j,k] for j=1:N for k=1:K if i != j) == 1)

#02: load on the pickup route does not exit the truck capacity
@constraint(cd_modl, truck_cap_pickup[k=1:K],
                    sum(d_i[i]*x_ijk[i,j,k] for j=vcat(1:P,[5,6]) for i=1:P if i != j) <= Q)

#03: load on the pickup route does not exit the truck capacity
@constraint(cd_modl, truck_cap_delivery[k=1:K],
                    sum(d_i[i]*x_ijk[i,j,k] for i=P+1:P+D for j=vcat(P+1:P+D,[7,8]) if i != j) <= Q)

#04: each vehicle’s pickup route must depart from o1 and delivery route must leave from o3
# pickup side
@constraint(cd_modl, node_exit_P[h=[5],k=1:K],
                    sum(x_ijk[h,j,k] for j=1:P if h != j) == 1)
# delivery side
@constraint(cd_modl, node_exit_D[h=[7],k=1:K],
                    sum(x_ijk[h,j,k] for j=P+1:P+D if h != j) == 1)

#05: flow conservation constraints.
@constraint(cd_modl, flow_conserv[k=1:K, h=1:P+D],
                    sum(x_ijk[i,h,k] for i=1:N)
                    - sum(x_ijk[h,j,k] for j=1:N) == 0)

#06: each vehicle to return to o2 on its pickup route and return to o4 on its delivery route
# Pickup side
@constraint(cd_modl, node_return_P[h=6,k=1:K],
                    sum(x_ijk[j,h,k] for j=1:P if j != h) == 1)

# Delivery side
@constraint(cd_modl, node_return_D[h=8,k=1:K],
                    sum(x_ijk[j,h,k] for j=P+1:P+D if j != h) == 1)

#07: travelling time between two nodes if they are visited consecutively by the same vehicle
@constraint(cd_modl, tt[i=1:N, j=1:N, k=1:K; i != j],
                    s_ik[j,k]
                    >= s_ik[i,k]
                        + c_ij[i,j]
                        - M*(1-x_ijk[i,j,k]) )

#08: each node is visited within its time window and the whole operation is completed within the time horizon.
@constraint(cd_modl, time_window[i=1:N,k=1:K],
                    a_i[i] <= s_ik[i,k] <= b_i[i])

#09: Link between pickup and delivery part (??)
@constraint(cd_modl, link_pd1[i=1:P, k=1:K],
                    u_ik[i,k] - r_ik[i,k]
                    == sum(x_ijk[i,j,k] for j=vcat(1:P,6))
                        - sum(x_ijk[i,j,k] for j=vcat(P+1:P+D,8))
            )

#10: Link between pickup and delivery part (??)
@constraint(cd_modl, link_pd_cd2[i=1:P, k=1:K],
                    u_ik[i,k] + r_ik[i,k]
                    <= 1
            )

#11:
@constraint(cd_modl, unld_cd1[k=1:K],
                    sum((u_ik[i,k] / M) for i=1:P)
                    <= g_k[k]
            )

@constraint(cd_modl, unld_cd2[k=1:K],
                    g_k[k]
                    <= sum(u_ik[i,k] for i=1:P)
            )

#12

#13

#14

#15

#16

#17
