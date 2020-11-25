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
supplier = Set{String}(["p1", "p2"])
# customers are the delivery nodes
customer = Set{String}(["d1", "d2"])
# request
# cross-docks
cross_pickup = Set{String}(["o1", "o2"])
cross_delivery = Set{String}(["o3", "o4"])
# vehicles
veh = Set{String}(["v1", "v2"])

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

Random.seed!(314)
# ci j = the travel time between node i and node j ((i, j ) ∈ E);
c_ij = rand(10:60, N,N)
# assign 0 to the diagonals
c_ij[diagind(c_ij)] .= 0

# [ai , bi] = the time window for node i (i ∈ N);
a_i = zeros(N)
b_i = rand(5:10, N)
# di = the amount of demand of request i (i ∈ P);
di = rand(5:10, P)
# Q = the vehicle capacity;
Q = 50
# A = the fixed time for unloading and reloading at the cross-dock;
A = rand(5:10, O)
# B = the time for unloading and reloading a pallet.
A = rand(5:10, O)

# ****************************************
# Create a JuMP model
# ****************************************
cd_modl = Model(GLPK.Optimizer)

# ****************************************
# decision variables
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

# Each node visited once by the vehicle
@constraint(cd_modl, node_visit[j=1:P+D],sum(x_ijk[i,j,k] for i=1:N for k=1:K) == 1)
# load on the pickup route does not exit the truck capacity
@constraint(cd_modl, truck_cap_pickup[k=1:K],sum(di[i]*x_ijk[i,j,k] for i=1:P for j=1:N) <= Q)
# load on the pickup route does not exit the truck capacity
@constraint(cd_modl, truck_cap_delivery[k=1:K],sum(di[i]*x_ijk[i,j,k] for i=1:D for j=1:N) <= Q)
#each vehicle’s pickup route must depart from o1 and delivery route must leave from o3
@constraint(cd_modl, node_exit[h=[5,7],k=1:K],sum(x_ijk[h,j,k] for j=1:N if h != j) == 1)
# flow conservation constraints.
@constraint(cd_modl, flow_conserv[k=1:K, h=1:P+D],sum(x_ijk[i,h,k] for i=1:N) - sum(x_ijk[h,j,k] for j=1:N) == 1)
# each vehicle to return to o2 on its pickup route and return to o4 on its delivery route
@constraint(cd_modl, node_return[h=[6,8],k=1:K],sum(x_ijk[j,h,k] for j=1:N if h != j) == 1)
