P_un = [1]
K_S_new = [5]
K_S_exst = [2]
K_S = vcat(K_S_new, K_S_exst)
O_K_S_exst = [2]
current_time = 30
sch_arr = [0,100]


# ****************************************
# Create a JuMP model
# ****************************************
cd_reroute_mod = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************

# BINARY
# ^^^^^^^^^
# see Nikolopoulou et al. (2016)
# x_ijk: 1 if veh k travels from node i to node j
@variable(cd_reroute_mod, x_ijk[i= vcat(cd_pick_start,P,cd_pick_end) ,j= vcat(cd_pick_start,P,cd_pick_end),k=K_S], Bin)
@variable(cd_reroute_mod, z_ijk[i= vcat(cd_del_start, D, cd_del_end) ,j= vcat(cd_del_start, D, cd_del_end),k=K_D], Bin)

# y_ik: 1 if node i is serviced by veh k
@variable(cd_reroute_mod, y_ik[i= vcat(cd_pick_start,P,cd_pick_end), k=K_S], Bin)
@variable(cd_reroute_mod, v_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# atp_ik :: Time at which vehicle k leaves pickup node i (i [ S, k [ KS)
# atd_hl :: Time at which vehicle l leaves delivery node h (h [ D, l [ KD)
# ld_ik = tardiness (tardiness/lateness = upper TW (b_i) - arrival time (s_ik)) at node i by vehicle k ==> this would be a non-linear constraint
# ld_i = tardiness at node i
@variable(cd_reroute_mod, atp_ik[i=vcat(cd_pick_start,P,cd_pick_end), k=K_S] >=0)
@variable(cd_reroute_mod, atd_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D] >=0)
@variable(cd_reroute_mod, ld_ik[i=vcat(cd_del_start,D,cd_del_end), k = K_D] >= 0)
@variable(cd_reroute_mod, l_cdk[i=cd_pick_end, k = K_S] >= 0)


# RT_k :: Release time of pickup vehicle k at the cross-dock (k [ KS)
# DT_l :: Starting time of delivery vehicle l at the cross dock (l [ KD)
@variable(cd_reroute_mod, rt_k[k=K_S] >=0)
@variable(cd_reroute_mod, dt_k[k=K_D] >=0)


# ****************************************
# Constraints
#*****************************************

# PICKUP PROCESS -------------------------
#*****************************************
# --------------------[new vehicles starting from depot]-----------------------------------------
# [new vehicle:: K_S_new]
# entry to nodes
@constraint(cd_reroute_mod, p_re_nd_entry[j= P_un, k= K_S_new],
                    sum(x_ijk[i,j,k] for i=vcat(cd_pick_start,P_un) if i != j)
                        == y_ik[j,k]
)

# [new vehicle:: K_S_new]
# exit from node
@constraint(cd_reroute_mod, p_re_nd_exit[i= P_un, k=K_S_new],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P_un) if i != j)
                        == y_ik[i,k]
)

# [new vehicle:: K_S_new]::
# flow conservation
@constraint(cd_reroute_mod, p_re_mvmt[l= P_un, k=K_S_new],
                    sum(x_ijk[i,l,k] for i=vcat(cd_pick_start,P_un) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P_un) if l != j)
                    == 0
)

# [new vehicle:: K_S_new]::
# if vehicle is used its route starts from CD
@constraint(cd_reroute_mod, p_re_newveh_start[i= cd_pick_start, k = K_S_new],
                    sum(x_ijk[i,j,k] for j=P_un) <= y_ik[i,k]
)

# [new vehicle:: K_S_new]::
# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(cd_reroute_mod, p_re_newveh_end[j= cd_pick_end, k = K_S_new],
                sum(x_ijk[i,j,k] for i=P_un) <= y_ik[j,k]
)

# [new vehicle:: K_S_new]::
# Capacity Constraints
@constraint(cd_reroute_mod, p_re_newveh_cap[k=K_S_new],
                    sum(p_i[i] * y_ik[i,k] for i=P_un) <= Q
)


# [pickup process:: time constraints] ^^^^^^^^^^^^^^^^^^^^^^
# [new vehicle:: K_S_new]:: ****current_time added******************
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(cd_reroute_mod, p_re_newveh_arr[i=vcat(cd_pick_start,P_un), j=vcat(P_un,cd_pick_end), k=K_S_new; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )

# [new vehicle:: K_S_new]
# Compute the departure time to zero if not visited
@constraint(cd_reroute_mod, p_re_init_dep_time[i=vcat(cd_pick_start, P, cd_pick_end), k=K_S],
            atp_ik[i,k] <= M*y_ik[i,k])



# --------------------[update information of existing vehicle]------------------
# routes already taken
@constraint(cd_reroute_mod, p_update_route[i = [5], j= [2] ,k =K_S_exst],
            x_ijk[i,j,k] ==1
)

# nodes already visited by the existing vehicles
@constraint(cd_reroute_mod, p_update_node_start[i= cd_pick_start ,k =K_S_exst],
            y_ik[i,k] ==1
)

@constraint(cd_reroute_mod, p_update_node[i= [2] ,k =K_S_exst],
            y_ik[i,k] ==1
)

@constraint(cd_reroute_mod, p_update_node_end[i= cd_pick_end ,k =K_S_exst],
            y_ik[i,k] ==1
)

# arrival time of the NEW vehicle to cross-dock
@constraint(cd_reroute_mod, p_update_cdarr[i= cd_pick_start ,k =K_S_new],
            atp_ik[i,k] == current_time
)

# arrival time of the existing vehicle to the latest node
@constraint(cd_reroute_mod, p_update_arr[i= [2] ,k =K_S_exst],
            atp_ik[i,k] == 35
)


# --------------------[existing vehicle starting from latest arrival node]------
# [existing vehicle:: K_S_exist]
# entry to the node
@constraint(cd_reroute_mod, p_re_exstveh_entry[j= P_un, k= K_S_exst],
                    sum(x_ijk[i,j,k] for i=vcat([2],P_un) if i != j)
                        == y_ik[j,k]
)

# [existing vehicle:: K_S_exist]
# exit to the node
@constraint(cd_reroute_mod, p_re_exstveh_exit[i= P_un, k=K_S_exst],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P_un) if i != j)
                        == y_ik[i,k]
)

# [existing vehicle:: K_S_exist]::
# flow conservation
@constraint(cd_reroute_mod, p_re_exstveh_mvmt[l= P_un, k=K_S_exst],
                    sum(x_ijk[i,l,k] for i=vcat([2],P_un) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P_un) if l != j)
                    == 0
)

# [existing vehicle:: K_S_exist]::
# if vehicle is used its route starts from starting location of the vehicle
@constraint(cd_reroute_mod, p_re_exstveh_start[i= [2], k = K_S_exst],
                    sum(x_ijk[i,j,k] for j=P_un) <= y_ik[i,k]
)

# [existing vehicle:: K_S_exist]::
# Vehicle already en-route, so the vehicles route SHOULD end at cross-dock
@constraint(cd_reroute_mod, p_re_exstveh_end[j= cd_pick_end, k = K_S_exst],
                sum(x_ijk[i,j,k] for i=vcat([2],P_un)) == y_ik[j,k]
)

# [existing vehicle:: K_S_exist]::
# Capacity Constraints
@constraint(cd_reroute_mod, p_exstveh_cap[k=K_S_exst],
                    sum(p_i[i] * y_ik[i,k] for i=P_un) <= Q - sum(p_i[i] * y_ik[i,k] for i =[2])
)

# -----[pickup process:: time constraints]
# [existing vehicle:: K_S_exist]::
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(cd_reroute_mod, p_re_exstveh_arr[i=vcat([2],P_un), j=vcat(P_un,cd_pick_end), k=K_S_exst; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )
# --------------------[dependency between new and existing vehicle]-------------
# Each node is visited once by one vehicle, either new or existing
@constraint(cd_reroute_mod, p_re_one_veh[i= P_un],
                    sum(y_ik[i,k] for k=K_S) == 1)


# DELIVERY PROCESS -----------------------
#*****************************************
# Each node is visited once by one vehicle
@constraint(cd_reroute_mod, d_one_veh[i= D],
                    sum(v_ik[i,k] for k=K_D) == 1)

@constraint(cd_reroute_mod, d_cd_one_veh[i= vcat(cd_del_start,cd_del_end), k=K_D],
                    sum(v_ik[i,k] ) <= length(K_D))

# only one vehicle enters each node
@constraint(cd_reroute_mod, d_nd_entry[j= D, k=K_D],
                    sum(z_ijk[i,j,k] for i=vcat(cd_del_start,D) if i != j)
                        == v_ik[j,k]
            )

# only one vehicle leaves each node
@constraint(cd_reroute_mod, d_nd_exit[i= D,k=K_D],
                    sum(z_ijk[i,j,k] for j=vcat(cd_del_end,D) if i != j)
                        == v_ik[i,k]
            )

# flow conservation
@constraint(cd_reroute_mod, d_mvmt[l= D, k=K_D],
                    sum(z_ijk[i,l,k] for i=vcat(cd_del_start,D) if i != l)
                    - sum(z_ijk[l,j,k] for j=vcat(cd_del_end,D) if l != j)
                    == 0
            )

# if vehicle is used its route starts from CD
@constraint(cd_reroute_mod, d_veh_start[i= cd_del_start, k = K_D],
                    sum(z_ijk[i,j,k] for j=D) <= v_ik[i,k])


# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(cd_reroute_mod, d_veh_end[j= cd_del_end, k = K_D],
                sum(z_ijk[i,j,k] for i=D) <= v_ik[j,k])


# Capacity Constraints
@constraint(cd_reroute_mod, d_veh_cap[k=K_D],
                    sum(d_i[i] * v_ik[i,k] for i=D) <= Q)


# AT CROSS DOCK -----------------------
#**************************************
# -----------[At Cross Dock:: time constraints]
# release time of veh k (in pickup process)
@constraint(cd_reroute_mod, release_time[i=P, j= cd_pick_end, k=K_S],
                rt_k[k] >= atp_ik[i,k] + t_ij[i,j] + u - M*(1-x_ijk[i,j,k]) )


# delivery vehicle (l) cannot start unless the pickup vehicle (k)
# is unloaded (greater than release time of pickup vehicle)
# delivery route (h) starts if pickup vehicle brings the product from (i)
@constraint(cd_reroute_mod, dept_time[k=K_S,l=K_D, i=P,h=D],
                dt_k[l]*CV[i,h] >= rt_k[k] - M*(2 - y_ik[i,k] - CV[i,h]*v_ik[h,l]))

# calculate the tardiness of existing vehicle at the cross dock of pickup process
@constraint(cd_reroute_mod, tard_cd[j=cd_pick_end, k=K_S_exst],
    l_cdk[j,k]
    >= atp_ik[j,k] - sch_arr[k]
    - M*(1-y_ik[j,k])
    )


# DELIVERY PROCESS -----------------------
#*****************************************
# -----------[delivery process:: time constraints]
# set the departure time of vehicles from cross dock to be dt_k + loading time
@constraint(cd_reroute_mod, d_depCD[i= cd_del_start, j=D, k=K_D],
                    atd_ik[j,k]
                    >= dt_k[k] + l + t_ij[i,j] - M*(1-z_ijk[i,j,k])
            )

# time a delivery vehicle k arrives at delivery node j after leaving the cross dock
@constraint(cd_reroute_mod, d_veh_arr[i=D, j=D, k=K_D; i !=j],
                    atd_ik[j,k]
                    >= atd_ik[i,k] + sp_i + t_ij[i,j]  - M*(1-z_ijk[i,j,k])
            )

@constraint(cd_reroute_mod, d_veh_CDarr[i=D, j=cd_del_end, k=K_D; i !=j],
                atd_ik[j,k]
                >= atd_ik[i,k] + sp_i + t_ij[i,j]  - M*(1-z_ijk[i,j,k])
        )


@constraint(cd_reroute_mod, d_init_dep_time[i=vcat(cd_del_start, D, cd_del_end), k=K_D],
                atd_ik[i,k] <= M*v_ik[i,k])

# calculate the tardiness of vehicle at the customer location of delivery Process
@constraint(cd_reroute_mod, tard_del[j=vcat(cd_del_start, D, cd_del_end), k=K_D],
    ld_ik[j,k]
    >= atd_ik[j,k] - b_i[j]
    - M*(1-v_ik[j,k])
    )



# ****************************************
# Objective
#*****************************************
transp_cost = 200
@objective(cd_reroute_mod, Min,
sum(transp_cost* t_ij[i,j] * x_ijk[i,j,k] for i=vcat(cd_pick_start,P) for j=vcat(P,cd_pick_end) for k=K_S)
+ sum(transp_cost * t_ij[i,j] * z_ijk[i,j,k] for i=vcat(cd_del_start,D) for j=vcat(D,cd_del_end) for k=K_D)
+ sum(5*transp_cost * ld_ik[i,k] for i=vcat(cd_del_start,D,cd_del_end) for k=K_D)
+ sum(10*transp_cost * l_cdk[i,k] for i=cd_pick_end for k=K_S)
)


print("------------------------starting optimization--------------------------")
optimize!(cd_reroute_mod)


if termination_status(cd_reroute_mod) == MOI.OPTIMAL
    optimal_solution = value.(x_ijk)
    optimal_objective = objective_value(cd_reroute_mod)
    print("Objective value: ", optimal_objective)
elseif termination_status(cd_reroute_mod) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x_ijk)
    suboptimal_objective = objective_value(cd_reroute_mod)
else
    error("The model was not solved correctly.")
end


# OUTPUT from the model---------------------------------------------------------


using DataFrames

p_re_node_visit = DataFrame(node =[], veh =[], visit = [])
for k in keys(y_ik)
    key = Tuple(k[:])
    node = key[1]
    veh = key[2]
    visit = JuMP.value.(y_ik[key...])
    print("Node: ", node, "\t", "veh\t", veh, "visit ", visit, "\n")
    push!(p_re_node_visit, vcat(node, veh, visit))
end

d_re_node_visit = DataFrame(node =[], veh =[], visit = [])
for k in keys(v_ik)
    key = Tuple(k[:])
    node = key[1]
    veh = key[2]
    visit = JuMP.value.(v_ik[key...])
    print("Node: ", node, "\t", "veh\t", veh, "visit ", visit, "\n")
    push!(d_re_node_visit, vcat(node, veh, visit))
end

re_node_visit = vcat(p_re_node_visit, d_re_node_visit)

p_re_arr_time = DataFrame(node =[], veh =[], arr = [])
for k in keys(atp_ik)
    key = Tuple(k[:])
    node = key[1]
    veh = key[2]
    re_arr_time = JuMP.value.(atp_ik[key...])
    print("Node: ", node, "\t", "veh\t", veh, "\tarr ", re_arr_time, "\n")
    push!(p_re_arr_time, vcat(node, veh, re_arr_time))
end

d_re_arr_time = DataFrame(node =[], veh =[], arr = [])
for k in keys(atd_ik)
    key = Tuple(k[:])
    node = key[1]
    veh = key[2]
    arr_time = JuMP.value.(atd_ik[key...])
    print("Node: ", node, "\t", "veh\t", veh, "\tarr ", arr_time, "\n")
    push!(d_re_arr_time, vcat(node, veh, arr_time))
end

re_arr_time = vcat(p_re_arr_time, d_re_arr_time)

re_node_output = innerjoin(re_node_visit, re_arr_time, on = [:node,:veh])
