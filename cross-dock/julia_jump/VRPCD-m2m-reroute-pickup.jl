P= [1,2]
P_un = [1]
K_S_new = [5]
K_S_exst = [2]
K_S = vcat(K_S_new, K_S_exst)
O_K_S_exst = [2]
current_time = 30






# ****************************************
# Create a JuMP_un model
# ****************************************
pick_reroute_mod = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************

# BINARY
# ^^^^^^^^^
# see Nikolopoulou et al. (2016)
# x_ijk: 1 if veh k travels from node i to node j
@variable(pick_reroute_mod, x_ijk[i= vcat(cd_pick_start,P,cd_pick_end) ,j= vcat(cd_pick_start,P,cd_pick_end),k=K_S], Bin)

# y_ik: 1 if node i is serviced by veh k
@variable(pick_reroute_mod, y_ik[i= vcat(cd_pick_start,P,cd_pick_end), k=K_S], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# atp_ik :: Time at which vehicle k leaves pickup node i (i [ S, k [ KS)
@variable(pick_reroute_mod, atp_ik[i=vcat(cd_pick_start,P,cd_pick_end), k=K_S] >=0)

# RT_k :: Release time of pickup vehicle k at the cross-dock (k [ KS)
# DT_l :: Starting time of delivery vehicle l at the cross dock (l [ KD)
@variable(pick_reroute_mod, rt_k[k=K_S] >=0)

# ****************************************
# Constraints
#*****************************************

# --------------------[P_unick-up P_unrocess]-----------------------------------------

# --------------------[new vehicles starting from depot]-----------------------------------------
# [new vehicle:: K_S_new]
@constraint(pick_reroute_mod, p_re_nd_entry[j= P_un, k= K_S_new],
                    sum(x_ijk[i,j,k] for i=vcat(cd_pick_start,P_un) if i != j)
                        == y_ik[j,k]
)

# [new vehicle:: K_S_new]
@constraint(pick_reroute_mod, p_re_nd_exit[i= P_un, k=K_S_new],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P_un) if i != j)
                        == y_ik[i,k]
)

# [new vehicle:: K_S_new]:: flow conservation
@constraint(pick_reroute_mod, p_re_mvmt[l= P_un, k=K_S_new],
                    sum(x_ijk[i,l,k] for i=vcat(cd_pick_start,P_un) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P_un) if l != j)
                    == 0
)

# [new vehicle:: K_S_new]::
# if vehicle is used its route starts from CD
@constraint(pick_reroute_mod, p_re_newveh_start[i= cd_pick_start, k = K_S_new],
                    sum(x_ijk[i,j,k] for j=P_un) <= y_ik[i,k]
)

# [new vehicle:: K_S_new]::
# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(pick_reroute_mod, p_re_newveh_end[j= cd_pick_end, k = K_S_new],
                sum(x_ijk[i,j,k] for i=P_un) <= y_ik[j,k]
)

# [new vehicle:: K_S_new]::
# Capacity Constraints
@constraint(pick_reroute_mod, p_re_newveh_cap[k=K_S_new],
                    sum(p_i[i] * y_ik[i,k] for i=P_un) <= Q
)

# [new vehicle:: K_S_new]:: ****current_time added******************
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(pick_reroute_mod, p_re_newveh_arr[i=vcat(cd_pick_start,P_un), j=vcat(P_un,cd_pick_end), k=K_S_new; i !=j],
                    atp_ik[j,k]
                    >= current_time + atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )
# --------------------[update information of existing vehicle]------------------
# routes already taken
@constraint(pick_reroute_mod, p_update_route[i = [5], j= [2] ,k =K_S_exst],
            x_ijk[i,j,k] ==1
)

# nodes already visited by the existing vehicles
@constraint(pick_reroute_mod, p_update_node_start[i= cd_pick_start ,k =K_S_exst],
            y_ik[i,k] ==1
)

@constraint(pick_reroute_mod, p_update_node[i= [2] ,k =K_S_exst],
            y_ik[i,k] ==1
)

@constraint(pick_reroute_mod, p_update_node_end[i= cd_pick_end ,k =K_S_exst],
            y_ik[i,k] ==1
)

# arrival time of the existing vehicle to the node
@constraint(pick_reroute_mod, p_update_arr[i= [2] ,k =K_S_exst],
            atp_ik[i,k] == 35
)

# --------------------[existing vehicle starting from latest arrival node]------
# [existing vehicle:: K_S_exist]
@constraint(pick_reroute_mod, p_re_exstveh_entry[j= P_un, k= K_S_exst],
                    sum(x_ijk[i,j,k] for i=vcat([2],P_un) if i != j)
                        == y_ik[j,k]
)

# [existing vehicle:: K_S_exist]
@constraint(pick_reroute_mod, p_re_exstveh_exit[i= P_un, k=K_S_exst],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P_un) if i != j)
                        == y_ik[i,k]
)

# [existing vehicle:: K_S_exist]:: flow conservation
@constraint(pick_reroute_mod, p_re_exstveh_mvmt[l= P_un, k=K_S_exst],
                    sum(x_ijk[i,l,k] for i=vcat([2],P_un) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P_un) if l != j)
                    == 0
)

# [existing vehicle:: K_S_exist]::
# if vehicle is used its route starts from starting location of the vehicle
@constraint(pick_reroute_mod, p_re_exstveh_start[i= [2], k = K_S_exst],
                    sum(x_ijk[i,j,k] for j=P_un) <= y_ik[i,k]
)

# [existing vehicle:: K_S_exist]::
# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(pick_reroute_mod, p_re_exstveh_end[j= cd_pick_end, k = K_S_exst],
                sum(x_ijk[i,j,k] for i=vcat([2],P_un)) == y_ik[j,k]
)

# [existing vehicle:: K_S_exist]::
# Capacity Constraints
@constraint(pick_reroute_mod, p_exstveh_cap[k=K_S_exst],
                    sum(p_i[i] * y_ik[i,k] for i=P_un) <= Q - sum(p_i[i] * y_ik[i,k] for i =[2])
)

# [existing vehicle:: K_S_exist]::
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(pick_reroute_mod, p_re_exstveh_arr[i=vcat([2],P_un), j=vcat(P_un,cd_pick_end), k=K_S_exst; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )
# --------------------[dependency between new and existing vehicle]-------------
# Each node is visited once by one vehicle, either new or existing
@constraint(pick_reroute_mod, p_re_one_veh[i= P_un],
                    sum(y_ik[i,k] for k=K_S) == 1)


# ****************************************
# Objective
#*****************************************
P=[1,2]
@objective(pick_reroute_mod, Min,
sum(t_ij[i,j] * x_ijk[i,j,k] for i=vcat(cd_pick_start,P) for j=vcat(P,cd_pick_end) for k=K_S)
#    +sum(t_ij[h,f] * z_hfl[h,f,l] for h=D for f=D for l=K_D)
)


print("------------------------starting optimization--------------------------")
optimize!(pick_reroute_mod)

if termination_status(pick_reroute_mod) == MOI.OPTIMAL
    optimal_solution = value.(x_ijk)
    optimal_objective = objective_value(pick_reroute_mod)
    print("objective value: ",optimal_objective)
elseif termination_status(pick_reroute_mod) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x_ijk)
    suboptimal_objective = objective_value(pick_reroute_mod)
else
    error("The model was not solved correctly.")
end


for k=K_S
    print("\nveh: ", k, "\t")
    for i=vcat(cd_pick_start,P)
        for j=vcat(P,cd_pick_end)
            if value.(x_ijk[i,j,k]) == 1
                print(i, "-->", j, "\t")
            end
        end
    end
end
