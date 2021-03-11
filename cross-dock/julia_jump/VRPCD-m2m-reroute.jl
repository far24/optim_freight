P_un = [1]
K_S_new = [5]
K_S_exst = [2]
K_S = vcat(K_S_new, K_S_exst)
O_K_S_exst = [2]



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

# RT_k :: Release time of pickup vehicle k at the cross-dock (k [ KS)
# DT_l :: Starting time of delivery vehicle l at the cross dock (l [ KD)
@variable(cd_reroute_mod, rt_k[k=K_S] >=0)
@variable(cd_reroute_mod, dt_k[k=K_D] >=0)


# ****************************************
# Constraints
#*****************************************

# --------------------[Pick-up Process]-----------------------------------------
# --------------------[new vehicles starting from depot]-----------------------------------------
# set P as the unserved nodes for pickup processing
P = P_un
# Each node is visited once by one vehicle
@constraint(cd_reroute_mod, p_one_veh[i= P],
                    sum(y_ik[i,k] for k=K_S) == 1)

# [new vehicle:: K_S_new]
@constraint(cd_reroute_mod, p_re_nd_entry[j= P, k= K_S_new],
                    sum(x_ijk[i,j,k] for i=vcat(cd_pick_start,P) if i != j)
                        == y_ik[j,k]
            )
# [new vehicle:: K_S_new]
@constraint(cd_reroute_mod, p_re_nd_exit[i= P, k=K_S_new],
                    sum(x_ijk[i,j,k] for j=vcat(cd_pick_end,P) if i != j)
                        == y_ik[i,k]
            )

# [new vehicle:: K_S_new]:: flow conservation
@constraint(cd_reroute_mod, p_re_mvmt[l= P, k=K_S_new],
                    sum(x_ijk[i,l,k] for i=vcat(cd_pick_start,P) if i != l)
                    - sum(x_ijk[l,j,k] for j=vcat(cd_pick_end,P) if l != j)
                    == 0
            )

# [new vehicle:: K_S_new]::
# if vehicle is used its route starts from CD
@constraint(cd_reroute_mod, p_re_newveh_start[i= cd_pick_start, k = K_S_new],
                    sum(x_ijk[i,j,k] for j=P) <= y_ik[i,k])

# [new vehicle:: K_S_new]::
# if vehicle starts from CD, vehicle routes ends at CD after routing
@constraint(cd_reroute_mod, p_re_newveh_end[j= cd_pick_end, k = K_S_new],
                sum(x_ijk[i,j,k] for i=P) <= y_ik[j,k])

# Capacity Constraints
@constraint(cd_reroute_mod, p_newveh_cap[k=K_S_new],
                    sum(p_i[i] * y_ik[i,k] for i=P) <= Q)




# -----------[pickup process:: time constraints]
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(cd_reroute_mod, p_veh_arr[i=vcat(cd_pick_start,P), j=vcat(P,cd_pick_end), k=K_S; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )

# Compute the departure time to zero if not visited
@constraint(cd_reroute_mod, p_init_dep_time[i=vcat(cd_pick_start, P, cd_pick_end), k=K_S],
                atp_ik[i,k] <= M*y_ik[i,k])


# --------------------[delivery Process]-----------------------------------------
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



# -----------[At Cross Dock:: time constraints]
# release time of veh k (in pickup process)
@constraint(cd_reroute_mod, release_time[i=P, j= cd_pick_end, k=K_S],
                rt_k[k] >= atp_ik[i,k] + t_ij[i,j] + u - M*(1-x_ijk[i,j,k]) )


# delivery vehicle (l) cannot start unless the pickup vehicle (k)
# is unloaded (greater than release time of pickup vehicle)
# delivery route (h) starts if pickup vehicle brings the product from (i)
@constraint(cd_reroute_mod, dept_time[k=K_S,l=K_D, i=P,h=D],
                dt_k[l]*CV[i,h] >= rt_k[k] - M*(2 - y_ik[i,k] - CV[i,h]*v_ik[h,l]))


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
)


print("------------------------starting optimization--------------------------")
optimize!(cd_reroute_mod)


if termination_status(cd_reroute_mod) == MOI.OPTIMAL
    optimal_solution = value.(x_ijk)
    optimal_objective = objective_value(cd_reroute_mod)
elseif termination_status(cd_reroute_mod) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x_ijk)
    suboptimal_objective = objective_value(cd_reroute_mod)
else
    error("The model was not solved correctly.")
end
