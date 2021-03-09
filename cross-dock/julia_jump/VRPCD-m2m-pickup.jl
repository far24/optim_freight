K_S =[1,2]

# ****************************************
# Create a JuMP model
# ****************************************
pick_mod = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************

# BINARY
# ^^^^^^^^^
# see Nikolopoulou et al. (2016)
# x_ijk: 1 if veh k travels from node i to node j
@variable(pick_mod, x_ijk[i= vcat(cd_pick_start,P,cd_pick_end) ,j= vcat(cd_pick_start,P,cd_pick_end),k=K_S], Bin)

# y_ik: 1 if node i is serviced by veh k
@variable(pick_mod, y_ik[i= vcat(cd_pick_start,P,cd_pick_end), k=K_S], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# atp_ik :: Time at which vehicle k leaves pickup node i (i [ S, k [ KS)
# td_hl :: Time at which vehicle l leaves delivery node h (h [ D, l [ KD)
# RT_k :: Release time of pickup vehicle k at the cross-dock (k [ KS)
# DT_l :: Starting time of delivery vehicle l at the cross dock (l [ KD)
@variable(pick_mod, atp_ik[i=vcat(cd_pick_start,P,cd_pick_end), k=K_S] >=0)


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

# -----------[pickup process:: time constraints]
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(pick_mod, p_veh_arr[i=vcat(cd_pick_start,P), j=vcat(P,cd_pick_end), k=K_S; i !=j],
                    atp_ik[j,k]
                    >= atp_ik[i,k] + t_ij[i,j] - M*(1-x_ijk[i,j,k])
            )

# Compute the departure time to zero if not visited
@constraint(pick_mod, p_init_dep_time[i=vcat(cd_pick_start, P, cd_pick_end), k=K_S],
                atp_ik[i,k] <= M*y_ik[i,k])


# ****************************************
# Objective
#*****************************************
@objective(pick_mod, Min,
sum(t_ij[i,j] * x_ijk[i,j,k] for i=vcat(cd_pick_start,P) for j=vcat(P,cd_pick_end) for k=K_S)
#    +sum(t_ij[h,f] * z_hfl[h,f,l] for h=D for f=D for l=K_D)
)


print("------------------------starting optimization--------------------------")
optimize!(pick_mod)

if termination_status(pick_mod) == MOI.OPTIMAL
    optimal_solution = value.(x_ijk)
    optimal_objective = objective_value(pick_mod)
elseif termination_status(pick_mod) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x_ijk)
    suboptimal_objective = objective_value(pick_mod)
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
