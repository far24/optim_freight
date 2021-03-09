# ****************************************
# Create a JuMP model
# ****************************************
del_mod = Model(GLPK.Optimizer)

# ****************************************
# Decision variables
#*****************************************

# BINARY
# ^^^^^^^^^
# see Nikolopoulou et al. (2016)
# z_ijk: 1 if veh k travels from node i to node j
@variable(del_mod, z_ijk[i= vcat(cd_del_start, D, cd_del_end) ,j= vcat(cd_del_start, D, cd_del_end),k=K_D], Bin)

# v_ik: 1 if node i is serviced by veh k
@variable(del_mod, v_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D], Bin)

# CONTINUOUS
# ^^^^^^^^^^^^
# atd_ik :: Time at which vehicle k leaves pickup node i (i [ S, k [ KS)
@variable(del_mod, atd_ik[i=vcat(cd_del_start, D, cd_del_end), k=K_D] >=0)
# l_ik = tardiness (tardiness/lateness = upper TW (b_i) - arrival time (s_ik)) at node i by vehicle k ==> this would be a non-linear constraint
# l_i = tardiness at node i
@variable(del_mod, l_ik[i=vcat(cd_del_start,D,cd_del_end), k = K_D] >= 0)


# ****************************************
# Constraints
#*****************************************

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

# -----------[delivery process:: time constraints]
# time a pickup vehicle k arrives at pickup node j after visiting pickup node i
@constraint(del_mod, d_veh_arr[i=vcat(cd_del_start,D), j=vcat(D,cd_del_end), k=K_D; i !=j],
                    atd_ik[j,k]
                    >= atd_ik[i,k] + sp_i + t_ij[i,j]  - M*(1-z_ijk[i,j,k])
            )

# Compute the departure time to zero if not visited
@constraint(del_mod, d_init_dep_time[i=D, k=K_D],
                atd_ik[i,k] <= M*v_ik[i,k])

b_i = [0,0,50,40,100,100,100,100]
# calculate the tardiness of vehicle arrival at node i
@constraint(del_mod, tard_del[j=vcat(cd_del_start, D, cd_del_end), k=K_D],
    l_ik[j,k]
    >= atd_ik[j,k] - b_i[j]
    - M*(1-v_ik[j,k])
    )



# ****************************************
# Objective
#*****************************************
@objective(del_mod, Min,
sum(t_ij[i,j] * z_ijk[i,j,k] for i=vcat(cd_del_start,D) for j=vcat(D,cd_del_end) for k=K_D)
   +sum(9999 * l_ik[i,k] for i=vcat(cd_del_start,D,cd_del_end) for k=K_D)
)


print("------------------------starting optimization--------------------------")
optimize!(del_mod)

if termination_status(del_mod) == MOI.OPTIMAL
    optimal_solution = value.(z_ijk)
    optimal_objective = objective_value(del_mod)
elseif termination_status(del_mod) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(z_ijk)
    suboptimal_objective = objective_value(del_mod)
else
    error("The model was not solved correctly.")
end



for k=K_D
    print("\nveh: ", k, "\t")
    for i=vcat(cd_del_start,D)
        for j=vcat(D,cd_del_end)
            if value.(z_ijk[i,j,k]) == 1
                print(i, "-->", j, "\t")
            end
        end
    end
end
