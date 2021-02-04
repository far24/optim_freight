using JuMP, GLPK
model = Model(with_optimizer(GLPK.Optimizer))
@variable(model, x[1:2, 1:2], Bin)
optimize!(model)
filter(k -> value(x[k]) > 0.5, keys(x))



Containers.DenseAxisArray:

using JuMP, GLPK
model = Model(with_optimizer(GLPK.Optimizer))
@variable(model, x[1:2, [:A, :B]], Bin)
optimize!(model)
# The following doesn't work for some reason
filter(k -> value(x[k]) > 0.5, keys(x))
# Here is a work-around
filter(k -> value(x[k]) > 0.5, collect(keys(x)))







Containers.SparseAxisArray:

using JuMP, GLPK
model = Model(with_optimizer(GLPK.Optimizer))
@variable(model, x[i=1:2, j=i:2], Bin)
optimize!(model)
# The following doesn't work: it appears to be a bug.
filter(k -> value(x[k]) > 0.5, keys(x))
# Here is a work-around
filter(k -> value(x[k]) > 0.5, collect(keys(x.data)))
