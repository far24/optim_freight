---
title: "Replicating VRPCD Paper"
author: "FA"
date: "11-23-2020"
---



# VRPCD
This scripts shows how to replicate a optimization paper.

````julia

using JuMP
using GLPK

model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
@variable(model, y >= 0)
@constraint(model, 6x + 8y >= 100)
@constraint(model, 7x + 12y >= 120)
@objective(model, Min, 12x + 20y)

optimize!(model)

@show value(x);
@show value(y);
@show objective_value(model);
````


````
value(x) = 14.999999999999993
value(y) = 1.2500000000000047
objective_value(model) = 205.0
````


