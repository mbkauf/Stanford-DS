#run the following only if this is the first time using this script/Julia
#this will be slow
using Pkg
Pkg.add("CSV")
Pkg.add("DifferentialEquations")
Pkg.add("DataFrames")
Pkg.add("StructArrays")
Pkg.add("Tables")
Pkg.add("Gadfly")

#load packages - this will be slow the 1st time you run it (and even a little slower than you might expect after that)
using CSV
using DifferentialEquations
using DataFrames
using StructArrays
using Tables
using Gadfly

params = DataFrame(beta=1.5, gamma=1, omega=0.5) #dataframes are mutable but all columns must be the same length
#if loading from file: params=CSV.file("params.csv") |> DataFrame
params_tuple=Tables.rowtable(params)[1] #named tuples are immutable but can have differing "column" lengths (good if some but not all parameters are age-stratified)

S_0 = 0.99
I_0 = 0.01
R_0 = 0.0
initial = vcat(S_0, I_0, R_0)

function disease_model!(du, u, p, t)
  S=u[1]
  I=u[2]
  R=u[3]
  du[:]=zeros(size(u))

  #1. new infections from S compartment to I
  du[1] = -p.beta*S*I + du[1]
  du[2] = p.beta*S*I + du[2]

  #2. recovery to R compartments from I
  du[3] = p.gamma*I + du[3]
  du[2] = -p.gamma*I + du[2]

  #3. waning of immunity from R to S
  du[3] = -p.omega*R + du[3]
  du[1] = p.omega*R + du[1]
end

#run ODE solver
t = 100 #100 cycles
ts = range(0, stop=t, step=1)
prob = ODEProblem(disease_model!, initial, (0.0, t), params_tuple)
sol = solve(prob, dense=false, save_timeseries=false, tstops=ts, saveat=ts)

#clean up output from ODE solver
out = DataFrame()
namelist = [:S, :I, :R]
for (i, name) in enumerate(namelist)
    out[!,name] = [sol.u[j][i] for j in 2:length(sol.u)]
end
out.time = Array(range(1, stop=t, step=1))
out_graph = stack(out, namelist, [:time], variable_name=:compartment)

#graph % of pop in each compartment over time (graphing is slow the 1st time you run it each time you laod the script)
plot(out_graph, x=:time, y=:value, color=:compartment, Geom.line)
