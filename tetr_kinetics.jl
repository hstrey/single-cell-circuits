using Random
using DifferentialEquations
using DiffEqCallbacks
using Plots
using BenchmarkTools

# TetR binding kinetic constants
k_on = 6e4 # M^-1 s^-1
k_off =0.015 # s^-1

# activator binding kinetic constants
Ak_on = 0.0042*60*60 # h^-1
Ak_off = 0.0038*60*60 # s^-1

km1 = 0.0
km2 = 250.0
kt = 10.0
m = 0.8
g = 0.3

function bursting!(du,u,p,t)
    k, m, kt, g = p # transcription, RNA decay, translation, protein decay
    rna, gfp = u
    du[1] = k - m*rna
    du[2] = kt*rna - g*gfp
end

global state = [0]
global time = [0.0]

delta_t = 0.001
p11 = exp(-Ak_on*delta_t)
p22 = exp(-Ak_off*delta_t)
state_time = 0.001:0.001:300

function affect!(integrator)
    if state[end]==0
        if rand()>p11
            push!(state,1)
        else
            push!(state,0)
        end
    else
        if rand()>p22
            push!(state,0)
        else
            push!(state,1)
        end
    end
    integrator.p[1] = state[end]*km2
    push!(time, integrator.t)
end

cb = PresetTimeCallback(state_time,affect!)

u0 = [0.0,0.0]
p = [km1,m,kt,g]
prob = ODEProblem(bursting!,u0,(0.0,300.0),p)

sol = solve(prob,Tsit5(),callback = cb)

plot(sol,xlim=(20,30))
plot!(time, state*4000)

gfp = sol[1,:]

histogram(gfp[50000:end])

