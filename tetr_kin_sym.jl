using Random
using DifferentialEquations
using ModelingToolkit
using Plots

# TetR binding kinetic constants
k_on = 6e4 # M^-1 s^-1
k_off =0.015 # s^-1

# activator binding kinetic constants
Ak_on = 0.0042*60*60 # h^-1
Ak_off = 0.0038*60*60 # s^-1

@parameters t k m kt g
D = Differential(t)
@variables rna(t) gfp(t)
eqs = [D(rna) ~ k - m*rna,
       D(gfp) ~ kt*rna - g*gfp]

delta_t = 0.001
p11 = exp(-Ak_on*delta_t)
p22 = exp(-Ak_off*delta_t)
state = [0]
time = [0.0]

function affect!(integrator,u,p,ctx)
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

    integrator.p[p.k] = state[end]*250.0
    push!(time, integrator.t)
end

@named burst = ODESystem(eqs,t,[rna,gfp],[k,m,kt,g];discrete_events=[0.001 => (affect!,[],[k],nothing)])

u0 = [rna => 0.0, gfp => 0.0]
p = [k => 250, m => 0.8, kt => 10.0, g => 0.3]
prob = ODEProblem(burst,u0,(0.0,300.0),p)

sol = solve(prob,Tsit5())

plot(sol,xlim=(0,200))
savefig("tetr_kin_sym.png")
plot(time, state,xlim=(0,5),label="state")
savefig("tetr_kin_sym_state.png")
gfp = sol[2,:]

histogram(gfp[50000:end])
