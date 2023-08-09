using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions
using ModelingToolkit, Random
using Distributions: Geometric

@parameters k b deg_R deg_G kT
@variables t
@species RNA(t) GFP(t)
@register_symbolic Geometric(b)
m = rand(Geometric(1/b)) + 1

burstyrxs = [Reaction(k, nothing, [RNA], nothing, [m]),
			 Reaction(deg_R, [RNA], nothing),
             Reaction(deg_G, [GFP], nothing),
             Reaction(kT, [RNA],[RNA,GFP])]

@named burstygene = ReactionSystem(burstyrxs, t)

tspan = (0.0, 300.0)
u0 = [:RNA => 0, :GFP => 0]
p = [:b => 0.0042*60*60, :k => 1/(1/(0.0042*60*60)+1/(0.0038*60*60)), :kT => 10.0, :deg_R => 0.8, :deg_G => 0.3]

dprob = DiscreteProblem(burstygene, u0, tspan, p)
jprob = JumpProblem(burstygene, dprob, Direct())
@time sol = solve(jprob, SSAStepper())

plot(sol)