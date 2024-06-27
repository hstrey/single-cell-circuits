using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions

@parameters kOn kOff k kT deg_R deg_G
@variables t
@species A(t) DNA(t) A_DNA(t) RNA(t) GFP(t)

rxs = [(@reaction kOn, A + DNA --> A_DNA),
       (@reaction kOff, A_DNA --> A + DNA),
       (@reaction k, A_DNA --> A_DNA + RNA),
       (@reaction kT, RNA --> RNA + GFP),
       (@reaction deg_R, RNA --> 0),
       (@reaction deg_G, GFP --> 0)]

discrete_events = [[t == 100.0] => [k ~ 500]]

tspan = (0.0, 300.0)
u0 = [A => 1, DNA => 1, A_DNA => 0, RNA => 0, GFP => 0]
p = [kOn => 0.0042*60*60, kOff => 0.00038*60*60, k => 250, kT => 10.0, deg_R => 0.8, deg_G => 0.3]

@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, RNA, GFP], [kOn, kOff, k, kT, deg_R, deg_G]; discrete_events)

dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())

sol = solve(jprob, SSAStepper(); tstops = [100.0])

plot(sol)