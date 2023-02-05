using Catalyst, DifferentialEquations, Plots

rn = @reaction_network begin
    (kOn,kOff), A + DNA <--> A_DNA
    k, A_DNA --> A_DNA + RNA
    kT, RNA --> RNA + GFP
    deg_R, RNA --> 0
    deg_G, GFP --> 0
end kOn kOff k kT deg_R deg_G

tspan = (0.0, 300.0)
u0 = [:A => 1, :DNA => 1, :A_DNA => 0, :RNA => 0, :GFP => 0]
p = [:kOn => 0.0042*60*60, :kOff => 0.0038*60*60, :k => 250, :kT => 10.0, :deg_R => 0.8, :deg_G => 0.3]

dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())
@time sol = solve(jprob, SSAStepper())

plot(sol)
plot(sol,idxs=1, xlimit=(0.0,1.0))
plot(sol,idxs=2, xlimit=(0.0,1.0))
plot(sol,idxs=3, xlimit=(0.0,1.0))
plot(sol,idxs=4)
plot(sol,idxs=5)

gfp = sol[5,:]
histogram(gfp[50000:end])