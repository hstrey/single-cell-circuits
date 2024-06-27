using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions

rn = @reaction_network begin
    (kOn,kOff), A + DNA <--> A_DNA
    k, A_DNA --> A_DNA + RNA
    kT, RNA --> RNA + GFP
    deg_R, RNA --> 0
    deg_G, GFP --> 0
end

tspan = (0.0, 300.0)
u0 = [:A => 1, :DNA => 1, :A_DNA => 0, :RNA => 0, :GFP => 0]
p = [:kOn => 0.0042*60*60, :kOff => 0.00038*60*60, :k => 250, :kT => 10.0, :deg_R => 0.8, :deg_G => 0.3]

dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())
@time sol = solve(jprob, SSAStepper())

plot(sol)
plot(sol,idxs=1, xlimit=(0.0,1.0))
plot(sol,idxs=2, xlimit=(0.0,1.0))
plot(sol,idxs=3, xlimit=(0.0,1.0))
plot(sol,idxs=4)
gfp_plot = plot(sol,idxs=5)
savefig(gfp_plot,"gfp_catalyst.png")

t_inter = 0:0.001:300
gfp_int = linear_interpolation(sol.t, sol[5,:])
gfp = gfp_int.(t_inter)
rna_int = linear_interpolation(sol.t, sol[4,:])
rna = rna_int.(t_inter)

plot(rna)

rna_mean = mean(rna[25000:end])
rna_std = std(rna[25000:end])
rna_norm = Normal(rna_mean,rna_std)

histogram(rna[25000:end], normed=true)
plot!(80:0.1:350,pdf.(rna_norm,80:0.1:350))

plot(t_inter,gfp,xlim=(100,110))

