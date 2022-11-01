using Turing, CSV, XLSX, Plots, DataFrames, Statistics

data = DataFrame(XLSX.readtable("Autocorrelation_2022_10_25.xlsx","constitutive_time_trace"))
data2 = DataFrame(XLSX.readtable("Autocorrelation_2022_10_25.xlsx","negative_feedback_time_trace"))
# Ornstein-Uhlenbeck process
@model ou(rn,T,delta_t) = begin
    ampl ~ Uniform(0.0,4000.0)
    tau ~ Uniform(0.0,200.0)
    
    b = exp(-delta_t/tau)
    
    rn[1] ~ Normal(0,sqrt(ampl))
    
    for i=2:T
        rn[i] ~ Normal(rn[i-1]*b,sqrt(ampl*(1-b^2)))
    end
end

const_names = names(data)
dt = data[:,:timestamp][2]-data[:,:timestamp][1]
ampl_list = []
tau_list = []
ampl_std_list = []
tau_std_list = []
for name in const_names[2:end]
    mydata = data[!,name] .- mean(data[!,name])
    oumodel = ou(mydata,201,dt)
    s = sample(oumodel, NUTS(0.65),2000)
    push!(ampl_list,mean(s[:ampl]))
    push!(tau_list,mean(s[:tau]))
    push!(ampl_std_list,std(s[:ampl]))
    push!(tau_std_list,std(s[:tau]))
end

histogram(tau_list,bins=0:20)
mean(tau_list)

neg_names = names(data2)
dt = data2[:,:timestamp][2]-data2[:,:timestamp][1]
ampl2_list = []
tau2_list = []
ampl2_std_list = []
tau2_std_list = []
for name in neg_names[2:end]
    mydata = data2[!,name] .- mean(data2[!,name])
    oumodel = ou(mydata,201,dt)
    s = sample(oumodel, NUTS(0.65),2000)
    push!(ampl2_list,mean(s[:ampl]))
    push!(tau2_list,mean(s[:tau]))
    push!(ampl2_std_list,std(s[:ampl]))
    push!(tau2_std_list,std(s[:tau]))
end

histogram(tau2_list,bins=30)
mean(tau2_list)

dfconst = DataFrame(tau=tau_list,tau_std=tau_std_list,ampl=ampl_list,ampl_std=ampl_std_list)
CSV.write("const.csv",dfconst)

dfneg = DataFrame(tau=tau2_list,tau_std=tau2_std_list,ampl=ampl2_list,ampl_std=ampl2_std_list)
CSV.write("neg.csv",dfneg)
