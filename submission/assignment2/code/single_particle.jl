include("utils.jl")
using Plots
using Random

# Random.seed!(10)  # Used for debugging

end_time = 1000
T = 80
#D = 2E-4
D = 3.25e-4
alpha = 0.2
dt = get_max_dt(D, alpha, 0.05)


t_range = 0:dt:end_time
x = Vector{Float64}(undef, length(t_range))

x[1] = 0
for (i, t) in enumerate(t_range[1:end-1])
    x[i+1] = euler_step(x[i], t, dt, D, T=T, alpha=0.2)
end

plot(t_range, x)
