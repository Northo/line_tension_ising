include("utils.jl")
using Plots
using Printf
using Random
using Statistics

#Random.seed!(10)  # Used for debugging

D = 3.3E-4
T = 230

alpha = 0.2
dt = get_max_dt(D, alpha, 0.05)
t_end = 7000.0
t = 0.0:dt:t_end

@printf("dt computed to %.2e\nNumber of steps: %i", dt, length(t))

x = Vector{Float64}(undef, length(t))
x[1] = 0

for i in eachindex(t[1:end-1])
    x[i+1] = euler_step(
        x[i],
        t[i],
        dt,
        D,
        alpha=alpha,
        T=T,
    )
end

#pyplot()
p = plot(t, x, title=@sprintf("D = %.2f, T = %.1f", D, T))
hline!([0.2:1:5], color="gray")
vline!(1:T:t_end, color="gray")
