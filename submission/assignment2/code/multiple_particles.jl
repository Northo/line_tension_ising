include("utils.jl")
using Plots
using Random
using Statistics
using Printf

#Random.seed!(10)  # Used for debugging

D = 2E-4
T = 110

alpha = 0.2
t_end = 7000
dt = get_max_dt(D, alpha, 0.05)
t_space = 0:dt:t_end

N = 40

@printf("dt computed to %.2e\nNumber of steps: %i\n", dt, length(t_space))

x = Array{Float64, 2}(undef, N, length(t_space))
x[:, 1] .= 0

for i in eachindex(t_space[1:end-1])
    t = mod(t_space[i], T)
    if t < 3*T/4
        x[:, i+1] = x[:, i] + sqrt(2*D*dt)*randn(N)
    else
        x[:, i+1] = x[:, i] + sqrt(2*D*dt)*randn(N) - potential_force_time_independent.(x[:, i], alpha=alpha)*dt
    end
end
