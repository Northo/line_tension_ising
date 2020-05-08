include("utils.jl")
using PyPlot
using Random
using Statistics
using Printf

#Random.seed!(10)  # Used for debugging

D = 3.25e-4
T = 70

alpha = 0.2
t_end = 1000
num_particles = 5000
dt = get_max_dt(D, alpha, 0.1)
t_space = 0:dt:t_end

setup_info = format_info(Dict([
    ("End time", t_end),
    ("dt", dt),
    ("alpha", alpha),
    ("D", D),
]))
print(setup_info)

x = zeros(Float16, length(t_space), num_particles)
for (i,t) in enumerate(t_space[1:end-1])
    x[i+1, :] = euler_step.(x[i, :], t, dt, D, alpha=alpha, T=T)
end

sub_tspace = t_space[1:10:end]
sub_x = x[1:10:end, :]
plt.plot(sub_tspace, std(sub_x, dims=2))
plt.plot(sub_stapce, 0.035*sqrt.(sub_tspace))  # Magic number to scale
plt.show()
