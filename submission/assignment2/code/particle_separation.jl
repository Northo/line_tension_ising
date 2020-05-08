include("utils.jl")
using PyPlot
using Random
using Statistics
using Printf

#Random.seed!(10)  # Used for debugging

D = 3.25e-4

T = 70

alpha = 0.2
t_end = 10000
num_particles = 1000  # Number of iterations for each period
dt = get_max_dt(D, alpha, 0.1)
t_space = 0:dt:t_end

setup_info = format_info(Dict([
    ("End time", t_end),
    ("dt", dt),
    ("alpha", alpha),
    ("D", D),
]))
print(setup_info)

function run_sim(T, t_end)
    x = zeros(num_particles)
    for t in 0:dt:t_end
        x = euler_step.(x, t, dt, D, alpha=alpha, T=T)
    end
    return x
end

# run the simulation for two different radius, r and 3r, ie scale time

x_small = run_sim(T, t_end)
x_big = run_sim(3*T, t_end/3)

plt.hist(x_small, density=true, label="Small particle")
plt.hist(x_big, density=true, label="Big particle")
plt.xlabel("Total distance")
plt.legend()
plt.savefig("particle_separation.pdf")
