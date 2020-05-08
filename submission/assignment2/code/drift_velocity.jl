include("utils.jl")
using PyPlot
using Random
using Statistics
using Printf

#Random.seed!(10)  # Used for debugging

#D = 2E-4
D = 3.25e-4
#Ts = [190, 200, 210]
Ts = 170:10:220

alpha = 0.2
t_end = 4000
num_iterations = 400  # Number of iterations for each period
dt = get_max_dt(D, alpha, 0.1)
t = 0:dt:t_end

setup_info = format_info(Dict([
    ("End time", t_end),
    ("steps each relization", length(t)),
    ("dt", dt),
    ("alpha", alpha),
    ("D", D),
    ("Ts", Ts),
    ("num iterations", num_iterations),
]))
print(setup_info)

x = Vector{Float64}(undef, length(t))
x[1] = 0

averages = Vector(undef, length(Ts))
stds = Vector(undef, length(Ts))

for (i_T, T) in enumerate(Ts)
    totals = Vector(undef, num_iterations)
    for iteration in 1:num_iterations
        @printf("\rT=%i [%s%s]", T, lpad('.', iteration, '.'), rpad(' ', 1+num_iterations-iteration)[2:end])
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
        totals[iteration] = x[end]
    end
    println()
    averages[i_T] = mean(totals)
    stds[i_T] = std(totals)
end

println(averages)
println(stds)
plt.errorbar(Ts, averages, yerr=stds, fmt="o")
plt.show()

# average_lengths = Vector(undef, length(Ts))
# for i in eachindex(Ts)
#     average_lengths[i] = mean(x[end, i:i+num_iterations])
#     println(average_lengths[i])
#     plt.scatter(average_lengths[i], 0, label=Ts[i])
# end

# plt.axvline(color="gray", linestyle="dotted")
# plt.legend()
# plt.plot()

