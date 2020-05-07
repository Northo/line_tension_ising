include("utils.jl")

#########
# Setup #
#########
N_sweeps = 2000000
N_sweeps_eq = 4000
N_sample = 3
N_resamples = 700  # Bootstrap
T = Tc
Nx = [2, 4, 6, 8, 10, 20, 25]
#system = :pp
system = :torus

plot = true
write_to_file = true
filename = string(
    "N_",
    N_sweeps,
    "_eq", N_sweeps_eq,
    "_sys", system,
)
println(filename)

################
# Calculations #
################
tau, tau_std = simulate_over_N(
    Nx,
    N_sweeps,
    N_sweeps_eq,
    T,
    T_hamil=:zero,
    t_sample=N_sample,
    bootstrap=true,
    system=:torus,
    N_resamples=N_resamples,
)

if write_to_file
    write_tau_X(Nx, tau, filename, tau_std)
end

if plot
    println("Importing PyPlot...")
    using PyPlot
    println("PyPlot imported.")
    plot_tau_X(Nx, tau, tau_std)
end
