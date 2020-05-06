include("utils.jl")

#########
# Setup #
#########
N_sweeps = 100000
N_sweeps_eq = 1000
N_sample = 3
N_resamples = 700 
Nx, Ny = 20, 20
T = collect(0.85:0.02:1.05) 
#T = vcat(0.2:0.2:0.8, 0.9:0.05:1.1, 1.2:0.1:1.4)  # Fraction of Tc

#system = :pp
system = :torus

plot = true
write_to_file = true
filename = string(
    "T_",
    N_sweeps,
    "_eq", N_sweeps_eq,
    "N", Nx,
    "_sys", system,
)
println(filename)

################
# Calculations #
################
if system == :pp
    H = get_pp_hamiltonian(Nx, Ny, T=:zero)
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
    difference_function = pp_pn_difference
elseif system == :torus
    H = get_random_hamiltonian(Nx, Ny, T=:zero)
    ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
    difference_function = torus_klein_difference
else
    throw(ArgumentError("Invalid system"))
end

# Initial energies
H_0 = calculate_energy(H, ir, il, iu, id)

tau, tau_std = simulate_over_T!(
    T * Tc,
    H,
    H_0,
    Nx,
    Ny,
    N_sweeps,
    N_sweeps_eq,
    ir,
    il,
    iu,
    id,
    difference_function=difference_function,
    bootstrap=true,
    t_sample=N_sample,
    N_resamples=N_resamples,
)

if write_to_file
    write_tau_X(T, tau, filename, tau_std)
end

if plot
    println("Importing PyPlot...")
    using PyPlot
    println("PyPlot imported.")
    plot_tau_X(T, tau, tau_std)
end
