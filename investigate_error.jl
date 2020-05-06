include("utils.jl")
using Statistics  # Mean
#using PyPlot

#########
# Setup #
#########

#system = :pp
system = :torus

N_sweeps = 10000
N_sweep_eq = 5000
N_sample = 10  # Measure per N_sample sweep
N_resamples = 100  # Number of random samples to use in Bootstrap
Nx, Ny = 32, 32
T = 1 * Tc

if system == :pp
    H = get_pp_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
    difference_method = pp_pn_difference
elseif system == :torus
    H = get_random_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
    difference_method = torus_klein_difference
else
    throw(ArgumentError("Invalid system"))
end


# M is the difference between the two systems
_, m = simulate!(
    H,
    Nx,
    Ny,
    T,
    N_sweeps,
    ir, il, iu, id,
    difference_function=difference_method
)

eq_samples = m[N_sweep_eq:N_sample:end]

bootstrap_tau(eq_samples, T, N_resamples)
