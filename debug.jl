include("utils.jl")

#########
# Setup #
#########
N_sweeps = 100000
N_sweeps_eq = 1000
N_resamples = 200
Nx, Ny = 30, 30
T = 1

system = :pp
#system = :torus


################
# Calculations #
################
if system == :pp
    H = get_pp_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
    difference_function = pp_pn_difference
elseif system == :torus
    H = get_random_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
    difference_function = torus_klein_difference
else
    throw(ArgumentError("Invalid system"))
end

# Initial energies
@time H, m = simulate!(
    H,
    Nx, Ny,
    T,
    N_sweeps,
    ir, il, iu, id,
    difference_function=difference_function,
)

N_tau, N_tau_std = bootstrap_tau(m[N_sweeps_eq:end], T, N_resamples)

println(N_tau, " : ", N_tau_std)
