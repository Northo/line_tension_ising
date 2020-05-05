include("utils.jl")
using Statistics  # Mean

"""Original Mon Jasnow. ie.
Periodic in y, positive/negative boundaries in x

         /-----------.
  -------------      |
  +            +/-   |
  +            +/-   |
  +            +/-   |
Ny+            +/-   |
  +            +/-   |
  +            +/-   |
  -------------      |
       Nx   _________|
"""

#########
# Setup #
#########
Nx, Ny = 30, 30
N = Nx*Ny
Tc = 2.269  # Tc = 2.269 from analytical
N_sweeps = 10000
N_sweep_eq = 4000

H_pp = get_random_Hamiltonian(Nx+2, Ny)
H_pn = copy(H_pp)

# Set second last column to 1 and last column to -1
H_pp[:, end-1] .= +1
H_pp[:, end] .= +1

H_pn[:, end-1] .= +1
H_pn[:, end] .= -1

# Fix boundary conditions
ir, il, iu, id = get_pp_pn_index_vectors(Nx, Ny)

# Initial energies
H_0_pp = calculate_energy(H_pp, ir, il, iu, id)
H_0_pn = calculate_energy(H_pn, ir, il, iu, id)

#######
# Run #
#######
T = collect(0.2:0.2:1.5)
tau = simulate_over_T!(
    T * Tc,
    H_pp,
    H_0_pp,
    H_0_pn,
    N,
    N_sweeps,
    ir,
    il,
    iu,
    id
)

println("Importing PyPlot...")
using PyPlot
println("PyPlot imported.")
plt.plot(T, tau)
plt.show()

# T = 0.3 * Tc
# delta_H_pp, delta_H_pn = simulate!(H_pp, N, T, N_sweeps, ir, il, iu, id)
# H_pp_time = H_0_pp .+ cumsum(delta_H_pp)*4
# H_pn_time = H_0_pn .+ cumsum(delta_H_pp + delta_H_pn)*4
