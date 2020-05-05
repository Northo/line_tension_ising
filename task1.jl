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

Tc = 2.26919 # Tc = 2.269 from analytical
N_sweeps = 50000
N_sweep_eq = 20000

# ## Simulate over T ##
# Nx, Ny = 32, 32
# N = Nx*Ny
# H_pp = get_pp_hamiltonian(Nx, Ny)
# # Fix boundary conditions
# ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
# # Initial energies
# H_0_pp = calculate_energy(H_pp, ir, il, iu, id)
# #T = collect(0.2:0.2:1.3)
# T = collect(0.1:0.3:1.7)
# tau = simulate_over_T!(
#     T * Tc,
#     H_pp,
#     H_0_pp,
#     Nx,
#     Ny,
#     N_sweeps,
#     ir,
#     il,
#     iu,
#     id
# )

# println("Importing PyPlot...")
# using PyPlot
# println("PyPlot imported.")
# plt.plot(T, tau)
# plt.show()

## Simulate over N ##
Nx = [2, 4, 6, 8, 10, 20, 30, 40]
T = Tc
tau = simulate_over_N(Nx, N_sweeps, T)

println("Importing PyPlot...")
using PyPlot
println("PyPlot imported.")
plt.plot(Nx, tau)
plt.show()
