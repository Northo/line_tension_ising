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
T = 0.9 * 2.269  # Tc = 2.269 from analytical
println("T: $T")
N_steps = 200000

H_pp = get_random_Hamiltonian(Nx+2, Ny)
H_pn = copy(H_pp)

# Set second last column to 1 and last column to -1
H_pp[:, end-1] .= +1
H_pp[:, end] .= +1

H_pn[:, end-1] .= +1
H_pn[:, end] .= -1

# Fix boundary conditions
ir, il, iu, id = get_index_vectors(Nx+2, Ny)
il[1] = Nx+1  # positive column
il[Nx+1] = Nx+1
il[Nx+2] = Nx
ir[Nx] = Nx+2  # positive/negative column
ir[Nx+2] = Nx+2
ir[Nx+1] = 1
iu[Ny] = 1
id[1] = Ny

# Initial energies
H_0_pp = calculate_energy(H_pp, ir, il, iu, id)
H_0_pn = calculate_energy(H_pn, ir, il, iu, id)

#######
# Run #
#######
delta_H_pp, delta_H_pn = simulate!(H_pp, N, T, N_steps, ir, il, iu, id)

H_pp_time = H_0_pp .+ cumsum(delta_H_pp)*4
H_pn_time = H_0_pn .+ cumsum(delta_H_pp + delta_H_pn)*4

N_tau = calculate_tau(H_pp_time, H_pn_time, T, 100000)


#tau = N_tau / Ny
#println(tau, " : ", N_tau)
