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
Nx, Ny = 25, 25
N = Nx*Ny
T = 0.2

H_pp = get_random_Hamiltonian(Nx+2, Ny)
H_pn = get_random_Hamiltonian(Nx+2, Ny)

exponent_lookup = get_exponent_lookup(T)

# Set second last column to 1 and last column to -1
H_pp[:, end-1] .= +1
H_pp[:, end] .= +1

H_pn[:, end-1] .= +1
H_pn[:, end] .= -1

# Fix boundary conditions
ir, il, iu, id = get_index_vectors(Nx, Ny)
il[1] = Nx+1  # positive column
ir[Nx] = Nx+2  # positive/negative column
iu[Ny] = 1
id[1] = Ny

# Initial energies
H_0_pp = calculate_energy(H_pp, N, ir, il, iu, id)
H_0_pn = calculate_energy(H_pn, N, ir, il, iu, id)

#######
# Run #
#######
t_space = 1:200000
delta_H_pp = zero(t_space)
delta_H_pn = zero(t_space)

for i in t_space
    delta_H_pp[i] = step!(H_pp, N, T, exponent_lookup, ir, il, iu, id)
    delta_H_pn[i] = step!(H_pn, N, T, exponent_lookup, ir, il, iu, id)
end

H_pp = H_0_pp .+ (cumsum(delta_H_pp)*4)
H_pn = H_0_pn .+ (cumsum(delta_H_pn)*4)

party_ratio = exp.((H_pp - H_pn)/T)
party_ratio_exp_val = mean(party_ratio[150000:end])
tau = - T/Ny * log.(party_ratio_exp_val)
println(tau)
