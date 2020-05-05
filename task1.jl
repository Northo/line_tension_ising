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
N_sweeps = 10000
N_sweep_eq = 1000

#system = :pp
system = :torus

## Simulate over T ##
Nx, Ny = 32, 32
N = Nx*Ny

if system == :pp
    H = get_pp_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
elseif system == :torus
    H = get_random_hamiltonian(Nx, Ny)
    ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
else
    throw(ArgumentError("Invalid system"))
end

# Initial energies
H_0 = calculate_energy(H, ir, il, iu, id)

#T = collect(0.2:0.2:1.3)
T = collect(0.2:0.2:1.2)
tau = simulate_over_T!(
    T * Tc,
    H,
    H_0,
    Nx,
    Ny,
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
plt.axhline(y=0, linestyle=":", color="gray", linewidth=0.3)
plt.axvline(x=1, linestyle=":", color="gray", linewidth=0.3)
plt.show()

# ## Simulate over N ##
# Nx = [2, 4, 6, 8, 10, 20, 30, 40]
# T = Tc
# tau = simulate_over_N(Nx, N_sweeps, T, T_hamil=:zero, t_sample=100)

# println("Importing PyPlot...")
# using PyPlot
# println("PyPlot imported.")
# plt.plot(Nx, tau)
# plt.show()
