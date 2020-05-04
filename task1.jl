include("utils.jl")

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
Nx, Ny = 5,5
N = Nx*Ny
T = 1
H = get_random_Hamiltonian(Nx+2, Ny)
exponent_lookup = get_exponent_lookup(T)

# Set second last column to 1 and last column to -1
H[:, end-1] .= 1
H[:, end] .= -1

# Fix boundary conditions
ir, il, iu, id = get_index_vectors(Nx, Ny)
il[1] = Nx+1  # positive column
iu[Ny] = 1
id[1] = Ny

#######
# Run #
#######
for i in 1:100
    step!(H, N, T, exponent_lookup, ir, il, iu, id)
end
