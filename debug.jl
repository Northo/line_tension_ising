include("utils.jl")

#########
# Setup #
#########
N_sweeps = 100000
Nx, Ny = 20, 20
T = 1

#system = :pp
system = :torus


################
# Calculations #
################
if system == :pp
    system = get_pp_system(Nx, Ny)
elseif system == :torus
    system = get_torus_system(Nx, Ny)
else
    throw(ArgumentError("Invalid system"))
end

# Initial energies
@time simulate!(
    system,
    T,
    N_sweeps,
)
