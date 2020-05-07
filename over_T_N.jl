include("utils.jl")

#########
# Setup #
#########
N_sample = 3
N_resamples = 700
#system = :pp
system = :torus

T_N_pairs = [
    (0.5*Tc, 10),
    (0.9*Tc, 10),
]

N_sweeps = fill(1000, length(T_N_pairs))
N_sweeps_eq = fill(100, length(T_N_pairs))

over_T_N(
    T_N_pairs,
    N_sweeps,
    N_sweeps_eq,
    system,
    N_resamples,
)
