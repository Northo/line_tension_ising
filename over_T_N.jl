include("utils.jl")

#########
# Setup #
#########
datafile = "datadir/T_N_thursday_evening_fast.dat"

N_sample = 3
N_resamples = 700
#system = :pp
system = :torus

N_sweeps_tau_T = 100000  # Used for tau(T), does not need so much
N_sweeps_eq_tau_T = 1000

N_sweeps_tau_N = 500000  # Used for tau(N) at Tc, need alot
N_sweeps_eq_tau_N = 50000

############################
## Build pairs of T and N ##
############################
# T_N_10 = [(T, 10) for T in (0.2:0.05:1.15) * Tc]
# T_N_30 = [(T, 30) for T in (0.2:0.05:1.15) * Tc]
# N_sweeps_tau_T = fill(N_sweeps_tau_T, length(T_N_10)+length(T_N_30))
# N_sweeps_eq_tau_T = fill(N_sweeps_eq_tau_T, length(T_N_10)+length(T_N_30))

#T_N_pairs = [(Tc, N) for N in [2, 4, 6, 8, 10, 15, 20]]
#N_sweeps = Iterators.repeated(N_sweeps_tau_N)
#N_sweeps_eq = Iterators.repeated(N_sweeps_eq_tau_N)
T_N_pairs = [(Tc, 40), (Tc, 50)]
N_sweeps = Iterators.repeated(1000000)
N_sweeps_eq = Iterators.repeated(100000)

# T_N_pairs = vcat(T_N_10, T_N_30, Tc_N)
# N_sweeps = vcat(N_sweeps_tau_T, N_sweeps_tau_N)
# N_sweeps_eq = vcat(N_sweeps_eq_tau_T, N_sweeps_eq_tau_N)

# ### Mega ###
# T_N_pairs = [(Tc, 40), (Tc, 50), (Tc, 65)]
# N_sweeps = fill(1500000, 3)
# N_sweeps_eq = fill(20000, 3)

### Series of t ###
function T_from_t(t)
    return Tc*(1-t)
end

#T_N_pairs = [
    # (T_from_t(0.1), 20),
    # (T_from_t(0.1), 10),
    # (T_from_t(0.05), 30),
    # (T_from_t(0.05), 20),
    # (T_from_t(0.01), 20),
    # (T_from_t(0.01), 30),
    # (T_from_t(0.02), 20),
    # (T_from_t(0.02), 20),
    # (T_from_t(0.02), 30),
    #
    # (T_from_t(0.05), 30),
    # (T_from_t(0.1), 40),
    # (T_from_t(0.1), 60),
    # (T_from_t(0.05), 60),
    # (T_from_t(0.1), 100),
    # (T_from_t(0.05), 100),
    # (T_from_t(0.02), 120),
    # (T_from_t(0.02), 30),
    # (T_from_t(0.02), 20),
    # (T_from_t(0.02), 20),
    # (T_from_t(0.01), 18),
#]

#N_sweeps = Iterators.repeated(400000)
# N_sweeps = [
#     400000,
#     1000000,
#     400000,
#     800000,
# ]
# N_sweeps_eq = Iterators.repeated(1000)


##########
## Run! ##
##########
over_T_N(
    T_N_pairs,
    N_sweeps,
    N_sweeps_eq,
    system,
    N_resamples,
    datafile=datafile,
)
