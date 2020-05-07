include("utils.jl")

#########
# Setup #
#########
datafile = "datadir/T_N_fresh_run.dat"

N_sample = 3
N_resamples = 700
#system = :pp
system = :torus

N_sweeps_tau_T = 20000  # Used for tau(T), does not need so much
N_sweeps_eq_tau_T = 1000

N_sweeps_tau_N = 500000  # Used for tau(N) at Tc, need alot
N_sweeps_eq_tau_N = 1000

############################
## Build pairs of T and N ##
############################
T_N_10 = [(T, 10) for T in (0.2:0.05:1.15) * Tc]
T_N_30 = [(T, 30) for T in (0.2:0.05:1.15) * Tc]
N_sweeps_tau_T = fill(N_sweeps_tau_T, length(T_N_10)+length(T_N_30))
N_sweeps_eq_tau_T = fill(N_sweeps_eq_tau_T, length(T_N_10)+length(T_N_30))

Tc_N = [(Tc, N) for N in [2, 4, 8, 10, 20, 30]]
N_sweeps_tau_N = fill(N_sweeps_tau_N, length(Tc_N))
N_sweeps_eq_tau_N = fill(N_sweeps_eq_tau_N, length(Tc_N))

T_N_pairs = vcat(T_N_10, T_N_30, Tc_N)
N_sweeps = vcat(N_sweeps_tau_T, N_sweeps_tau_N)
N_sweeps_eq = vcat(N_sweeps_eq_tau_T, N_sweeps_eq_tau_N)

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
