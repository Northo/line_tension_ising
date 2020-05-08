include("utils.jl")

#########
# Setup #
#########
datafile = "datadir/T_N_friday_morning_fig2_v2.dat"

N_sample = 1
N_resamples = 700

############################
## Build pairs of T and N ##
############################


# # Fig. 1
# T_range = [0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3] * Tc
# N_range = [5, 30, 80]
# T_N_pairs = Iterators.product(T_range, N_range)
# N_sweeps = Iterators.repeated(1000000)
# N_sweeps_eq = Iterators.repeated(100000)
# systems = Iterators.repeated(:pp)

# # Fig. 2
# T_range = [Tc]
# N_range = [2,4,6,8,10,15,20,25,30]
# T_N_pairs = Iterators.product(T_range, N_range)
# N_sweeps = Iterators.repeated(1000000)
# N_sweeps_eq = Iterators.repeated(100000)
# systems = Iterators.repeated(:pp)

# Fig. 3
### Series of t ###
function T_from_t(t)
    return Tc*(1-t)
end
# t_range = [0.02, 0.05, 0.10]
# one_over_Nt_range = [1, 1.5, 2, 2.5, 3]

# T_N_pairs = []
# for t in t_range, Nt_inv in one_over_Nt_range
#     push!(T_N_pairs, (
#         T_from_t(t),
#         round(Integer, 1/(t*Nt_inv))
#     ))
# end

# N_sweeps = Iterators.repeated(1000000)
# N_sweeps_eq = Iterators.repeated(100000)

# t=0.1, 1/Nt = 1, 1.5, 2
# t=0.05, 1/Nt = 1,2,2.5
# t=0.02, 1/Nt = 2, 2.5, 3, 4
T_N_pairs = [
    # (T_from_t(0.1), 5),
    # (T_from_t(0.1), 7),
    # (T_from_t(0.1), 10),
    # #
    # (T_from_t(0.05), 10),
    # (T_from_t(0.05), 16),
    # (T_from_t(0.05), 5),
    # #
    # (T_from_t(0.05), 20),
    # (T_from_t(0.05), 18),
    # (T_from_t(0.05), 7),
    # (T_from_t(0.02), 16),
    # (T_from_t(0.02), 13),
    # (T_from_t(0.02), 10),
    #
    # (T_from_t(0.05), 8),
    # (T_from_t(0.05), 10),
    # (T_from_t(0.05), 20),
    #
    # (T_from_t(0.02), 5),
    # (T_from_t(0.02), 12),
    # (T_from_t(0.02), 20),
    # (T_from_t(0.02), 25),
    #
    # Thread 1
    # (T_from_t(0.02), 7), (T_from_t(0.02), 8), (T_from_t(0.02), 9), (T_from_t(0.02), 10), (T_from_t(0.02), 12),
    #
    # Thread 2
    #
    # (T_from_t(0.1), 20),
    # Thread 3
    #
    #(T_from_t(0.05), 6), (T_from_t(0.05), 8), (T_from_t(0.05), 11),
    # Thread 1
    # (Tc, 2), (Tc, 4), (Tc, 6), (Tc, 8), (Tc, 10),
    # Thread 2
    (Tc, 10), (Tc, 15), (Tc, 20), (Tc, 30),
]

N_sweeps = Iterators.repeated(10000000)
N_sweeps_eq = Iterators.repeated(500000)
systems = Iterators.repeated(:torus)

##########
## Run! ##
##########
over_T_N(
    T_N_pairs,
    N_sweeps,
    N_sweeps_eq,
    systems,
    N_resamples,
    datafile=datafile,
    N_sample=N_sample,
)
