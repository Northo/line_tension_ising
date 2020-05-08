include("utils.jl")

#########
# Setup #
#########
datafile = "datadir/T_N_friday_morning_fig1.dat"

N_sample = 3
N_resamples = 700

############################
## Build pairs of T and N ##
############################


# Fig. 1
T_range = [0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3] * Tc
N_range = [5, 30, 80]
T_N_pairs = Iterators.product(T_range, N_range)
N_sweeps = Iterators.repeated(1000000)
N_sweeps_eq = Iterators.repeated(100000)
systems = Iterators.repeated(:pp)

# # Fig. 2
# T_range = [Tc]
# N_range = [2,4,6,8,10,15,20,25,30]
# T_N_pairs = Iterators.product(T_range, N_range)
# N_sweeps = Iterators.repeated(1000000)
# N_sweeps_eq = Iterators.repeated(100000)
# systems = Iterators.repeated(:pp)

# # Fig. 3
# ### Series of t ###
# function T_from_t(t)
#     return Tc*(1-t)
# end
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
# systems = Iterators.repeated(:torus)

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
)
