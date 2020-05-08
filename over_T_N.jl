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
T_range = [0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3] * Tc
N_range = [5, 30, 80]
T_N_pairs = Iterators.product(T_range, N_range)

N_sweeps = Iterators.repeated(1000000)
N_sweeps_eq = Iterators.repeated(100000)
systems = Iterators.repeated(:pp)

### Series of t ###
function T_from_t(t)
    return Tc*(1-t)
end

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
