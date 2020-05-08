""" Metropolis Monte Carlo using the Mon-Jasnow algorithm
Written by: Thorvald M. Ballestad, May 2020.
This was written as a solution the exam in Computational Physics at NTNU.

The code is structured as follows:
Section one: helper functions used in the simulation itself
Section two: the heavy lifters of the simulation, ie. main part
Section three: helper functions such as writing to file, benchmark, analytical solutions, ...

Note that PyPlot is not included, because it takes too much time to import.
When using functions that require PyPlot, import it first with 'using PyPlot'.


Naming conventions:
 H - spin array
 ir, il, iu, id - neighbor index vectors, right, left, up, down respectively
 pp - original Mon-Jasnow algorithm, inspired by positive-positive
 torus - extended Mon-Jasnow algorithm

Central functions:
 simulate!
    Main function in the simulation.
    Takes some parameters that describes the system and returns the energy and energy difference as a function of time.
    Measures the difference in energy between the two systems (positive/negative or torus/klein) for every sweep
 sweep!
    Carries out one sweep, ie. NxN flips.
 step!
    Performs an actual flip attempt
 bootstrap!
    Calculates tau and uncertainty in tau given an energy difference vector
"""

using Statistics  # Mean and std

# Global constants
Tc = 2.26919 # from Onsager's solution

#####################################
## Section one -- Helper functions ##
#####################################
function get_random_hamiltonian(Nx::Int, Ny::Int)
    H = rand([1, -1], (Ny, Nx))
    return H
end


function get_zero_T_hamiltonian(Nx::Int, Ny::Int)
    H = ones(Int, (Ny, Nx))
    return H
end


function get_pp_hamiltonian(Nx, Ny; T=:inf)
    if T==:inf
        H_pp = get_random_hamiltonian(Nx+1, Ny)
    elseif T==:zero
        H_pp = get_zero_T_hamiltonian(Nx+1, Ny)
    else
        throw(ArgumentError("T must be :inf or :zero, got $(repr(T))"))
    end

    # Set fixed positive column
    H_pp[:, end] .= +1

    return H_pp
end


function get_torus_hamiltonian(Nx, Ny; T=:inf)
    if T==:inf
        H_torus = get_random_hamiltonian(Nx, Ny)
    elseif T==:zero
        H_torus = get_zero_T_hamiltonian(Nx, Ny)
    else
        throw(ArgumentError("T must be :inf or :zero, got $(repr(T))"))
    end

    return H_torus
end


function neighbor_interaction(H, i_y, i_x, ir, il, iu, id)
    """Given a spin array H, position indices (i_y, i_x), and index vectors, returns a vector of neighbor spins"""
    return [
        H[i_y, ir[i_x]],
        H[i_y, il[i_x]],
        H[iu[i_y], i_x],
        H[id[i_y], i_x],
    ]
end


function neighbor_interaction_sum(H, i_y, i_x, ir, il, iu, id)
    """Given a spin array H, position indices (i_y, i_x), and index vectors, returns a the sum of neighbor spins"""
    return H[i_y, ir[i_x]] +
        H[i_y, il[i_x]] +
        H[iu[i_y], i_x] +
        H[id[i_y], i_x]
end


function calculate_energy(H, ir, il, iu, id)
    """Finds total energy of spin array H"""
    total_energy = 0

    for i in CartesianIndices(H)
        s_i_y, s_i_x = Tuple(i)
        neighbors = neighbor_interaction_sum(H, s_i_y, s_i_x, ir, il, iu, id)
        total_energy -= 1/2 * flipsign(neighbors, H[i])
    end
    return total_energy
end


function pp_pn_difference(H, Ny, Nx)
    """Calculates the difference between a positive-positive and positive-negative system, ie. H_+- - H_++"""
    return 2*sum(H[:, Nx])
end


function torus_klein_difference(H, Ny, Nx)
    """Calculates the difference between a torus and klein system, ie. H_k - H_t"""
    # sum = 0
    # for y in 1:Ny
    #     sum -= H[y, 1]*(H[y, Nx] + H[Ny+1-y, Nx])
    # end
    return sum(
        H[:, 1] .*
        (H[:, Nx] .+ Iterators.reverse(H[:, Nx]))
    )
end


function get_index_vectors(Nx, Ny)
    """Create index lookup vectors.
    Note that the ends simple go out of bound,
    conditions such as periodic boundaries
    must be added"""

    ir = Vector{Integer}(undef, Nx)
    il = Vector{Integer}(undef, Nx)
    iu = Vector{Integer}(undef, Ny)
    id = Vector{Integer}(undef, Ny)

    for index in CartesianIndices((1:Ny, 1:Nx))
        y,x = Tuple(index)
        ir[x] = x+1
        il[x] = x-1
        iu[y] = y+1
        id[y] = y-1
    end

    return ir, il, iu, id
end


function get_pp_index_vectors(Nx, Ny)
    """Get index vector for the original Mon-Jasnow system of
    size Nx+2 x Ny.

    Note that the acutal spin matrix H that corresponds to such a system is only Nx+1 x Ny.
    Using the index vectors from this function, make the rightmost column appear as to be both right and left.
    """

    ir, il, iu, id = get_index_vectors(Nx+1, Ny)
    il[1] = Nx+1     # Positive column left
    ir[Nx] = Nx+1    # Positive column right
    ir[Nx+1] = Nx+1  # Links to itself
    iu[Ny] = 1       # Periodic in y
    id[1] = Ny

    return ir, il, iu, id
end


function get_torus_index_vectors(Nx, Ny)
    """Get index vectors for the torus system from the extended Mon-Jasnow"""

    ir, il, iu, id = get_index_vectors(Nx, Ny)
    # Connect the edges cyclic
    ir[Nx] = 1
    il[1] = Nx
    iu[Ny] = 1
    id[1] = Ny
    return ir, il, iu, id
end

function get_exponent_lookup(T)
    """Since e^(-beta Delta_H) can only have
    certain values, calculate the values first.

    Possible values for delta H are 4 and 8.
    We look up delta_H/4, so that is 1 and 2"""
    table = zeros(2)
    table[1] = exp(-4/T)
    table[2] = exp(-8/T)
    return table
end


######################################
## Section two -- The heavy lifters ##
######################################
function step!(
    H::AbstractArray,
    N::Int,
    T,
    exponent_lookup::AbstractVector,
    ir::AbstractVector,
    il::AbstractVector,
    iu::AbstractVector,
    id::AbstractVector,
)
    """Performs one step of the Minnapolis algorithm.
    Attempt to flip one spin
    """

    # This section is not as trivial as it seems
    # H may be bigger than Nx x Ny, to accomodate, for example,
    # fixed columns at the edges and such.
    # Julia has column oriented arrays, so
    # in linear coordinates, they go column first, then row.
    # Thus, so long as these "extra" columns, that make H bigger,
    # are only at the edge, Nx+1, Nx+2 ..., this will work.
    # Ny, must however be the correct size in H.
    s_i = rand(1:N)
    s_i_cartesian = CartesianIndices(H)[s_i]
    s_i_y, s_i_x = s_i_cartesian[1], s_i_cartesian[2]

    neighbor_sum = neighbor_interaction_sum(H, s_i_y, s_i_x, ir, il, iu, id)
    spin = H[s_i]

    # Delta_H can only take certain values
    # They are -8 -4 0 4 8.
    # -8, 4, and 0 give certain flip. 4 and 8 might give flip
    # To avoid calculating exp(-deltaH/T) each time,
    # we look up in exponent_lookup.
    # There, the index is the fourth of delta_H, so
    # delta_H = 8 corresponds to exponent_lookup[2]
    delta_H_fourth = flipsign(
        div(neighbor_sum, 2),  # div is integer division
        spin
    )
    accept = false
    if delta_H_fourth > 0
        r = rand()  # Random number r in [0,1)
        acceptance_criterion = exponent_lookup[delta_H_fourth]
        if acceptance_criterion > r
            accept = true
        end
    else
        accept = true
    end

    if accept == true
        H[s_i] = -spin
        return delta_H_fourth
    else
        return 0
    end
end


function sweep!(
    H,
    delta_H,
    N,
    T,
    exponent_lookup,
    ir, il, iu, id
)
    """Performs a sweep, more of a semantic function than anything else"""
    for i in 1:N
        delta_H += step!(H, N, T, exponent_lookup, ir, il, iu, id)
    end
    return delta_H*4  # Multiply by 4, because result from step is /4
end



function simulate!(
    H,
    Nx, Ny,
    T,
    N_sweeps,
    ir, il, iu, id;
    difference_function=pp_pn_difference
)
    """The acutal simulation"""
    N = Nx*Ny  # Used in random number generation inside step!
    delta_H = Vector{Int}(undef, N_sweeps)
    m = Vector{Int}(undef, N_sweeps)
    H_0 = calculate_energy(H, ir, il, iu, id)
    exponent_lookup = get_exponent_lookup(T)

    # Used per sweep, allocate now
    delta_H_sweep::Int = 0

    for i in 1:N_sweeps
        m[i] = difference_function(H, Ny, Nx)
        delta_H[i] = sweep!(
            H,
            delta_H_sweep,
            N,
            T,
            exponent_lookup,
            ir, il, iu, id,
        )
    end

    H_time = H_0 .+ cumsum(delta_H)
    return H_time, m
end


####################################
## Section three -- miscellaneous ##
####################################
function compare_initial_H(Nx, Ny, T, N_sweeps)
    """Compares a having an initial H with T=infty and T=zero,
    ie. having a random spin configuration
    and a completely ordered.
    Nice for investigating equilibration time"""
    ## Setup ##
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
    H_inf = get_pp_hamiltonian(Nx, Ny, T=:inf)
    H_zero = get_pp_hamiltonian(Nx, Ny, T=:zero)

    ## Simulate ##
    H_inf_time, m_inf = simulate!(H_inf, Nx, Ny, T, N_sweeps, ir, il, iu, id)
    H_zero_time, m_zero = simulate!(H_zero, Nx, Ny, T, N_sweeps, ir, il, iu, id)

    return H_inf_time, H_zero_time, m_inf, m_zero
end

function compare_and_plot_initial_H(Nx, Ny, T, N_sweeps)
    """Runs compare_initial_H and plots the results"""
    H_inf_time, H_zero_time, m_inf, m_zero = compare_initial_H(Nx, Ny, T, N_sweeps)
    plt.plot(H_inf_time, label="H++ inf")
    plt.plot(H_inf_time + m_inf, label="H+- inf")
    plt.plot(H_zero_time, label="H++ zero")
    plt.plot(H_zero_time + m_zero, label="H+- zero")
    plt.legend()
    plt.show()

    return H_inf_time, H_zero_time, m_inf, m_zero
end


function over_T_N(
    T_N_pairs,
    N_sweeps_array,
    N_sweeps_eq_array,
    systems,
    N_resamples;
    datafile,
    N_sample=3
)
    """Helper function for running simulate! for multiple values fo T and N
    Params:
     T_N_pairs - iterable with pairs of T and N
     N_sweeps_array - iterable with N_sweeps for each pair of T and N
     N_sweeps_eq_array - iterable with N_sweeps_eq for each pair of T and N
     systems - iterable with system type for each pair of T and N"""
    for (pair, N_sweeps, N_sweeps_eq, system) in zip(
        T_N_pairs,
        N_sweeps_array,
        N_sweeps_eq_array,
        systems,
    )
        T, N = pair

        if system == :pp
            H = get_pp_hamiltonian(N, N)
            ir, il, iu, id = get_pp_index_vectors(N, N)
            difference_function = pp_pn_difference
        elseif system == :torus
            H = get_random_hamiltonian(N, N)
            ir, il, iu, id = get_torus_index_vectors(N, N)
            difference_function = torus_klein_difference
        else
            throw(ArgumentError("Invalid system"))
        end

        H_time, m = simulate!(H, N, N, T, N_sweeps, ir, il, iu, id, difference_function=difference_function)

        print(".")  # Indicate simulate finished, only boot left
        N_tau, N_tau_std = bootstrap_tau(m[N_sweeps_eq:N_sample:end], T, N_resamples)
        print(".")
        tau = N_tau / N
        tau_std = N_tau_std / N

        write_T_N(
            T, N, tau, tau_std,
            N_sweeps=N_sweeps, N_sweeps_eq=N_sweeps_eq,
            t_sample=1,
            N_resamples=N_resamples,
            system=string(system),
            calculation_method="bootstrap",
            filename=datafile,
        )
        println("Wrote T=$T, N=$N")
    end
end


function benchmark(N_sweeps, T, Nx, Ny)
    """Used for benchmarking"""
    ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
    H = get_pp_hamiltonian(Nx, Ny)
    simulate!(H, Nx, Ny, T, N_sweeps, ir, il, iu, id)
end


function calculate_tau(diff, T, t_eq=1; t_sample=1)
    """Finds Tau*Ny"""
    # Ratio between partition functions
    party_ratio = exp.(-diff[t_eq:t_sample:end]/T)
    # Find expectation value
    party_ratio_mean = mean(party_ratio)
    tau = -T * log(party_ratio_mean)
    return tau
end


function bootstrap_tau(H_diff, T, N_resamples)
    """Applies the bootstrap on tau calculations
    H_diff is energy difference between two systems at equilibrium, with appropriate sample time.
    N_resamples is the number of resamples that are carried out in the bootstrap method.

    Returns:
     E<tau>, std(tau)
    """
    n = length(H_diff)
    party_ratio = exp.(-H_diff/T)
    taus = Vector(undef, N_resamples)
    tau_stds = Vector(undef, N_resamples)
    for i in 1:N_resamples
        party_ratio_sample = party_ratio[rand(1:n, n)]
        party_ratio_mean = mean(party_ratio_sample)
        taus[i] = party_ratio_mean
    end
    logs = log.(taus)
    tau = -T * mean(logs)
    tau_std = T * std(logs)

    return tau, tau_std
end


function write_T_N(
    T,
    N,
    tau,
    tau_std;
    N_sweeps,
    N_sweeps_eq,
    t_sample,
    N_resamples,  # If appliccable
    system,
    calculation_method,
    filename = "datadir/T_N_new.dat",
)
    # File format:
    # T N tau tau_std N_sweeps N_sweeps_eq t_sample system calculation_method N_resamples
    open(filename, "a") do file
        write(file, join([
            T, N, tau, tau_std,
            N_sweeps, N_sweeps_eq, t_sample,
            system,
            calculation_method, N_resamples,
        ], "\t"), "\n")
    end
end


function onsager(T)
    """Gives onsagers solution to tau(T/Tc)"""
    return 2 - Tc*T*log(coth(1/(T*Tc)))
end
