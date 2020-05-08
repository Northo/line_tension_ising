"""Simulation flow:
 . Given a nx by ny array of elements with -1 or +1.
 . Select a element at random
 . For that element, look up Delta H
 . Compare Delta H to acceptance function.
    - If accepted do change, if not reject
"""

using Statistics  # Mean
using DelimitedFiles

# Global constants
Tc = 2.26919 # Tc = 2.269 from analytical


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
        H_pp = get_random_hamiltonian(Nx+2, Ny)
    elseif T==:zero
        H_pp = get_zero_T_hamiltonian(Nx+2, Ny)
    else
        throw(ArgumentError("T must be :inf or :zero, got $(repr(T))"))
    end

    # Set fixed positive/negative columns
    H_pp[:, [end-1, end]] .= +1

    return H_pp
end


function get_torus_hamiltonian(Nx, Ny; T=:inf)
    if T==:inf
        H_pp = get_random_hamiltonian(Nx, Ny)
    elseif T==:zero
        H_pp = get_zero_T_hamiltonian(Nx, Ny)
    else
        throw(ArgumentError("T must be :inf or :zero, got $(repr(T))"))
    end

    return H_pp
end


function neighbor_interaction(H, i_y, i_x, ir, il, iu, id)
    return [
        H[i_y, ir[i_x]],
        H[i_y, il[i_x]],
        H[iu[i_y], i_x],
        H[id[i_y], i_x],
    ]
end


function neighbor_interaction_sum(H, i_y, i_x, ir, il, iu, id)
    return H[i_y, ir[i_x]] +
        H[i_y, il[i_x]] +
        H[iu[i_y], i_x] +
        H[id[i_y], i_x]
end


function calculate_energy(H, ir, il, iu, id)
    total_energy = 0

    for i in CartesianIndices(H)
        s_i_y, s_i_x = Tuple(i)
        neighbors = neighbor_interaction(H, s_i_y, s_i_x, ir, il, iu, id)
        total_energy -= 1/2 * flipsign(sum(neighbors), H[i])
    end
    return total_energy
end


function pp_pn_difference(H, Ny, Nx)
    """Calculates H_pn - H_pp"""
    return 2*sum(H[:, Nx])
end


function torus_klein_difference(H, Ny, Nx)
    """Calculates H_k - H_t"""
    # sum = 0
    # for y in 1:Ny
    #     sum -= H[y, 1]*(H[y, Nx] + H[Ny+1-y, Nx])
    # end
    return sum(
        H[:, 1] .*
        (H[:, Nx] .+ Iterators.reverse(H[:, Nx]))
    )
end


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

    ir, il, iu, id are lookups for index right, left, up, down

    TODO:
     . Consider if one should change order of operations,
    to reduce operations. Ie. calculate acceptance_criterion
    before r, then throw away at once if DH < 0.
    """

    # This section is not as trivial as it seems
    # H may be bigger than Nx x Ny, to accomodate, for example.
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
    for i in 1:N
        delta_H += step!(H, N, T, exponent_lookup, ir, il, iu, id)
    end
    return delta_H*4  # Multiply by 4, because result from step is /4
end



function simulate!(H, Nx, Ny, T, N_sweeps, ir, il, iu, id; difference_function=pp_pn_difference)
    N = Nx*Ny
    delta_H = zeros(N_sweeps)
    H_0 = calculate_energy(H, ir, il, iu, id)
    m = zeros(N_sweeps)
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


function compare_initial_H(Nx, Ny, T, N_sweeps)
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
    H_inf_time, H_zero_time, m_inf, m_zero = compare_initial_H(Nx, Ny, T, N_sweeps)
    plt.plot(H_inf_time, label="H++ inf")
    plt.plot(H_inf_time + m_inf, label="H+- inf")
    plt.plot(H_zero_time, label="H++ zero")
    plt.plot(H_zero_time + m_zero, label="H+- zero")
    plt.legend()
    plt.show()

    return H_inf_time, H_zero_time, m_inf, m_zero
end


function simulate_over_T!(
    T_range,
    H,
    Nx,
    Ny,
    N_sweeps,
    N_sweep_eq,
    ir, il, iu, id;
    N_resamples=100,
    t_sample=1,  # Step for selecting measurements
    difference_function=pp_pn_difference,
    bootstrap=false,
)
    tau_list = zero(T_range)
    if bootstrap
        tau_std_list = zero(T_range)
    end

    for (i, T) in enumerate(T_range)
        H_time, m = simulate!(H, Nx, Ny, T, N_sweeps, ir, il, iu, id, difference_function=difference_function)

        if bootstrap
            N_tau, N_tau_std = bootstrap_tau(m[N_sweep_eq:t_sample:end], T, N_resamples)
            tau_std = N_tau_std/Ny
            tau_std_list[i] = tau_std
        else
            N_tau = calculate_tau(m, T, N_sweep_eq, t_sample=t_sample)
        end
        tau = N_tau / Ny
        tau_list[i] = tau
        println("T: $T")
        if bootstrap
            println(" .tau: $tau, Ntau: $N_tau, std: $tau_std")
        else
            println(" .tau: $tau, Ntau: $N_tau")
        end
    end

    if bootstrap
        return tau_list, tau_std_list
    else
        return tau_list
    end
end


function simulate_over_N(
    Nx_range,
    N_sweeps,
    N_sweeps_eq,
    T;
    T_hamil=:inf,
    t_sample=1,
    system=:pp,
    bootstrap=false,
    N_resamples=100,  # For bootstrap
)
    N_tau_list = zeros(length(Nx_range))
    if bootstrap
        N_tau_std_list = zeros(length(Nx_range))
    end

    for (i, Nx) in enumerate(Nx_range)
        Ny = Nx
        N = Nx*Ny

        if system == :pp
            H = get_pp_hamiltonian(Nx, Ny)
            ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
            difference_function = pp_pn_difference
        elseif system == :torus
            H = get_random_hamiltonian(Nx, Ny)
            ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
            difference_function = torus_klein_difference
        else
            throw(ArgumentError("Invalid system"))
        end

        H_time, m = simulate!(H, Nx, Ny, T, N_sweeps, ir, il, iu, id, difference_function=difference_function)
        # Give feedback that simulation has finished and tau calculation begins
        print(".")
        if bootstrap
            N_tau, N_tau_std = bootstrap_tau(m[N_sweeps_eq:t_sample:end], T, N_resamples)
            N_tau_std_list[i] = N_tau_std
        else
            N_tau = calculate_tau(m, T, N_sweeps_eq, t_sample=t_sample)
        end

        tau = N_tau / Ny
        N_tau_list[i] = N_tau
        println("Nx: $Nx")
        println(" .tau: $tau, Ntau: $N_tau")
    end

    if bootstrap
        return N_tau_list, N_tau_std_list
    else
        return N_tau_list
    end
end


function over_T_N(
    T_N_pairs,
    N_sweeps_array,
    N_sweeps_eq_array,
    systems,
    N_resamples;
    datafile
)
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
        N_tau, N_tau_std = bootstrap_tau(m[N_sweeps_eq:end], T, N_resamples)
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


function get_index_vectors(Nx, Ny)
    """Create index lookup tables.
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
    """Get index vector for the ++/+- system of
    size Nx+2 x Ny.
    """

    ir, il, iu, id = get_index_vectors(Nx+2, Ny)
    il[1] = Nx+1     # Positive column
    il[Nx+1] = Nx+1  # Leftmost column links to itself
    il[Nx+2] = Nx    # Rightmost column
    ir[Nx] = Nx+2    # Positive/negative column
    ir[Nx+2] = Nx+2  # Rightmost column links to itself
    ir[Nx+1] = 1     # Leftmost column
    iu[Ny] = 1
    id[1] = Ny

    return ir, il, iu, id
end


function get_torus_index_vectors(Nx, Ny)
    ir, il, iu, id = get_index_vectors(Nx, Ny)
    # Connect the edges
    ir[Nx] = 1
    il[1] = Nx
    iu[Ny] = 1
    id[1] = Ny
    return ir, il, iu, id
end


function get_klein_index_vectors(Nx, Ny)
    ir, il, iu, id = get_index_vectors(Nx, Ny)
    # Connect the Mobius band
    for y in 1:Ny
        ir[y, Nx] = (Ny+1-y, 1)
        il[y, 1] = (Ny+1-y, Nx)
    end
    # Connect the periodic BC in y
    iu[Ny, :] = [(1, x) for x in 1:Nx]
    id[1, :] = [(Ny, x) for x in 1:Nx]
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


function bootstrap_tau(H_diff, T, N_resamples)
    """Applies the bootstrap on tau calculations
    H_diff is energy difference between two systems at equilibrium, with appropriate sample time.
    N_resamples is the number of resamples that are carried out in the bootstrap method.

    Returns:
     E<tau>, std(tau)
    """
    n = length(H_diff)
    party_ratio = exp.(-H_diff/T)
    party_ratio_samples = party_ratio[rand(1:n, (N_resamples, n))]
    party_ratio_mean = mean(party_ratio_samples, dims=2)
    log_of_mean = log.(party_ratio_mean)

    tau = -T * mean(log_of_mean)
    tau_std = T * std(log_of_mean)

    return tau, tau_std
end


function plot_tau_X(X, tau, tau_std=undef)
    if tau_std!=undef
        plt.errorbar(X, tau, yerr=tau_std)
    else
        plt.plot(X, tau)
    end
    plt.axhline(y=0, linestyle=":", color="gray", linewidth=0.3)
    plt.axvline(x=1, linestyle=":", color="gray", linewidth=0.3)
    plt.show()
end


function write_tau_X(X, tau, filename, tau_std=undef)
    if tau_std==undef
        writedlm(string("datadir/", filename, ".dat"), [X, tau])
    else
        writedlm(string("datadir/", filename, "_std.dat"), [X, tau, tau_std])


    end
end


function read_tau_X(filename, tau_std=false)
    if tau_std
        data = readdlm(string("datadir/", filename, "_std.dat"))
        X = data[1, :]
        tau = data[2, :]
        tau_std = data[3, :]
        return X, tau, tau_std
    else
        data = readdlm(string("datadir/", filename, ".dat"))
        X = data[1, :]
        tau = data[2, :]
        return X, tau
    end
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

# Nx = 5
# Ny = 5
# N = Nx * Ny
# T = 1

# H = get_random_Hamiltonian(Nx, Ny)
# exponent_lookup = get_exponent_lookup(T)
# ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
# for i in 1:100
#     step!(H, N, T, exponent_lookup, ir, il, iu, id)
# end
