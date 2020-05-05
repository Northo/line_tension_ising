"""Simulation flow:
 . Given a nx by ny array of elements with -1 or +1.
 . Select a element at random
 . For that element, look up Delta H
 . Compare Delta H to acceptance function.
    - If accepted do change, if not reject
"""

function get_random_hamiltonian(Nx::Int, Ny::Int)
    H = rand([1, -1], (Ny, Nx))
    return H
end


function get_pp_hamiltonian(Nx, Ny)
    H_pp = get_random_hamiltonian(Nx+2, Ny)
    # Set fixed positive/negative columns
    H_pp[:, [end-1, end]] .= +1

    return H_pp
end


function calculate_energy(H, ir, il, iu, id)
    total_energy = 0

    for i in CartesianIndices(H)
        s_i_y, s_i_x = Tuple(i)
        neighbors = [H[idy, idx] for (idy, idx) in [
            (s_i_y, ir[s_i_x]),  # Right
            (s_i_y, il[s_i_x]),  # Left
            (iu[s_i_y], s_i_x),  # Up
            (id[s_i_y], s_i_x),  # Down
        ]]
        total_energy -= 1/2 * flipsign(sum(neighbors), H[i])
    end
    return total_energy
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
     . Instead of evaluating exp(-DH/T), lookup preevaluated
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
    s_i = rand(1:N)  # Spin to flip
    s_i_cartesian = CartesianIndices(H)[s_i]
    s_i_y, s_i_x = Tuple(s_i_cartesian)
    neighbors = [H[idy, idx] for (idy, idx) in [
        (s_i_y, ir[s_i_x]),  # Right
        (s_i_y, il[s_i_x]),  # Left
        (iu[s_i_y], s_i_x),  # Up
        (id[s_i_y], s_i_x),  # Down
    ]]

    spin = H[s_i]
    neighbor_sum = sum(neighbors)

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
    delta_H_pp,
    N,
    T,
    exponent_lookup,
    ir, il, iu, id
)
    for i in 1:N
        delta_H_pp += step!(H, N, T, exponent_lookup, ir, il, iu, id)
    end
    return delta_H_pp*4  # Multiply by 4, because result from step is /4
end



function simulate!(H, Nx, Ny, T, N_sweeps, ir, il, iu, id)
    N = Nx*Ny
    delta_H_pp = zeros(N_sweeps)
    m = zeros(N_sweeps)
    exponent_lookup = get_exponent_lookup(T)

    # Used per sweep, allocate now
    delta_H_pp_sweep::Int = 0

    for i in 1:N_sweeps
        m[i] = 2*sum(H[:, Nx])
        delta_H_pp[i] = sweep!(
            H,
            delta_H_pp_sweep,
            N,
            T,
            exponent_lookup,
            ir, il, iu, id,
        )
    end

    return delta_H_pp, m
end


function simulate_over_T!(T_range, H, H_0_pp, Nx, Ny, N_sweeps, ir, il, iu, id)
    tau_list = zero(T_range)

    for (i, T) in enumerate(T_range)
        delta_H_pp, m = simulate!(H_pp, Nx, Ny, T, N_sweeps, ir, il, iu, id)
        H_pp_time = H_0_pp .+ cumsum(delta_H_pp)

        N_tau = calculate_tau(m, T, N_sweep_eq)
        tau = N_tau / Ny
        tau_list[i] = tau
        println("T: $T")
        println(" .tau: $tau, Ntau: $N_tau")

        H_0_pp = H_pp_time[end]
    end

    return tau_list
end


function simulate_over_N(Nx_range, N_sweeps, T)
    N_tau_list = zeros(length(Nx_range))
    for (i, Nx) in enumerate(Nx_range)
        Ny = Nx
        N = Nx*Ny
        ir, il, iu, id = get_pp_index_vectors(Nx, Ny)
        H_pp = get_pp_hamiltonian(Nx, Ny)

        H_0_pp = calculate_energy(H_pp, ir, il, iu, id)

        delta_H_pp, m = simulate!(H_pp, Nx, Ny, T, N_sweeps, ir, il, iu, id)
        H_pp_time = H_0_pp .+ cumsum(delta_H_pp)*4

        N_tau = calculate_tau(m, T, N_sweep_eq)
        tau = N_tau / Ny
        N_tau_list[i] = N_tau
        println("Nx: $Nx")
        println(" .tau: $tau, Ntau: $N_tau")
    end

    return N_tau_list
end


function calculate_tau(diff, T, t_eq)
    """Finds Tau*Ny"""
    # Ratio between partition functions
    party_ratio = exp.(-diff[t_eq:end]/T)
    # Find expectation value
    party_ratio_mean = mean(party_ratio)
    tau = -T * log.(party_ratio_mean)
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
    il[1] = Nx+1  # Positive column
    il[Nx+1] = Nx+1  # Leftmost column links to itself
    il[Nx+2] = Nx  # Rightmost column
    ir[Nx] = Nx+2  # Positive/negative column
    ir[Nx+2] = Nx+2  # Rightmost column links to itself
    ir[Nx+1] = 1  # Leftmost column
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
