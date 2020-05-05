"""Simulation flow:
 . Given a nx by ny array of elements with -1 or +1.
 . Select a element at random
 . For that element, look up Delta H
 . Compare Delta H to acceptance function.
    - If accepted do change, if not reject
"""

function get_random_Hamiltonian(Nx::Int, Ny::Int)
    H = rand([1, -1], (Ny, Nx))
    return H
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
        if s_i_x == Nx
            return delta_H_fourth, -4*spin
        else
            return delta_H_fourth, 0
        end
    else
        return 0, 0
    end
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
