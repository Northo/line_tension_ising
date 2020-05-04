"""Simulation flow:
 . Given a nx by ny array of elements with -1 or +1.
 . Select a element at random
 . For that element, look up Delta H
 . Compare Delta H to acceptance function.
    - If accepted do change, if not reject
"""

function get_random_Hamiltonian(Nx::Int, Ny::Int)
    H = rand(Int, (Ny, Nx))
    return H
end


function step!(
    H::AbstractArray,
    N::Int,
    T,
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
    s_i = rand(1:N)  # Spin to flip
    s_i_cartesian = CartesianIndices(H)[s_i]
    s_i_x, s_i_y = Tuple(s_i_cartesian)
    neighbors = [H[idx, idy] for (idx, idy) in [
        (ir[s_i_x], s_i_y),  # Right
        (il[s_i_x], s_i_y),  # Left
        (s_i_x, iu[s_i_y]),  # Up
        (s_i_x, id[s_i_y]),  # Down
    ]]
    r = rand()  # Random number r in [0,1)

    spin = H[s_i]
    # TODO: Looup exponential
    delta_H = 2*spin*sum(neighbors)
    acceptance_criterion = exp(-delta_H/T)
    if acceptance_criterion > r
        H[s_i] = -spin
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
        x,y = Tuple(index)
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

Nx = 5
Ny = 5
N = Nx * Ny
T = 1

H = get_random_Hamiltonian(Nx, Ny)
ir, il, iu, id = get_torus_index_vectors(Nx, Ny)
for i in 1:100
    step!(H, N, T, ir, il, iu, id)
end
