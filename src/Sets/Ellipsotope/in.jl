function ∈(x::AbstractVector, E::Ellipsotope)
    c, G, p, I, A, b = E.center, E.generators, E.p_norm, E.index_set, E.A, E.b
    n, m = size(G)

    # define new constraints and solve emptiness
    A_int = vcat(A, G)
    b_int = vcat(b, x - c)

    E_int = Ellipsotope(c, G, p, I, A_int, b_int)

    empty = is_empty(E_int)
    return !empty
end