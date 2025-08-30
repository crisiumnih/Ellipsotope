struct Ellipsotope{N<:AbstractFloat,
    VN<:AbstractVector{N},
    MN<:AbstractMatrix{N},
    PN<:Real,
    IV<:AbstractVector{<:AbstractVector{Int}}} <: AbstractCentrallySymmetric{N}

    center::VN
    generators::MN
    p_norm::PN
    index_set::IV
    A::MN  # Constraint matrix A
    b::VN  # Constraint vector b

    function Ellipsotope(c::VN, G::MN, pn::PN, index_set::IV, A::MN, b::VN) where {
                        N<:AbstractFloat,
                        VN<:AbstractVector{N},
                        MN<:AbstractMatrix{N},
                        PN<:Real,
                        IV<:AbstractVector{<:AbstractVector{Int}}}
        n = size(G, 1)
        m = size(G, 2)

        @assert length(c) == n "Center vector dimension ($(length(c))) " *
                                "does not match generator matrix row dimension ($n)."
        @assert pn >= 1 "p must be >= 1, got p = $pn"
        return new{N, VN, MN, PN, IV}(c, G, pn, index_set, A, b)
    end
end

# creating an ellipsotope when you only have c, G, p, and I.
function Ellipsotope(c::AbstractVector{N}, G::AbstractMatrix{N}, p::Real, index_set::AbstractVector{<:AbstractVector{Int}}) where {N<:AbstractFloat}

    A = Matrix{N}(undef, 0, size(G, 2))
    b = Vector{N}(undef, 0)
    return Ellipsotope(c, G, p, index_set, A, b)

end

#  "Basic" Ellipsotope Constructor
function Ellipsotope(c::AbstractVector{N}, G::AbstractMatrix{N}, p::Real) where {N<:AbstractFloat}
    m = size(G, 2)
    default_index_set = m > 0 ? [collect(1:m)] : Vector{Vector{Int}}()
    return Ellipsotope(c, G, p, default_index_set)
end

"""
    Ellipsotope(Z::Zonotope)

Convert a Zonotope to an Ellipsotope.
"""
function Ellipsotope(Z::Zonotope)
    c = center(Z)
    G = genmat(Z)
    m = size(G, 2)
    
    # J = {{1}, {2}, ..., {m}}
    J_zonotope = [[i] for i in 1:m]
    p_default = 2.0 # p doesn't matter, can be anything >= 1
    
    return Ellipsotope(c, G, p_default, J_zonotope)
end

"""
    Ellipsotope(E::Ellipsoid)

Convert a LazySets.Ellipsoid to a generator-based Ellipsotope.
An ellipsoid is the affine map of a Euclidean ball (l_2 ball).
"""
function Ellipsotope(E::Ellipsoid)
    c = center(E)
    Q = E.shape_matrix
    G = cholesky(Q).L # Or any other matrix square root
    
    return Ellipsotope(c, G, 2.0) # p=2, and default index set {{1...m}}
end