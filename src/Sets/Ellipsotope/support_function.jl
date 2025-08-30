"""
    ρ(d::AbstractVector, E::Ellipsotope)

Return the support function of an ellipsotope.

### Input

- `d` -- direction
- `E` -- ellipsotope

### Output

The support function of the ellipsotope in the given direction. The p-norm
used is specified by the ellipsotope's `p_norm` field.
"""
function ρ(d::AbstractVector, E::Ellipsotope)
    p = E.p_norm
    q = if p == 1
        Inf
    elseif isinf(p)
        1.0
    else
        p / (p - 1)
    end
    return dot(d, E.center) + sum(norm(E.generators[:, I]' * d, q) for I in E.index_set)
end
