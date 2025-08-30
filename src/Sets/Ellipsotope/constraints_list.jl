"""
    constraints_list(E::Ellipsotope)

Return the list of constraints of `E`.
"""
function constraints_list(E::Ellipsotope)
    return (E.A, E.b)
end
