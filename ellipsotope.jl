using ..LazySets.JuMP: Model, @variable, @objective, @constraint, optimize!,
                      set_silent, termination_status, objective_value, value,
                      SecondOrderCone, OPTIMAL, ALMOST_OPTIMAL, fix, MOI

                      
"""
    _plot_ellipsotope_raytracing_2d(E::Ellipsotope; N=100)

Generate vertices for a 2D plot of an Ellipsotope using the ray-tracing method.

### Input

- `E` -- The Ellipsotope to plot.
- `N` -- (optional, default: 100) The number of rays to trace.

### Output

A list of 2D vertices on the boundary of the ellipsotope.
"""
function _plot_ellipsotope_raytracing_2d(E::Ellipsotope; N=100)
    @assert dim(E) == 2 "This ray-tracing method is for 2D ellipsotopes."
    if isempty(E)
        return Vector{Vector{Float64}}()
    end

    c, G, p, I, A, b = E.center, E.generators, E.p_norm, E.index_set, E.A, E.b
    m = ngens(E)

    # `x_feasible` inside E
    model_feas = Model(SCS.Optimizer)
    set_silent(model_feas)
    @variable(model_feas, β_feas[1:m])
    if !isempty(A)
        @constraint(model_feas, A * β_feas .== b)
    end
    for J_k in I
        if isempty(J_k)
            continue
        end
        β_Jk = β_feas[J_k]
        if p == 2
            @constraint(model_feas, [1.0; β_Jk] in SecondOrderCone())
        elseif isinf(p)
            @constraint(model_feas, -1.0 .<= β_Jk .<= 1.0)
        else 
            len_Jk = length(J_k)
            @variable(model_feas, r[1:len_Jk] >= 0)
            @constraint(model_feas, sum(r) <= 1.0)
            for i in 1:len_Jk
                @constraint(model_feas, [r[i], 1.0, β_Jk[i]] in MOI.PowerCone(1/p))
            end
        end
    end
    optimize!(model_feas)

    if termination_status(model_feas) ∉ (OPTIMAL, ALMOST_OPTIMAL)
        @warn "Could not find a feasible starting point for plotting."
        return Vector{Vector{Float64}}()
    end

    β_start = value.(β_feas)
    x_feasible = c + G * β_start

    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, β[1:m])
    @variable(model, λ >= 0)
    @objective(model, Max, λ)
    if !isempty(A)
        @constraint(model, A * β .== b)
    end
    for J_k in I
        if isempty(J_k)
            continue
        end
        β_Jk = β[J_k]
        if p == 2
            @constraint(model, [1.0; β_Jk] in SecondOrderCone())
        elseif isinf(p)
            @constraint(model, -1.0 .<= β_Jk .<= 1.0)
        else # 1 < p < Inf
            len_Jk = length(J_k)
            @variable(model, r[1:len_Jk] >= 0)
            @constraint(model, sum(r) <= 1.0)
            for i in 1:len_Jk
                @constraint(model, [r[i], 1.0, β_Jk[i]] in MOI.PowerCone(1/p))
            end
        end
    end

    # ray constraint: c + Gβ = x_feasible + λg
    @variable(model, g[1:2])
    @constraint(model, c + G * β .== x_feasible + λ * g)

    # loop through directions, solve, collect vertices
    vertices = Vector{Vector{Float64}}()
    angles = range(0, 2π, length=N)

    for θ in angles
        g_val = [cos(θ), sin(θ)]
        fix.(g, g_val)
        optimize!(model)

        if termination_status(model) ∈ (OPTIMAL, ALMOST_OPTIMAL)
            λ_opt = objective_value(model)
            push!(vertices, x_feasible + λ_opt * g_val)
        else
            @warn "Optimizer failed for one of the directions."
        end
    end

    return vertices
end

"""
    plot_recipe(E::Ellipsotope, ε)

Plot an Ellipsotope.

### Input

- `E` -- The Ellipsotope to plot
- `ε` -- (optional) The precision, currently ignored for Ellipsotopes.

### Output

A tuple of (x,y) coordinates for plotting.
"""
function plot_recipe(E::Ellipsotope, ε)
    if dim(E) != 2
        throw(ArgumentError("Can only plot 2-dimensional Ellipsotopes, but got dimension $(dim(E))."))
    end
    
    vertices_list = _plot_ellipsotope_raytracing_2d(E)

    if isempty(vertices_list)
        return [], []
    end

    v_matrix = hcat(vertices_list...)
    x_coords = v_matrix[1, :]
    y_coords = v_matrix[2, :]

    push!(x_coords, x_coords[1])
    push!(y_coords, y_coords[1])
    
    return x_coords, y_coords
end
