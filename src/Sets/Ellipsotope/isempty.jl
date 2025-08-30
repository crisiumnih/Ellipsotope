function load_SCS_is_empty()
    return quote
        function isempty(E::Ellipsotope)
            c, G, p, I, A, b = E.center, E.generators, E.p_norm, E.index_set, E.A, E.b
            m = ngens(E)

            if isempty(A)
                return false
            end

            model = Model(SCS.Optimizer)
            set_silent(model)
            
            @variable(model, β[1:m])
            
            @objective(model, Min, sum((A * β - b).^2))
            
            for J_k in I
                if !isempty(J_k)
                    @constraint(model, [1.0; β[J_k]] in SecondOrderCone())
                end
            end
            
            optimize!(model)
            
            status = termination_status(model)
            if status == OPTIMAL || status == ALMOST_OPTIMAL
                min_sq_error = objective_value(model)
                is_empty_result = min_sq_error > 1e-8
                return is_empty_result
            else
                return false
            end
        end
    end
end
