using LazySets, Test
using LinearAlgebra

@testset "Ellipsotope Tests" begin
    # Type
    c = [1.0, 2.0]
    G = [1.0 0.0; 0.0 1.0]
    p = 2.0
    I = [[1], [2]]
    E = Ellipsotope(c, G, p, I)
    @test E isa Ellipsotope{Float64, Vector{Float64}, Matrix{Float64}, Float64, Vector{Vector{Int}}}

    # dim 
    @test dim(E) == 2

    # center 
    @test center(E) == c

    # p_norm 
    @test LazySets.norm(E) == p

    # genmat 
    @test genmat(E) == G

    # ngens 
    @test ngens(E) == 2

    # # index_set function
    # @test LazySets.   (E) == I

    # constraints 
    A = [1.0 1.0]
    b = [1.0]
    E_constrained = Ellipsotope(c, G, p, I, A, b)
    @test constraints_list(E_constrained) == (A, b)

    # isempty 
    @test !isempty(E)

    # in 
    x = [1.5, 2.5]
    @test x ∈ E
    x_out = [3.0, 3.0]
    @test x_out ∉ E
end
