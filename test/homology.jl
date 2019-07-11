@testset "Homology" begin
    cplx = SimplicialComplex(Simplex(1,2,3),
                        Simplex(2,4),
                        Simplex(3,4),
                        Simplex(5,4),
                        Simplex(6))

    h = homology(cplx, Int)
    @test eltype(h) == Tuple{Int64,Int64,Int64}
    @test grouptype(supertype(typeof(h))) == Nothing
    @test grouptype(typeof(h)) == Int
    @test length(h) == 3

    b, t, _ = group(h,0);
    @test b == 2
    @test t == 0

    itm, st = iterate(h)
    @test itm[1] == 0
    @test itm[2] == 2
    @test st[1] == 1

    itm, st = iterate(h, st)
    @test itm[1] == 1
    @test itm[2] == 1
    @test st[1] == 2

    itm, st = iterate(h, st)
    @test itm[1] == 2
    @test itm[2] == 0
    @test st[1] == 3

    @test iterate(h, st) === nothing

    @test ComputationalHomology.betti(h) == [2,1,0]
    @test ComputationalHomology.euler(h) == 1

    g = withgenerators(h)
    @test eltype(g) == Tuple{Int64,Int64,Int64,Dict{Chain,Int64}}
    @test length(g) == 3

    itm, st = iterate(g)
    @test itm[1] == 0
    @test itm[2] == 2
    @testset for (ch, d) in itm[4]
        @test d == 0
        @test ch[1][1] == 1
        @test ch[1][2] in [6,5]
    end

    itm, st = iterate(g, st)
    @test itm[1] == 1
    @test itm[2] == 1
    @test first(itm[4])[2] == 0
    @testset for (a,b) in zip( simplify(first(itm[4])[1]), Chain(1, Int) + (1, 4) + (1, 2) + (-1,3) + (-1,5) )
        @test a == b
    end

    itm, st = iterate(g, st)
    @test itm[1] == 2
    @test itm[2] == 0
    @test length(itm[4]) == 0

    @test iterate(g, st) === nothing

    @test length(generators(g)) == 3
end

@testset "Homology2" begin
    # data will contain all points {0,1}^d
    d = 3
    data = zeros(d,2^d)
    for i = 0:2^d-1
        s = bitstring(i)[end-d+1:end]
        s = split(s,"")
        s = parse.(Int,s)
        data[:,i+1] = s
    end

    cplx, w = vietorisrips(data, sqrt(3), maxoutdim=3)
    h = homology(cplx, Int)
    g = withgenerators(h)
    @test length(generators(g)[0]) == 1
    @test length(generators(g)[1]) == 0
    @test length(generators(g)[2]) == 0
    @test length(generators(g)[3]) == 0
end

@testset "Homology3" begin
    cplx = SimplicialComplex(Simplex(1,2), Simplex(2,3), Simplex(3,1))
    h = homology(cplx, Int)
    g = withgenerators(h)
    @test length(generators(g)[0]) == 1
    @test length(generators(g)[1]) == 1
end

@testset "Homology4" begin
    d = 2
    square = zeros(d,2^d)
    for i = 0:2^d-1
        s = bitstring(i)[end-d+1:end]
        s = split(s,"")
        s = parse.(Int,s)
        square[:,i+1] = s
    end
    square1 = [square square .+ [1;0] square .+ [2;0] square .+ [3;0]]#
    square2 = [square .+ [0;1] square .+ [3;1] square .+ [0;2] square .+ [3;2]]
    square3 = square1 .+ [0;3]
    data = [square1 square2 square3]
    data = unique(data, dims=2)

    cplx,w = vietorisrips(data, sqrt(2), maxoutdim=2)
    h = homology(cplx,Int)
    g = withgenerators(h)
    @test ComputationalHomology.betti(h) == [1;1;0]
end
