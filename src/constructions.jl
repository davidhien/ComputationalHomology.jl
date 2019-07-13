"""Find all neighbors of vertex `u`` within `G` that precede it in the given ordering"""
function lowernbrs(cplx, u, E)
    # get all edges of u
    uval = first(u[:values])
    uidx = u[:index]
    # select {v} s.t. u > v for all edges {u,v}
    !any(E[:,uval]) && return celltype(cplx)[]
    return filter(v->v[:index] < uidx && E[first(v[:values]),uval], cells(cplx,0))
end

"""Incremental construction of VR complex from neighborhood graph"""
function inductive!(cplx, k, E)
    for i in 1:k
        cls = cells(cplx, i)
        cls === nothing && continue
        for τ in cells(cplx, i)
            N = celltype(cplx)[]
            τvals = τ[:values]
            for (j, uval) in enumerate(τvals)
                uidx = cplx[Simplex(uval), 0]
                u = cplx[uidx, 0]
                N = j == 1 ? lowernbrs(cplx, u, E) : intersect(N, lowernbrs(cplx, u, E))
            end
            for v in N
                first(v[:values]) in τvals && continue
                σ = Simplex(τvals..., v[:values]...)
                cplx[σ] > size(cplx, i+1) && push!(cplx, σ)
            end
        end
    end
end

"""Add all simplices of the `k`-skeleton whose maximal vertex is `τ`"""
function addcofaces!(cplx, k, τ, N, E)
    # V ← V ∪ {τ}
    τdim = dim(τ)
    cplx[τ, τdim] > size(cplx, τdim) && push!(cplx, τ)
    # stop recursion
    dim(τ) >= k && return
    τvals = τ[:values]
    for v in N
        first(v[:values]) in τvals && continue
        # σ ← τ ∪ {v}
        σ = Simplex(τvals..., v[:values]...)
        M = intersect(N, lowernbrs(cplx, v, E))
        addcofaces!(cplx, k, σ, M, E)
    end
end

"""Incremental construction of VR complex from neighborhood graph"""
function incremental!(cplx, k, E)
    # construct VR complex
    for u in cells(cplx, 0)
        N = lowernbrs(cplx, u, E)
        addcofaces!(cplx, k, u, N, E)
    end
end

"""Weight function for VR complex"""
function weight(σ, w, cplx)
    k = dim(σ)
    k == 0 && return 0.
    k == 1 && return w[1][cplx[σ, 1]]
    return maximum([w[k-1][cplx[τ, k-1]] for τ in faces(σ)])
end

function expand(method, cplx, w, kmax, E)
    # perform nerve expansion
    if method == :incremental
        incremental!(cplx, kmax, E)
    elseif method == :inductive
        inductive!(cplx, kmax, E)
    else
        throw(ArgumentError("Invalid method name $(method)"))
    end

    if w !== nothing
        for k in 2:dim(cplx)
            w[k] = zeros(size(cplx,k))
            for σ in cells(cplx,k)
                w[k][σ[:index]] = weight(σ, w, cplx)
            end
        end
    end
    return
end

"""
    vietorisrips(X, ɛ[, weights])

Construction of a Vietoris–Rips complex

# Arguments
- `X::AbstractMatrix{T}`: the dataset where each column defines a data point.
- `ɛ::Real`: the maximal scale value.
- `weights::Bool`: this flag indicates if the filteration weight should be calculated. Default value is `true`.

This is an inefficient implementation (based on distance matrix) of the VRC algorithm from:

    Afra Zomorodian, Fast construction of the Vietoris-Rips complex, 2010
"""
function vietorisrips(X::AbstractMatrix{T}, ɛ::Real, weights = true;
                      expansion = :incremental,
                      distance  = Distances.Euclidean(),
                      maxoutdim = size(X,2)-1) where T <: Real

    # calculate cross distances for each point in dataset
    D = Distances.pairwise(distance, X)

    # determine maximal dimension
    kmax = min(maxoutdim, maximum(mapslices(c->count(d-> 0.0 < d ≤ ɛ, c), D, dims=1)))

	return vrfromdistances(D, ϵ, weights, expansion, maxoutdim)
end

function vrfromdistances(D::AbstractMatrix{T}, ɛ::Real, weights = true;
                      expansion = :incremental,
                      maxoutdim = size(D,1)-1) where T <: Real
	n = size(D,1)

	# add zero-dimensional simplices to complex
	splxs = map(Simplex, 1:n)
	cplx = SimplicialComplex(splxs...)

	# determine maximal dimension
	kmax = min(maxoutdim, maximum(mapslices(c->count(d-> 0.0 < d ≤ ɛ, c), D, dims=1)))

	# build 1-skeleton (neighborhood graph)
	E = spzeros(Bool, n, n) # adjacency matrix
	for σ in cells(cplx, 0)
		u = σ[:values]
		Du = view(D, :, u)
		edges = (LinearIndices(Du))[findall(d->d <= ɛ && d > 0., Du)]
		length(edges) == 0 && continue
		for i in edges
			v = cplx[i, 0][:values]
			s = Simplex(u...,v...)
			if cplx[s] > size(cplx, 1)
				# add simplex to complex
				push!(cplx, s)
				# fill adjacency matrix
				E[s[:values], s[:values]] .= true
			end
		end
	end

	# calculate weights of nerve
	w = nothing
	if weights
		w = Dict{Int, Vector{T}}()
		w[0] = zeros(size(cplx,0))
		if size(cplx, 1) > 0
			w[1] = zeros(size(cplx,1))
			for e in cells(cplx,1)
				w[1][e[:index]] = D[e[:values]...]
			end
		end
	end

	# perform expansion
	expand(expansion, cplx, w, kmax, E)

	return cplx, w
end

"""
Landmark selection for witness complex
"""
function landmarks(X::AbstractMatrix, l::Int; method = :minmax,
                   distance = Distances.Euclidean(), firstpoint = 0)
    d, n = size(X)
    if method == :random
        return Random.randperm(n)[1:l]
    elseif method == :minmax
        idxs = firstpoint > 0 ? [firstpoint] : rand(1:n, 1)
        Z = view(X, :, idxs)
        for i in 2:l
            maxv = -Inf
            maxvidx = 0
            for j in 1:n
                if j ∉ idxs
                    v, _ = findmin(Distances.colwise(distance, X[:,j], Z))
					if v > maxv
					    maxv = v
                        maxvidx = j
                    end
                end
            end
            push!(idxs, maxvidx)
        end
        return idxs
    else
        throw(ArgumentError("Invalid method name $(method)"))
    end
end

"""
    witness(X, l, ɛ[, weights])

Construction of a witness complex.

# Arguments
- `X::AbstractMatrix{T}`: the dataset where each column defines a data point.
- `l::Int`: the number of landmarks.
- `ɛ::Real`: the maximal scale value.
- `weights::Bool`: this flag indicates if the filteration weight should be calculated. Default value is `true`.

# Keyword arguments
For parameter ν = 0, 1, 2 determines size of landmark radius.
- If ν = 0, then for i = 1, 2, . . . , N define m_i = 0.
- If ν > 0, then for i = 1, 2, . . . , N define m_i to be the ν-th smallest entry of the i-th column of D.
"""
function witness(X::AbstractMatrix{T}, l::Int, ɛ::Real, weights = true;
                 landmark = :minmax, distance = Distances.Euclidean(),
                 expansion = :incremental, ν::Int = 2, maxoutdim = size(X,2)-1,
                 firstpoint = 0) where T <: Real

    # get landmarks
    L = landmarks(X, l, method = landmark, distance = distance, firstpoint=firstpoint)

    # get distances to landmarks
    D = Distances.pairwise(distance, X[:,L], X)

    # simplexes contain indexes to landamarks
    cplx, w = witness(D, ɛ, weights, expansion=expansion, ν=ν, maxoutdim=maxoutdim, firstpoint=firstpoint)

    return cplx, w
end

"""
    witness(D, ɛ[, weights])

Construction of a witness complex from the distance matrix `D` with a maximum scale value `ɛ`.
"""
function witness(D::AbstractMatrix, ɛ::Real, weights = true;
                 expansion = :incremental,
                 ν::Int = 2,
                 maxoutdim = size(D,1)-1,
                 firstpoint = 0) where T <: Real

    @assert ν < 3 "ν only can take values of 0, 1 or 2"
    @assert maxoutdim > 0 "Maximum dimension should be more then 0"

    l, n = size(D)

    # get landmark radius extension m_i
    m = zeros(eltype(D), n)
    if ν > 0
        idxs = zeros(Int, l)
        for i in 1:n
            sortperm!(idxs, D[:,i])
            m[i] = D[idxs[ν], i]
        end
    end

    # construct 1-seleton edges
    splxs = map(Simplex, 1:l)
    cplx = SimplicialComplex(splxs...)
    w = nothing
    if weights
        w = Dict{Int, Vector{Float64}}()
        w[0] = zeros(size(cplx,0))
        w[1] = zeros(0)
    end

    # the edge σ = [ab] belongs to W(D; ɛ, ν) iff there exists a witness i ∈ {1, 2, ..., N} such that:
    # max(D(a, i), D(b, i)) ≤ ɛ + m_i
    E = spzeros(Bool, l, l) # adjacency matrix
    for i in 1:l
        for j in i+1:l
            e = Inf
            @inbounds for k in 1:n
                d = max(D[i, k], D[j, k])
                d = d < m[k] ? 0.0 : d - m[k]
                if d < e
                    e = d
                end
            end

            if e ≤ ɛ
                push!(cplx, Simplex(i, j))
                weights && push!(w[1], e)
                E[i, j] = true
            end
        end
    end

    # determine maximal dimension
    Rmax = ɛ+maximum(m)
    kmax = min(maxoutdim, maximum(mapslices(c->count(d->0.0<d≤Rmax, c), D, dims=1)))

    # perform expansion
    expand(expansion, cplx, w, kmax, E)

    return cplx, w
end


"""
    welzl(P, R)

Computes the smallest enclosing circle of the union of `P` and `R`.
"""
function welzl(P, R)
    p = length(P)
    r = length(R)
    if p == 0 || r >= 3
        if r == 1
            return R[1], 0.0
        elseif r == 2
            p1 = R[1]
            p2 = R[2]
            m = (p1.+p2) ./ 2.
            d = sqrt(sum(abs2, (x .- y)))
            return m, d / 2.
        elseif true
            @warn "the points of R are cocircular"
            return [], 0.0
        else
            @error "Undefined"
            return [], -1.0
       end
    end
    i = rand(1:p)
    x = P[i]
    c, r = welzl(deleteat!(P, i), R)
    norm(x - c) <= r && return (c, r)
    return welzl(P, push!(R, x))
end

"""
    čech(X, ɛ[, weights])
    cech(X, ɛ[, weights])

Construction of a Čeck complex.

# Arguments
- `X::AbstractMatrix{T}`: the dataset where each column defines a data point.
- `ɛ::Real`: the maximal scale value.
- `weights::Bool`: this flag indicates if the filteration weight should be calculated. Default value is `true`.

"""
function čech(X::AbstractMatrix{T}, ɛ::Real, weights = true;
              distance  = Distances.Euclidean(),
              maxoutdim = size(X,2)-1) where T <: Real

    d, n = size(X)

    splxs = map(Simplex, 1:n) # construct 0-simplices
    cplx = SimplicialComplex(splxs...) # and 0-complex

    w = nothing
    if weights
        w = Dict{Int, Vector{Float64}}()
        w[0] = zeros(size(cplx,0))
    end

    for h in 1:maxoutdim
        idx = collect(h+1:-1:1)
        weights && setindex!(w, zeros(T, 0), h)

        while idx[end] <= n-h
            # find radius of bounding sphere
            c, r = BoundingSphere.boundingsphere([X[:, i] for i in idx])

            # add simplex in radius less then maximal filtration value
            if r <= ε
                s = Simplex(reverse(idx))
                addsimplex(cplx, s)
                weights && push!(w[h], r)
            end

            # iterate through all combinations of point indecies
            for i in 1:h+1
                idx[i] += 1
                if idx[i] <= n-i+1
                    for j in i-1:-1:1
                        idx[j] = idx[j+1] + 1
                    end
                    break
                end
            end
        end
    end

    return cplx, w
end
const cech = čech
