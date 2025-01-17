"""Filtration of an abstract complex

We call this sequence of complexes the **filtration** of `f` and
think of it as a construction by adding chunks of simplices at a time `t::FI`.
∅ = K0 ⊆ K1 ⊆ . . . ⊆ Kn = K.
"""
mutable struct Filtration{C<:AbstractComplex, FI}
    # underlying abstract cell complex
    complex::C
    # total order of simplices as array of (dim, simplex id, filtation value)
    total::Vector{Tuple{Int,Int,FI}}
    divisions::Number
end
Base.show(io::IO, flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = print(io, "Filtration($(complex(flt)), $FI)")
Base.valtype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = FI

Base.complex(flt::Filtration) = flt.complex
order(flt::Filtration) = flt.total
Base.minimum(flt::Filtration) = order(flt)[1][3]
Base.maximum(flt::Filtration) = order(flt)[end][3]
#
# Constructors
#
Filtration(::Type{C}, ::Type{FI}) where {C <: AbstractComplex, FI} =
    Filtration(C(), Vector{Tuple{Int,Int,FI}}(), Inf)
Base.similar(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = Filtration(C, FI)

"""Construct filtration from a cell complex using the order of their appearence in the complex"""
function filtration(cplx::AbstractComplex)
    idx = Vector{Tuple{Int,Int,Int}}()
    i = 1
    for d in 0:dim(cplx)
        for c in cells(cplx, d)
            push!(idx, (dim(c), c[:index], i))
            i += 1
        end
    end
    return Filtration(cplx, idx, Inf)
end

"""Construct filtration from a cell complex and a complex weight function"""
function filtration(cplx::C, w::Dict{Int,Vector{FI}}; divisions::Number=Inf) where {C<:AbstractComplex, FI}
    idx = Vector{Tuple{Int,Int,FI}}()
    for d in 0:dim(cplx)
        for c in cells(cplx, d)
            ci = c[:index]
            push!(idx, (d, ci, w[d][ci]))
        end
    end
    sort!(idx, by=x->(x[3], x[1])) # sort by dimension & filtration value
    return Filtration(cplx, idx, divisions)
end

function Base.push!(flt::Filtration{C,FI}, cl::AbstractCell, v::FI; recursive=false) where {C<:AbstractComplex, FI}
    cplx = complex(flt)
    @assert isa(cl, celltype(cplx)) "Complex $(cplx) does not accept $(typeof(cl))"
    cls = push!(cplx, cl, recursive=recursive)
    ord = order(flt)
    idx = length(ord) == 0 ? 1 : findlast(e->e[3]<=v, ord)
    for c in sort!(cls, by=s->dim(s))
        if idx == length(ord)
            push!(ord, (dim(c), c[:index], v))
        else
            insert!(ord, idx, (dim(c), c[:index], v))
        end
        idx += 1
    end
    return flt
end

"""Generate a combined boundary matrix from the filtration `flt` for the persistent homology calculations."""
function boundary_matrix(flt::Filtration; reduced=false)
    ridx = reduced ? 1 : 0
    # initialize boundary matrix
    cplx = complex(flt)
    bm = map(i->BitSet(), 1:sum(size(cplx))+ridx)
    # fill boundary matrix
    ord = order(flt)
    for (i, (d, ci, fv)) in enumerate(ord)
        if d > 0
            splx = cplx[ci, d]
            for face in faces(splx)
                fi = cplx[face, d-1]
                push!(bm[i+ridx], findfirst(e->e[1] == d-1 && e[2] == fi, ord)+ridx)
            end
        elseif reduced
            push!(bm[i+ridx], 1)
        end
    end
    return bm
end

function SparseArrays.sparse(∂::Vector{BitSet})
    m = length(∂)
    ret = spzeros(Int, m, m)
    for i in 1:m
        bm = ∂[i]
        for (l, j) in enumerate(bm)
            ret[j,i] = j # (-1)^((l-1)%2) # coefs require exact order of faces in provides simplex
        end
    end
    return ret
end

#
# I/O
#

function Base.write(io::IO, flt::Filtration)
    cplx = complex(flt)
    for (d, ci, fv) in order(flt)
        for k in cplx[ci,d][:values]
            write(io, "$k,")
        end
        write(io, "$fv\n")
    end
end

function Base.read(io::IO, ::Type{Filtration{C,FI}}) where {C <: AbstractComplex, FI}
    flt = Filtration(C,FI)
    ET = eltype(celltype(complex(flt)))
    while !eof(io)
        l = readline(io)
        vals = split(l, ',')
        svals = map(v->parse(ET, v), vals[1:end-1])
        fval = parse(FI, vals[end])
        push!(flt, Simplex(svals), fval)
    end
    return flt
end

"Write a combined boundary matrix to a text file"
function writeboundarymatrix(io::IO, bm::Vector, zeroindex = true)
    for smplx in bm
        if length(smplx) == 0
            write(io, "0")
        else
            write(io, "$(length(smplx)-1)")
        end
        for i in smplx
            write(io, " $(zeroindex ? i-1 : i)")
        end
        write(io, 0x0A)
    end
end

#
# Iterators
#

# Filtration simplicial complex iterator
Base.length(flt::Filtration) = isinf(flt.divisions) ? length(unique(e->e[3], order(flt))) : flt.divisions
Base.eltype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = C

"""Loop through the filtration `flt` producing growing simplicial complexes on every iteration"""
function Base.iterate(flt::Filtration, state=(1, -Inf, 0))
    ord = order(flt)
    idx = state[1]
    idx > length(ord) && return nothing # done
    if state[2] == -Inf # calculate initial state
        fval = ord[idx][3]
        incr = (maximum(flt)-minimum(flt)) / flt.divisions
    else
        fval = state[2]
        incr = state[3]
    end
    splxs = Tuple{Int,Int}[] #simplex dim & index
    while idx <= length(ord) && (fval+incr) >= ord[idx][3]
        push!(splxs, ord[idx][1:2])
        idx += 1
    end
    nextfval = fval+incr
    if idx <= length(ord) && isinf(flt.divisions)
        nextfval = ord[idx][3]
    end
    return (fval, splxs), (idx, nextfval, incr)
end

# Filtration simplex iterator
Base.length(splxs::Simplices{<:Filtration}) = length(splxs.itr)
Base.eltype(splxs::Simplices{<:Filtration}) = celltype(splxs.itr)

function Base.iterate(splxs::Simplices{<:Filtration},state=nothing)
    # call underlying iterator
    if state === nothing
        res = iterate(splxs.itr)
    else
        res = iterate(splxs.itr, state)
    end
    # final state
    res == nothing && return nothing
    # get complex
    cplx = complex(splxs.itr)
    ss = [cplx[i, d] for (d, i) in res[1][2]]
    return (res[1][1], ss), res[2] # state of filtration iterator
end
