module ComputationalHomology

import SmithNormalForm
import Distances

import Base: ==, -, +, *

export AbstractCell,
       dim, faces, volume,
       Simplex, Cube,

       AbstractChain,
       setdim!, simplify,
       Chain,

       AbstractComplex,
       boundary, coboundary, celltype, cells, boundary_matrix,
       SimplicialComplex,
       addsimplex, addsimplex!,
       vietorisrips,

       AbstractHomology,
       grouptype, group,
       Homology,
       homology, withgenerators, generators

global SNF = SmithNormalForm.snf

include("cells.jl")
include("chain.jl")
include("complex.jl")
include("simplicial.jl")
include("vietorisrips.jl")
include("homology.jl")

function setsnf!(f::Function)
    global SNF
    SNF = f
end

end # module