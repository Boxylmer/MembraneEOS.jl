# using MembraneEOS
using Test

@testset "MembraneEOS.jl" begin
    include("TestChemicalLookup.jl")
    # include("TestPR.jl")
    include("TestSL.jl")
end
