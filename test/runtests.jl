using Test
using MembraneBase
using MembraneEOS
using Measurements

@testset "MembraneEOS.jl" begin
    include("TestChemicalLookup.jl")
    include("TestPR.jl")
    include("TestSL.jl")
end

nothing
