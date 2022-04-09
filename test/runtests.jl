# using MembraneEOS
using Test
using MembraneEOS
using Random
using Measurements
using MembraneBase

@testset "MembraneEOS.jl" begin
    include("TestChemicalLookup.jl")
    include("TestPR.jl")
    include("TestSL.jl")
end
nothing