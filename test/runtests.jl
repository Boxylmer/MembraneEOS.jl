using MembraneEOS
using Test

@testset "MembraneEOS.jl" begin
    include("TestChemicalLookup.jl")
    # include("TestPengRobinson.jl")
    # include("TestSanchezLacombe.jl")
end
