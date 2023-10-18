using Test
using MembraneBase
using MembraneEOS
using Measurements

@testset "MembraneEOS.jl" begin
    include("TestChemicalLookup.jl")
    include("TestPR.jl")
    include("TestSL.jl")


    @testset "EOS Abstractions" begin
    # make sure all abstractions throw when called alone

        struct testmodel <: MembraneEOS.MEOSModel end
        m = testmodel()

        @test_throws ErrorException pressure(m, nothing, nothing, nothing)
        @test_throws ErrorException volume(m, nothing, nothing, nothing)
        @test_throws ErrorException fugacity(m, nothing, nothing, nothing)
        @test_throws ErrorException molecular_weight(m)
        @test_throws ErrorException compressibility_factor(m, nothing, nothing, nothing)
        @test_throws ErrorException VT_compressibility_factor(m, nothing, nothing, nothing) 
        @test_throws ErrorException mass_density(m, nothing) 
        @test_throws ErrorException VT_mass_density(m, nothing)
        @test_throws ErrorException chemical_potential(m, nothing) 
        @test_throws ErrorException chemical_potential_res(m, nothing)
        @test_throws ErrorException VT_chemical_potential(m, nothing)
        @test_throws ErrorException ρTω_chemical_potential(m, nothing)
        @test_throws ErrorException ρTω_chemical_potential_res(m, nothing)
        @test_throws ErrorException activity(m, nothing) 
        @test_throws ErrorException ρTω_activity(m, nothing) 
        @test_throws ErrorException density_upper_bound(m, nothing) 
        @test_throws ErrorException get_kij(m, "1", "2") 
        @test_throws ErrorException get_kij_matrix(m, ["1", "2"])
    end
end

nothing
