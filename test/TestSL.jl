@testset "TestSanchezLacombe.jl" begin


    
    model_pdms = MembraneEOS.SL([302.0], [476.0], [1.104], [1.00E+30])
    # @show MembraneEOS.mass_density(model_pdms, 1, 273.15)
    # @show MembraneEOS.chemical_potential(model_pdms, 1, 273.15)
    # @show MembraneEOS.chemical_potential(model_co2, 1, 273.15, [0.0, 0.2])
    # @show MembraneEOS.chemical_potential(model_co2_pure, 1, 273.15)
    
    model_co2_pure = MembraneEOS.SL([630.0], [300.0], [1.515], [44.0])
    vol = MembraneEOS.volume(model_co2_pure, 1., 273.15)
    pres = MembraneEOS.pressure(model_co2_pure, vol, 273.15)
    @test 1. ≈ pres
    
    model_co2 = MembraneEOS.SL([630.0, 630.0], [300.0, 300.0], [1.515, 1.515], [44.0, 44.0], [0 0; 0 0])
    vol = MembraneEOS.volume(model_co2, 1., 273.15, [0.5, 0.5])
    pres = MembraneEOS.pressure(model_co2, vol, 273.15, [0.5, 0.5])
    @test 1. ≈ pres
    
end
