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

    μ_pure = chemical_potential(model_co2_pure, 1, 273.15)[1]
    μ_psuedo_pure = chemical_potential(model_co2, 1, 273.15, [1, 0])[1]
    @test μ_pure ≈ μ_psuedo_pure

    fug_co2 = fugacity(model_co2_pure, 1, 273.15)
    
    v = volume(model_co2_pure, 1, 273.15)
    @show z_pt = compressibility_factor(model_co2_pure, 1, 273.15)
    @show z_vt = VT_compressibility_factor(model_co2_pure, v, 273.15)
    @test z_pt ≈ z_vt

    @show μ_pt = chemical_potential(model_co2_pure, 1, 273.15)
    @show μ_vt = VT_chemical_potential(model_co2_pure, v, 273.15)


end
