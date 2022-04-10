@testset "TestSanchezLacombe.jl" begin

    test_lookup_model = SL(["CO2", "CH4"], [0 0; 0 0])

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
    z_pt = compressibility_factor(model_co2_pure, 1, 273.15)
    z_vt = VT_compressibility_factor(model_co2_pure, v, 273.15)
    @test z_pt ≈ z_vt

    μ_pt = chemical_potential(model_co2_pure, 1, 273.15)
    μ_vt = VT_chemical_potential(model_co2_pure, v, 273.15)
    @test μ_pt ≈ μ_vt

    co2_ch4_model = SL(["CO2", "CH4"])
    z = [0.2, 0.8]
    p = 1; t = 273.15
    ω = mole_fractions_to_mass_fractions(z, molecular_weight(co2_ch4_model))
    ρ = mass_density(co2_ch4_model, p, t, z)  # g/cm3
    v = volume(co2_ch4_model, p, t, z)  # l/mol
    μ_ρtω = ρTω_chemical_potential(co2_ch4_model, ρ, t, ω)
    μ = VT_chemical_potential(co2_ch4_model, v, t, z)
    @test μ ≈ μ_ρtω

    @show μ_ρtω = ρTω_chemical_potential(co2_ch4_model, ρ, t, ω)
    @show μ_ρtω_res = ρTω_chemical_potential_res(co2_ch4_model, ρ, t, ω)
    @show a = ρTω_activity(co2_ch4_model, ρ, t, ω)
    @show a_res = ρTω_activity_res(co2_ch4_model, ρ, t, ω)

end
