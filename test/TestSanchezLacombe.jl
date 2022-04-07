@testset "TestSanchezLacombe.jl" begin


    model_co2 = MembraneEOS.SL([630.0, 630.0], [300.0, 300.0], [1.515, 1.515], [44.0, 44.0], [0 0; 0 0])
    model_co2_pure = MembraneEOS.SL([630.0], [300.0], [1.515], [44.0])

    model_pdms = MembraneEOS.SL([302.0], [476.0], [1.104], [1.00E+30])
    @show MembraneEOS.mass_density(model_pdms, 1, 273.15)
    @show MembraneEOS.chemical_potential(model_pdms, 1, 273.15)
    @show MembraneEOS.chemical_potential(model_co2, 1, 273.15, [0.0, 0.2])
    @show MembraneEOS.chemical_potential(model_co2_pure, 1, 273.15)
end

function test_sanchez_lacombe_pure_gas()
    O2 = ChemicalParameters(molecular_weight=32, characteristic_temperature=170.0, characteristic_pressure=280.0, characteristic_density=1.29,
        critical_pressure=49.8, critical_temperature=154.6, acentric_factor=0.02)
    CO2 = ChemicalParameters(molecular_weight=44, characteristic_temperature=300.0, characteristic_pressure=630.0, characteristic_density=1.515,
        critical_pressure=72.8, critical_temperature=304.19, acentric_factor=0.228)
    p_mpa = 0.101325
    p_atm = p_mpa / 0.101325
    t_k = 200.15
    # state_O2 = SanchezLacombeState([O2], p_mpa, t_k, [1.0], [0.][:,:])
    # prstate_O2 = PengRobinsonState(O2; p=p_atm, t=t_k)

    # state_CO2 = SanchezLacombeState([CO2], p_mpa, t_k, [1.0], [0.][:,:])
    # prstate_CO2 = PengRobinsonState(CO2; p=p_atm, t=t_k)
    # @show volume(state_O2), volume(state_CO2)
    # @show volume(prstate_O2), volume(prstate_CO2)
    
    # mixed gas O2 and CO2
    mole_fractions = [0.5, 0.5]
    mass_fractions = PolymerMembranes.mole_fractions_to_mass_fractions(
        mole_fractions, PolymerMembranes.molecular_weight.([O2, CO2]))
    state_O2_CO2 = SanchezLacombeState([O2, CO2], p_mpa, t_k, mass_fractions, [0. 0.; 0. 0.])
    prstate_O2_CO2 = PengRobinsonState([O2, CO2], mole_fractions, [0. 0.; 0. 0.]; p=p_atm, t=t_k)

    # make sure both interfaces work
    state_O2_CO2_standard_interface = SanchezLacombeState([O2, CO2], mole_fractions, [0. 0.; 0. 0.]; p=p_atm, t=t_k)
    state_CO2 = SanchezLacombeState([CO2], p_mpa, t_k, [1.0], [0.][:,:])
    state_CO2_monocomponent_interface = SanchezLacombeState(CO2, p_mpa, t_k)
    @test chemical_potentials(state_CO2) == chemical_potentials(state_CO2_monocomponent_interface) 
    # @show volume(state_O2_CO2)
    # @show volume(state_O2_CO2_standard_interface)
    # @show volume(prstate_O2_CO2)
    # @test volume(state_O2_CO2_standard_interface) == volume(state_O2_CO2)

    # make sure pressures are recoverable (the pressure equation and density equation match up)
    v_l_mol = volume(state_O2_CO2)
    state_O2_CO2_p = SanchezLacombeState([O2, CO2], mole_fractions, [0. 0.; 0. 0.]; v=v_l_mol, t=t_k)
    @test pressure(state_O2_CO2_p) ≈ p_atm * 0.101325
    
    # @show PolymerMembranes.compute_partial_pressures(prstate_O2_CO2)
    # @show μ_O2_CO2 = state_O2_CO2.chemical_potentials

    
end
function test_sanchez_lacombe_methods()
    O2 = ChemicalParameters(molecular_weight=32, characteristic_temperature=170.0, characteristic_pressure=280.0, characteristic_density=1.29,
        critical_pressure=49.8, critical_temperature=154.6, acentric_factor=0.02)
    CO2 = ChemicalParameters(molecular_weight=44, characteristic_temperature=300.0, characteristic_pressure=630.0, characteristic_density=1.515,
        critical_pressure=72.8, critical_temperature=304.19, acentric_factor=0.228)
    kijtable = [0 0.01; 0.01 0]
    model = SanchezLacombeModel([O2, CO2], kijtable)
    p_mpa = 0.101325
    p_atm = p_mpa / 0.101325
    t_k = 200.15
    mole_fractions = [0.5, 0.5]

    myvol = volume(model, p_mpa, t_k, mole_fractions)
    @test p_mpa ≈ pressure(model, myvol, t_k, mole_fractions)
    
    # mixed gas O2 and CO2
    mass_fractions = PolymerMembranes.mole_fractions_to_mass_fractions(
        mole_fractions, PolymerMembranes.molecular_weight.([O2, CO2]))
    density = PolymerMembranes.molar_volume_to_density(myvol, mole_fractions, PolymerMembranes.molecular_weight.([O2, CO2]))
    μ_pt = chemical_potential(model, p_mpa, t_k, mole_fractions)
    μ_ρt = ρTω_chemical_potential(model, density, t_k, mass_fractions)
    @test μ_pt ≈ μ_ρt

    
    
end
function test_sanchez_lacombe()
    @testset "TestSanchezLacombe.jl" begin
        test_sanchez_lacombe_pure_gas()
        test_sanchez_lacombe_methods()
    end
end
