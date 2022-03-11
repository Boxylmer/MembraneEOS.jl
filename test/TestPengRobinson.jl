function test_preos_pvt_solvers(;precision=10, speedtests=true)
    # nitrogen = PengRobinsonChemical(33.5, 126.2, 0.04, "N2")
    # methane = PengRobinsonChemical(45.79, 190.8, 0.012, "CH4")
    
    nitrogen = ChemicalParameters("N2", critical_pressure=33.5, critical_temperature=126.2, acentric_factor=0.04, fill_missing_with_database=false)
    methane = ChemicalParameters("CH4", critical_pressure=45.79, critical_temperature=190.8, acentric_factor=0.012,  fill_missing_with_database=false)
    chemical_components = [nitrogen, methane]
    mole_fractions = [0.4, 0.6]
    kij_table = initmatrix(chemical_components)
    kij_table[1, 2] = 0.03
    kij_table[2, 1] = 0.03

    # single components
    @test isapprox(PengRobinsonState([nitrogen], [1.], [0.][:,:]; v=20.08405, t=273.15).pressure, 1.1150630860000421)  
    # multiple components
    @test isapprox(PengRobinsonState(chemical_components, [0.4, 0.6], kij_table; v=20.08405, t=273.15).pressure, 1.1136140330622173) 
    # solving volume (mixtures)
    solved_vol = round(PengRobinsonState(chemical_components, [0.4, 0.6], kij_table; p=1.113612752274884, t=273.15).volume; digits=precision)
    target_vol = round(20.08407315; digits=precision)
    @test solved_vol == target_vol
    
    # solving fugacity 
    state = PengRobinsonState(chemical_components, [0.4, 0.6], kij_table; p=1.113612752274884, t=273.15)
    @test fugacities(state) == [0.44517995494867035, 0.6660438011048816]
    ammonia = ChemicalParameters(; critical_pressure=11.357, critical_temperature=405.5, acentric_factor=0.253)
    state_ammonia_1 = PengRobinsonState(ammonia; p=10 * 9.86923, t=473)
    state_ch4 = EOS(PengRobinson(), ChemicalParameters("CH4"); p=30, t=293.15)
    state_ch4_n2 = EOS(PengRobinson(), ["CH4", "N2"], [0.4, 0.6]; p=39.4769, t=293.15)
    
    # hydrogen = CubicModel()

    if speedtests
        # test volume solver time
        time_volume_solver =  @belapsed PengRobinsonState($chemical_components, [0.4, 0.6], $kij_table; p=1.113612752274884, t=273.15).volume
        println("Time (μs) to calculate multi-component PengRobinson volume: " * string(time_volume_solver * 10e6) * " (target: " * string(100) * ")")
        
        # ensure base calculations can be done in less than 10 microseconds
        target_time = 10e-6
        time_single_comp = @belapsed PengRobinsonState([$nitrogen], [1.], [0.][:,:]; v=20.08405, t=273.15)
        time_multi_comp  = @belapsed PengRobinsonState($chemical_components, [0.4, 0.6], $kij_table; v=20.08405, t=273.15)
        println("Time (μs) to calculate single-component PengRobinson pressure: " * string(time_single_comp * 10e6) * " (target: " * string(round(target_time * 10e6; digits=3)) * ")")
        println("Time (μs) to calculate multi-component PengRobinson pressure: " * string(time_multi_comp * 10e6) * " (target: " * string(round(target_time * 10e6; digits=3)) * ")")
        @test time_single_comp < target_time
        @test time_multi_comp < target_time
    end
end
function test_preos_pvt_solvers_methods(;precision=10, speedtests=true)
    # nitrogen = PengRobinsonChemical(33.5, 126.2, 0.04, "N2")
    # methane = PengRobinsonChemical(45.79, 190.8, 0.012, "CH4")
    
    nitrogen = ChemicalParameters("N2", critical_pressure=33.5, critical_temperature=126.2, acentric_factor=0.04, fill_missing_with_database=false)
    methane = ChemicalParameters("CH4", critical_pressure=45.79, critical_temperature=190.8, acentric_factor=0.012,  fill_missing_with_database=false)
    chemical_components = [nitrogen, methane]
    mole_fractions = [0.4, 0.6]
    kij_table = initmatrix(chemical_components)
    kij_table[1, 2] = 0.03
    kij_table[2, 1] = 0.03
    nitrogenmodel = CubicModel(PengRobinson(), [nitrogen])

    model = CubicModel(PengRobinson(), chemical_components, kij_table)

    # single component pressure
    @test isapprox(pressure(nitrogenmodel, 20.08405, 273.15), 1.1150630860000421)  
    # multiple component pressure
    @test isapprox(pressure(model, 20.08405, 273.15, mole_fractions), 1.1136140330622173) 
    
    # multicomponent volume
    solved_vol = round(volume(model, 1.113612752274884, 273.15, mole_fractions); digits=precision)
    target_vol = round(20.08407315; digits=precision)
    @test solved_vol == target_vol
    
    # single component fugacity
    ammonia = ChemicalParameters(; critical_pressure=11.357, critical_temperature=405.5, acentric_factor=0.253)
    ammoniamodel = CubicModel(PengRobinson(), [ammonia])
    @test fugacity(ammoniamodel, 1, 273.15)[1] ≈ 0.3890265134042765
    # mutlicomponent fugacity
    @test fugacity(model, 1.113612752274884, 273.15, mole_fractions) == [0.44517995494867035, 0.6660438011048816]
    
    # compressibility factor methods
    my_volume = volume(model, 1, 273.15, mole_fractions)
    @test compressibility_factor(model, 1, 273.15, mole_fractions) ≈ VT_compressibility_factor(model, my_volume, 273.15, mole_fractions)
    
end

function test_preos_measurements_integration(;precision=5)
    # define some working inputs
    nitrogen_with_err = ChemicalParameters("N2", critical_pressure=33.5 ± 1, critical_temperature=126.2 ± 1, acentric_factor=0.04 ± 0.005, fill_missing_with_database=false)
    methane_with_err = ChemicalParameters("CH4", critical_pressure=45.79 ± 10, critical_temperature=190.8 ± 1, acentric_factor=0.012 ± 0.003, fill_missing_with_database=false)
    chem_comps_with_err = [nitrogen_with_err, methane_with_err]
    kij_table_with_err = [(0. ± 0.) (0.03 ± 0.01); (0.03 ± 0.01) (0. ± 0.)]
    
    pres_1 = PengRobinsonState(chem_comps_with_err, [0.4 ± 0.1, 0.6 ± 0.1], kij_table_with_err; v=20.08405 ± 1.2, t=273.15 ± 1).pressure
    @test round(pres_1.err; digits=precision) == round(0.06652796; digits=precision)
    
    vol_1 = PengRobinsonState(chem_comps_with_err, [0.4 ± 0.1234567, 0.6 ± 0.1234567], kij_table_with_err; p=(1.113612 ± 0.1234567), t=273.15 ± 1.3).volume
    vol_1_val = round(vol_1.val; digits=precision)
    vol_1_err = round(vol_1.err; digits=precision)
    vol_1_err_ans = round(2.23349794; digits=precision)
    @test (vol_1_err == vol_1_err_ans)
    
    vol_2 = PengRobinsonState(chem_comps_with_err, [0.4 ± 0.1234567, 0.6 ± 0.1234567], kij_table_with_err; p=(1.113612 ± 0.1234567), t=273.15 ± 10.3).volume # should be 20.08405
    vol_2_val = round(vol_2.val; digits=precision)
    vol_2_err = round(vol_2.err; digits=precision)
    vol_2_err_ans = round(2.35830357; digits=precision)
    @test vol_2_err == vol_2_err_ans
    
    @test vol_1_val == vol_2_val
    @test vol_2_val == round(20.08408675; digits=precision)
end

function test_preos_measurements_validity(; precision=2)
    """
    From @CaseyKneale
    So I did a little thinking on this... 
    What I would do is, make a bunch of artificial P's and +/- P's using an EOS with a variety of coefficients (maybe 10k random examples). 
    Then back track to solve for V and V+/- using your approach. 
        Tabulate the error in your approach, specifically on the uncertainty calculations. 
        I think what you'll find is that what you're doing is OK. But if not, we'll do a deep dive.
        in theory I think what you're doing is fine, but put it into practice in a verifiable way. 
        Should probably do that anyways from the perspective of unit tests.
        
        
        VoidEyes — Yesterday at 4:50 PM
        i mean, if you have an starting point, (V0+V0err) and recalculate pressure, the error is obviously higher that the errors expected by V and T alone (the parameters). 
        but as a way to aproaching the real V+Verr, you could use P0err - Perr ≈ Param_err
        Kirory — Yesterday at 4:51 PM
        That's what's done, I think, save for that its squared
        VoidEyes — Yesterday at 4:51 PM
        yeah, thats ok
        
        """
        
    num_iters = 10000 # see above conversations
    Random.seed!(1234)  # repeatability
    error_array = zeros(num_iters)
    
    
    for idx in 1:num_iters
        # define some working inputs for a realistic system
        nitrogen_with_err = ChemicalParameters("N2", critical_pressure=randomize(33.5) ± randomize(1), critical_temperature=randomize(126.2) ± randomize(1), acentric_factor=randomize(0.04) ± randomize(0.005), fill_missing_with_database=false)
        methane_with_err = ChemicalParameters("CH4", critical_pressure=randomize(45.79) ± randomize(10), critical_temperature=randomize(190.8) ± randomize(1), acentric_factor=randomize(0.012) ± randomize(0.003), fill_missing_with_database=false)
        chem_comps_with_err = [nitrogen_with_err, methane_with_err]
        kij = randomize(0.03) ± randomize(0.01)
        kij_table_with_err = [(0. ± 0.) (kij); kij (0. ± 0.)]
        mfrac_1 = randomize(0.04) ± randomize(0.1)
        mfrac_2 = (1 - mfrac_1.val) ± randomize(0.1)
        mole_fractions = [mfrac_1, mfrac_2]
        known_v = randomize(20.08405) ± randomize(1.2)
        known_t = randomize(273.15) ± randomize(1)
        
        # given this system, calculate the pressure given V and T (the "trivial" case)
        known_p = PengRobinsonState(chem_comps_with_err, mole_fractions, kij_table_with_err; v=known_v, t=known_t, minimal_calculations=true).pressure
        
        # now that we know the whole system, lets go back and see if we can re-solve the known volume (specifically with error) given P and T
        solved_v = PengRobinsonState(chem_comps_with_err, mole_fractions, kij_table_with_err; p=known_p, t=known_t, minimal_calculations=true).volume
        
        error_array[idx] = abs((known_v.err - solved_v.err)/known_v.err) * 100  # units of percentage
    end
    sort!(error_array; rev=true)
    @test -log10(error_array[1]) >= precision  # quantify how precise we want to be (for 10k examples, not very)
end

function test_convenience_methods()
    # test creating matrices
    component_ids = ["CO2", "CH4"]
    get_kij_matrix(PengRobinson(), component_ids)

    missing_component = "ABCD1234<>?F:"
    @test get_kij(PengRobinson(), missing_component, missing_component) == 0
    @test ismissing(get_kij_matrix(PengRobinson(), [missing_component, "CO2"])[1, 2])

    # test easy constructors with database lookups for 4 cases:
    # Manual multicomponent, easy multicomponent, manual single, easy single
    # we can also test the general eos interface here


    # manual multicomponent
    multi_params = ChemicalParameters(component_ids)
    multiparam_kijs = get_kij_matrix(PengRobinson(), component_ids)
    state1 = EOS(PengRobinson(), multi_params, [0.4, 0.6], multiparam_kijs; p=1, t=273.15)
    # easy multicomponent
    state2 = EOS(PengRobinson(), ["CO2", "CH4"], [0.4, 0.6]; p=1, t=273.15)
    @test state1.volume == state2.volume

    # manual single component
    single_param = multi_params[1]
    state3 = EOS(PengRobinson(), single_param; p=1, t=273.15)
    # easy single component
    state4 = EOS(PengRobinson(), "CO2"; t=273.15, p=1)
    @test state3.volume == state4.volume
    
    # also see if forwarding arbitrary kwargs works here
    state5 = EOS(PengRobinson(), "CO2"; t=273.15, p=1, minimal_calculations=true)
    @test state4.volume == state5.volume
end

function test_abstract_getters()
    state = EOS(PengRobinson(), ["CO2", "CH4"], [0.4, 0.6]; p=1.0, t=273.15)
    @test volume(state) == state.volume
    @test pressure(state) == 1
    @test temperature(state) == 273.15
    @test component_names(state) == ["CO2", "CH4"]
    @test mole_fractions(state) == [0.4, 0.6]
    @test components(state) == state.components
end

function test_peng_robinson(; precision=8, speedtests=true)
    @testset "TestPengRobinson.jl" begin
        test_preos_pvt_solvers(; precision=precision, speedtests=speedtests)
        test_preos_pvt_solvers_methods(; precision=precision, speedtests=speedtests)
        
        test_preos_measurements_integration(; precision=precision)
        test_preos_measurements_validity()
        test_convenience_methods()
        test_abstract_getters()
    end
end
