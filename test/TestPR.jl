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
    function randomize(var::Number)
        # randn() gives a number with average 0 and standard deviation 1
        multiplier = abs(randn(1)[1])
        return abs(var * multiplier / 2) + var  # most things need to stay positive
    end

    num_iters = 10000 # see above conversations
    Random.seed!(4321)  # repeatability
    error_array = zeros(num_iters)
    
    
    for idx in 1:num_iters
        # define some working inputs for a realistic system
        nitrogen_with_err = MembraneEOS.CubicParameters(randomize(126.2) ± randomize(1), randomize(33.5) ± randomize(1), randomize(0.04) ± randomize(0.005), missing)
        methane_with_err = MembraneEOS.CubicParameters(randomize(190.8) ± randomize(1), randomize(45.79) ± randomize(3), randomize(0.012) ± randomize(0.003), missing)
        chem_comps_with_err = [nitrogen_with_err, methane_with_err]
        kij = randomize(0.03) ± randomize(0.01)
        kij_table_with_err = [(0. ± 0.) (kij); kij (0. ± 0.)]
        mfrac_1 = randomize(0.04) ± randomize(0.1)
        mfrac_2 = (1 - mfrac_1.val) ± randomize(0.1)
        mole_fractions = [mfrac_1, mfrac_2]
        known_v = randomize(30.08405) ± randomize(0.2)
        known_t = randomize(273.15) ± randomize(1)
        
        model = PR(chem_comps_with_err, kij_table_with_err)


        # given this system, calculate the pressure given V and T (the "trivial" case)
        known_p = pressure(model, known_v, known_t, mole_fractions)
        
        # now that we know the whole system, lets go back and see if we can re-solve the known volume (specifically with error) given P and T
        solved_v = volume(model, known_p, known_t, mole_fractions)
        
        error_array[idx] = abs((known_v.err - solved_v.err)/known_v.err) * 100  # units of percentage

        if error_array[idx] > 0.1
            "A single random example failed to return the correct uncertainty or value. (Did the PREoS get the right phase?)"
            @show known_p
            @show solved_v
            @show known_v
            @show error_array[idx]
            @show kij_table_with_err
            @show mole_fractions
        end

    end
    sort!(error_array; rev=true)
    @test -log10(error_array[1]) >= precision  # quantify how precise we want to be (for 10k examples, not very)
end

@testset "TestPR.jl" begin
    model_co2 = PR("CO2")
    model_co2_array = PR(["CO2"])
    @test pressure(model_co2, 4, 273.15) == pressure(model_co2_array, 4, 273.15)


    nitrogen = MembraneEOS.CubicParameters(126.2, 33.5, 0.04, missing)
    methane = MembraneEOS.CubicParameters(190.8, 45.79, 0.012, missing)

    
    chemical_components = [nitrogen, methane]
    mole_fractions = [0.4, 0.6]
    model = PR(chemical_components, MembraneEOS.initmatrix(chemical_components))
    @test ismissing(molecular_weight(model)[1])
    v_0 = 4
    p = pressure(model, v_0, 273.15, mole_fractions)
    v = volume(model, p, 273.15, mole_fractions)
    @test v ≈ 4 

    # test some constructors and make sure they get the results shown before
    nitrogen_methane_pc = MembraneEOS.critical_pressure.([nitrogen, methane])
    nitrogen_methane_tc = MembraneEOS.critical_temperature.([nitrogen, methane])
    nitrogen_methane_ω = MembraneEOS.acentric_factor.([nitrogen, methane])
    nitrogen_methane_mw = [28, 16.04]
    model_generic_no_mw = PR(nitrogen_methane_pc, nitrogen_methane_tc, nitrogen_methane_ω)
    model_generic_with_mw = PR(nitrogen_methane_pc, nitrogen_methane_tc, nitrogen_methane_ω, nitrogen_methane_mw)
    @test pressure(model_generic_no_mw, v_0, 273.15, mole_fractions) == 
          pressure(model_generic_with_mw, v_0, 273.15, mole_fractions) ==
          p 
    
    test_preos_measurements_validity()
end