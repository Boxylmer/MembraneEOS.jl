module MembraneEOS
    using MembraneBase
    using MembraneBase: R_ATM_L_K_MOL
    using CSV
    using Clapeyron
    using StaticArrays
    using LinearAlgebra

    # using Measurements  # causes stack overflow errors somehow??
    function initmatrix(components::AbstractVector; initial_value = 0.0)
        matrix_size = length(components)
        ones(matrix_size,matrix_size) * initial_value
    end
    
    include("ChemicalLookup.jl")
    
    include(joinpath("CubicEOS.jl"))
    export PR
    include(joinpath("SLClapeyronWrapper.jl"))
    export SL

    export mass_density
    export chemical_potential
    export pressure
    export volume
    
    export get_kij
    export get_kij_matrix

    export strip_measurement_to_value

    function __init__()
        initialize_chemical_lookup()
    end
end

# using Measurements
# model = MembraneEOS.PR([123 ± 1, 124], [563, 78], [0, 0.012±0.2])
# vol = MembraneEOS.volume(model, 1, 273.15, [0.5, 0.5])

# model = MembraneEOS.PR("CO2")
# vol = MembraneEOS.volume(model, 1, 273.15) 
# pres = MembraneEOS.pressure(model, vol, 273.15)

# model = MembraneEOS.PR(["CO2", "CH4"])
# vol = MembraneEOS.volume(model, 15±0.1, 273.155±0.1, [0.5±0.1, 0.5±0.1])


# polycarbonate T* = 755, P* = 534, r0 = 1.275, mw = 1.00E+09 
# using Measurements
