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
        Matrix(initial_value*I,matrix_size,matrix_size)
    end
    include("ChemicalLookup.jl")
    include(joinpath("original_eos_implementations", "CubicEOS.jl"))
    # include("PRClapeyronWrapper.jl")

    export ChemicalParameters
    export UnorderedChemicalPair

    export PR
    
    export get_kij
    export get_kij_matrix

    export strip_measurement_to_value

    function __init__()
        initialize_chemical_lookup()
    end
end

using Measurements
model = MembraneEOS.PR([123 ± 1, 124], [56, 78], [0, 0.012])
vol = MembraneEOS.volume(model, 1, 273.15, [0.5, 0.5])