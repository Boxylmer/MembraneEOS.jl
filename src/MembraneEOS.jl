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
    export CubicParameters
    export CubicModel
    
    include(joinpath("SLClapeyronWrapper.jl"))
    export SL

    # shared methods
    export pressure
    export volume
    export fugacity
    export molecular_weight
    export compressibility_factor
    export VT_compressibility_factor

    # lattice fluid methods
    export mass_density
    export VT_mass_density
    export chemical_potential
    export chemical_potential_res
    export VT_chemical_potential
    export ρTω_chemical_potential
    export ρTω_chemical_potential_res
    export activity
    export activity_res
    export ρTω_activity
    export ρTω_activity_res

    export get_kij
    export get_kij_matrix

    export strip_measurement_to_value

    function __init__()
        initialize_chemical_lookup()
    end
end
