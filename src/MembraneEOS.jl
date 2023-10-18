module MembraneEOS
    using MembraneBase
    using MembraneBase: R_ATM_L_K_MOL
    using Measurements
    using CSV
    import Clapeyron
    using Clapeyron: SingleParam, PairParam
    using StaticArrays
    using LinearAlgebra
    using Roots

    # using Measurements  # causes stack overflow errors somehow??
    """
        initmatrix(components::AbstractVector; [initial_value=0.0])
    Create an ideal KIJ matrix from a list of components or just a number of elements.
    """
    initmatrix(components::AbstractVector; initial_value = 0.0) = initmatrix(length(components); initial_value)
    initmatrix(n; initial_value = 0.0) = ones(n, n) * initial_value

    include("EOSAbstractions.jl")

    include("ChemicalLookup.jl")
    export get_kij
    export get_kij_matrix

    include(joinpath("CubicEOS.jl"))
    export PR
    export CubicParameters
    export CubicModel
    
    # include(joinpath("SLClapeyronWrapper.jl"))
    include(joinpath("SanchezLacombe.jl"))
    export SL
    export SanchezLacombeParameters

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
    export ρTω_activity
    export density_upper_bound
    
    export strip_measurement_to_value

    function __init__()
        initialize_chemical_lookup()
    end
end
