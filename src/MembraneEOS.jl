module MembraneEOS
    using MembraneBase
    using CSV
    using Clapeyron
    # using Measurements  # causes stack overflow errors somehow??

    include("ChemicalLookup.jl")
    include("ClapeyronWrapper.jl")

    export ChemicalParameters
    export UnorderedChemicalPair
    
    export get_kij
    export get_kij_matrix

    export strip_measurement_to_value

    function __init__()
        initialize_chemical_lookup()
    end
end
