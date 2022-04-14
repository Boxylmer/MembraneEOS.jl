"""
WHEN UPDATING CHEMICAL LOOKUP TABLES WITH NEW VALUE TYPES

The following must be updated

- Add the indices in ChemicalLookupHelper
- Add the values in the ChemicalParameters struct itself
- Update the ChemicalParameters constructor functions
- Update the docs table with the units of what was added
- Add the kwarg to the ChemicalParameters manual constructor
- Create the getter function for the parameters e.g., `characteristic_density(::ChemicalParameters)` 
- Update `strip_measurement_to_value`
- Update the tests :)
"""


"""Structs and globals"""
struct UnorderedChemicalPair; a::String; b::String; end  # courtesy of Gandalf from HoJ
Base.hash(x::UnorderedChemicalPair, h::UInt) = hash(x.a < x.b ? (x.a, x.b) : (x.b, x.a), h)
Base.isequal(x::UnorderedChemicalPair, y::UnorderedChemicalPair) = x.a == y.a ? x.b == y.b : (x.a == y.b && x.b == y.a)

const ChemicalNameToRowLookup = Dict{String, Int64}()
# const PolymerNameToRowLookup = Dict{String, Int64}()
const PRKijLookup = Dict{UnorderedChemicalPair, Number}()
const SLKijLookup = Dict{UnorderedChemicalPair, Number}()

"""
Helper functions specific only to looking up these CSV files 
    (e.g., @DIR only applies to the file it's written in, so we can't really move it to HelperFunctions.jl)
"""

# get a CSV file within the folder this function is written in
function read_chemical_database_csv_file(filename::String)
    read_file = read(
        joinpath(@__DIR__, filename)
    )
    return CSV.File(read_file; header=true)
end

# create a dict that will take all values present in the index_column and return the row they exist on
# assumes all data in the indexing column will be unique and a string 
function create_csv_indexing_dict!(lookup_dict::AbstractDict, csv_name::AbstractString, index_column::Integer)
    read_file = read(joinpath(@__DIR__, csv_name))
    csv_rows = CSV.Rows(read_file, header=true)
    for (row_id, row) in enumerate(csv_rows)
        name = row[index_column]
        lookup_dict[name] = row_id + 1
    end
    return nothing
end

function parse_measurement_from_csv_row(row, value_key, error_key)
    val = row[value_key]
    err = row[error_key]
    if !ismissing(err)
        val = val ± err
    end
    return val
end


"""
Chemical database functions and lookups
"""
function find_chemical_row_index(chemical_name::String)
    chemical_row = get!(ChemicalNameToRowLookup, chemical_name, -1)  # -1 meaning there is no row with this data
    if chemical_row == -1  # get! adds the name to the dict, so we need to remove it if it returned -1
        delete!(ChemicalNameToRowLookup, chemical_name)
        return missing
    end
    return chemical_row
end

function is_chemical_in_database(chemical_name::String)
    if ismissing(find_chemical_row_index(chemical_name)) 
        return false
    else 
        return true
    end
end

module ChemicalLookupHelper
    const POLYMER_PARAM_DB = joinpath("..", "ChemicalData", "polymer_parameters.csv")
    const CHEMICAL_PARAM_DB = joinpath("..", "ChemicalData", "chemical_parameters.csv")
    const PREOS_KIJ_DB = joinpath("..", "ChemicalData", "PR_binary_interactions.csv")
    const SLEOS_KIJ_DB = joinpath("..", "ChemicalData", "SL_binary_interactions.csv")
    
    # Chemical DB lookup cols
    const CHEM_NAME_COL = 1
    const CHEM_MW_COL = 2
    const CHEM_MW_ERR_COL = 3
    const CHEM_TSTAR_COL = 4
    const CHEM_TSTAR_ERR_COL = 5
    const CHEM_PSTAR_COL = 6
    const CHEM_PSTAR_ERR_COL = 7
    const CHEM_RHOSTAR_COL = 8
    const CHEM_RHOSTAR_ERR_COL = 9
    const CHEM_TC_COL = 10
    const CHEM_TC_ERR_COL = 11
    const CHEM_PC_COL = 12
    const CHEM_PC_ERR_COL = 13
    const CHEM_W_COL = 14
    const CHEM_W_ERR_COL = 15

    # PREoS KIJ DB lookup cols
    const PRKIJ_NAME_1 = 1
    const PRKIJ_NAME_2 = 2
    const PRKIJ_KIJVAL = 3
    const PRKIJ_UNCERT = 4

    # SLEoS KIJ DB lookup cols
    const SLKIJ_NAME_1 = 1
    const SLKIJ_NAME_2 = 2
    const SLKIJ_KIJVAL = 3
    const SLKIJ_UNCERT = 4
end

struct ChemicalParameters{MWT, CTT, CPT, AFT, CCTT, CCPT, CDT}
    name::Union{String, Missing}
    molecular_weight::MWT  # g/mol

    critical_temperature::CTT  # K
    critical_pressure::CPT  # atm
    acentric_factor::AFT  # no units

    characteristic_temperature::CCTT  # K
    characteristic_pressure::CCPT  # MPa
    characteristic_density::CDT  # g/cm3
end

function ChemicalParameters(chemical_name::String)
    chemical_row = find_chemical_row_index(chemical_name)
    if ismissing(chemical_row) return missing end 
    read_file = read(joinpath(@__DIR__, ChemicalLookupHelper.CHEMICAL_PARAM_DB))
    csv_row = CSV.File(read_file, header=true, skipto=chemical_row, limit=1)[1]
    name = chemical_name
    mw = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_MW_COL, ChemicalLookupHelper.CHEM_MW_ERR_COL)
    tc = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_TC_COL, ChemicalLookupHelper.CHEM_TC_ERR_COL)
    pc = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_PC_COL, ChemicalLookupHelper.CHEM_PC_ERR_COL)
    w  = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_W_COL , ChemicalLookupHelper.CHEM_W_ERR_COL )
    tstar = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_TSTAR_COL, ChemicalLookupHelper.CHEM_TSTAR_ERR_COL)
    pstar = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_PSTAR_COL, ChemicalLookupHelper.CHEM_PSTAR_ERR_COL)
    rhostar = parse_measurement_from_csv_row(csv_row, ChemicalLookupHelper.CHEM_RHOSTAR_COL, ChemicalLookupHelper.CHEM_RHOSTAR_ERR_COL)
    parsed_chemical = ChemicalParameters(name, mw, tc, pc, w, tstar, pstar, rhostar)
    return parsed_chemical
end

ChemicalParameters(chemical_names::AbstractVector{<:String}) = ChemicalParameters.(chemical_names)

# this constructor is necessary in the case that you want to partially replace a few parameters from the database default. 
# Thus, a name is required, otherwise we just want the default constructor.
"""
Define a set of parameters that describe a chemical. This method uses an optional (positional) name, and additional keyword arguments (with units) described below.
!!! note 
    It is highly reccomended to stay within the units commonly used through database lookups to avoid the risk of unit errors.

| Parameter                        |     unit    |
|----------------------------------|-------------|   
| `molecular_weight`               | g/mol       |
| `critical_temperature`           | K           |
| `critical_pressure`              | atm         |
| `acentric factor`                | no units    |
| `characteristic_temperature`     | K           |
| `characteristic_pressure`        | MPa         |
| `characteristic_density`         | g/cm3       |

"""
function ChemicalParameters(chemical_name::Union{String, Missing} = missing; fill_missing_with_database=false,
    molecular_weight=missing, critical_temperature=missing, critical_pressure=missing, acentric_factor=missing,
    characteristic_temperature=missing, characteristic_pressure=missing, characteristic_density=missing,
    )

    if fill_missing_with_database && !ismissing(chemical_name)
        default_parameters = ChemicalParameters(chemical_name)
        return ChemicalParameters(
            chemical_name,
            first_nonmissing_parameter(molecular_weight, default_parameters.molecular_weight),
            first_nonmissing_parameter(critical_temperature, default_parameters.molecular_weight),
            first_nonmissing_parameter(critical_pressure, default_parameters.critical_pressure),
            first_nonmissing_parameter(acentric_factor, default_parameters.acentric_factor),
            first_nonmissing_parameter(characteristic_temperature, default_parameters.characteristic_temperature),
            first_nonmissing_parameter(characteristic_pressure, default_parameters.characteristic_pressure),
            first_nonmissing_parameter(characteristic_density, default_parameters.characteristic_density)
        )
    else
        return ChemicalParameters(chemical_name, molecular_weight, critical_temperature, critical_pressure, acentric_factor,
        characteristic_temperature, characteristic_pressure, characteristic_density)
    end
end 

name(params::ChemicalParameters) = params.name
molecular_weight(params::ChemicalParameters) = params.molecular_weight
critical_temperature(params::ChemicalParameters) = params.critical_temperature
critical_pressure(params::ChemicalParameters) = params.critical_pressure
acentric_factor(params::ChemicalParameters) = params.acentric_factor
characteristic_temperature(params::ChemicalParameters) = params.characteristic_temperature
characteristic_pressure(params::ChemicalParameters) = params.characteristic_pressure
characteristic_density(params::ChemicalParameters) = params.characteristic_density

function MembraneBase.strip_measurement_to_value(obj::ChemicalParameters)
    return ChemicalParameters(obj.name, 
        strip_measurement_to_value(obj.molecular_weight), 
        strip_measurement_to_value(obj.critical_temperature), 
        strip_measurement_to_value(obj.critical_pressure), 
        strip_measurement_to_value(obj.acentric_factor), 
        strip_measurement_to_value(obj.characteristic_temperature), 
        strip_measurement_to_value(obj.characteristic_pressure), 
        strip_measurement_to_value(obj.characteristic_density)
    )
end
function MembraneBase.strip_measurement_to_value(obj::AbstractVector{<:ChemicalParameters})
    return strip_measurement_to_value.(obj)
end

# function contains_measurement_type(obj::ChemicalParameters)
#     if typeof(obj.molecular_weight) <: Measurement return true
#     elseif typeof(obj.critical_temperature) <: Measurement return true
#     elseif typeof(obj.critical_pressure) <: Measurement return true
#     elseif typeof(obj.acentric_factor) <: Measurement return true
#     elseif typeof(obj.characteristic_temperature) <: Measurement return true
#     elseif typeof(obj.characteristic_pressure) <: Measurement return true
#     elseif typeof(obj.characteristic_density) <: Measurement return true
#     end
#     return false
# end


"""
KIJ database reading functions and lookups
"""
function create_kij_dict_from_csv!(lookup_dict::AbstractDict, csv_name::AbstractString, 
    index_column_1::Integer, index_column_2::Integer, kij_val_column::Integer, kij_err_column::Integer)
    read_file = read(joinpath(@__DIR__, csv_name))
    csv_rows = CSV.File(read_file, header=true) #  CSV.Rows(read_file, header=true)
    for row in csv_rows
        name_1 = row[index_column_1]
        name_2 = row[index_column_2]
        chemical_name_pair = UnorderedChemicalPair(name_1, name_2)
        kij_measurement = parse_measurement_from_csv_row(row, kij_val_column, kij_err_column)
        lookup_dict[chemical_name_pair] = kij_measurement
    end
    return nothing
end

function kij_name_to_value_lookup(database::AbstractDict, chemical_pair::UnorderedChemicalPair)
    found_row = get!(database, chemical_pair, -1)  # -1 meaning there is no row with this data
    if found_row == -1  # get! adds the name to the dict, so we need to remove it if it returned -1
        delete!(database, chemical_pair)
        return missing
    end
    return found_row
end

function get_kij(database_dict::Dict{<:UnorderedChemicalPair, <:Number}, chemical_1::String, chemical_2::String; missing_value=missing)
    pairing = UnorderedChemicalPair(chemical_1, chemical_2)
    kij = kij_name_to_value_lookup(database_dict, pairing)
    if ismissing(kij) return missing_value end     
    return kij    
end

function get_kij_matrix(database_dict::Dict{<:UnorderedChemicalPair, <:Number}, chemical_strings::AbstractVector{<:String}; missing_value=0, ideal_value=0)
    matrix_size = length(chemical_strings)
    # kijmat = Matrix{Any}(ideal_value, matrix_size, matrix_size)
    kijmat = zeros(Union{Float64, Missing}, (matrix_size, matrix_size))
    for (idx1, comp_1) in enumerate(chemical_strings)
        for (idx2, comp_2) in enumerate(chemical_strings)
            if idx1 == idx2
                kijmat[idx1, idx2] = ideal_value
                continue
            end
            kijval = get_kij(database_dict, comp_1, comp_2; missing_value=missing_value)
            idx1, idx2
            kijmat[idx1, idx2] = kijval
            kijmat[idx2, idx1] = kijval
        end
    end
    return kijmat
end


# main call to initialize all of the data
# todo maybe move this to the top of the file and see if tests still run
function initialize_chemical_lookup()
    create_csv_indexing_dict!(ChemicalNameToRowLookup, ChemicalLookupHelper.CHEMICAL_PARAM_DB, ChemicalLookupHelper.CHEM_NAME_COL)
    create_kij_dict_from_csv!(
        PRKijLookup, ChemicalLookupHelper.PREOS_KIJ_DB, 
        ChemicalLookupHelper.PRKIJ_NAME_1, ChemicalLookupHelper.PRKIJ_NAME_2, 
        ChemicalLookupHelper.PRKIJ_KIJVAL, ChemicalLookupHelper.PRKIJ_UNCERT)
    create_kij_dict_from_csv!(
        SLKijLookup, ChemicalLookupHelper.SLEOS_KIJ_DB,
        ChemicalLookupHelper.SLKIJ_NAME_1, ChemicalLookupHelper.SLKIJ_NAME_2,
        ChemicalLookupHelper.SLKIJ_KIJVAL, ChemicalLookupHelper.SLKIJ_UNCERT
    )
end
