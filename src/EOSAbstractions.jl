
abstract type MEOSModel end

# shared methods
"""
    pressure(model, v, t, [mole_fractions])

Get the pressure of the state in MPa given a volume (L/mol), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
MembraneEOS.pressure(::MEOSModel, _, _, args...) = throw(ErrorException("Pressure not implemented for this model."))

"""
    volume(model, p, t, [mole_fractions])

Get the volume of the state in MPa given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
volume(::MEOSModel, _, _, args...) = throw(ErrorException("Volume not implemented for this model."))

"""
    fugacity(model, p, t, [mole_fractions])

Get a vector of fugacities for each component in the state in MPa given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
fugacity(::MEOSModel, _, _, args...) = throw(ErrorException("Fugacity not implemented for this model."))

"Get a vector of molecular weights for the system in g/cm3."
molecular_weight(::MEOSModel, args...) = throw(ErrorException("Molecular weight not implemented for this model."))

"""
    compressibility_factor(model, p, t, [mole_fractions])

Get the compressibility factor (z) of the state given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
compressibility_factor(::MEOSModel, _, _, args...) = throw(ErrorException("Compressibility factor not implemented for this model."))

"""
    compressibility_factor(model, v, t, [mole_fractions])

Get the compressibility factor (z) of the state given a volume (L/mol), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
VT_compressibility_factor(::MEOSModel, _, _, args...)  = throw(ErrorException("VT compressibility factor not implemented for this model."))



# Mainly lattice fluid methods
"""
    mass_density(model, p, t, [mole_fractions])

Get the mass density of the state in g/cm3 given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
mass_density(::MEOSModel, args...) = throw(ErrorException("Mass density not implemented for this model."))

"""
    VT_mass_density(model, v, t, [mole_fractions])

Get the density of the state in g/cm3 given a volume (L/mol), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
VT_mass_density(::MEOSModel, args...) = throw(ErrorException("VT Mass density not implemented for this model."))

"""
    chemical_potential(model, p, t, [mole_fractions])

Get the chemical potential of the state in J/mol given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
chemical_potential(::MEOSModel, args...) = throw(ErrorException("Chemical potential not implemented for this model."))

"""
    chemical_potential_res(model, p, t, [mole_fractions])

Get the residual chemical potential of the state in J/mol given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
chemical_potential_res(::MEOSModel, args...) = throw(ErrorException("Residual chemical potential not implemented for this model."))

"""
    VT_chemical_potential(model, v, t, [mole_fractions])

Get the chemical potential of the state in J/mol given a volume (L/mol), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
VT_chemical_potential(::MEOSModel, args...) = throw(ErrorException("VT chemical potential not implemented for this model."))

"""
    ρTω_chemical_potential(model, ρ, t, [mass_fractions])

Get the chemical potential of the state in J/mol given a density (g/cm3), temperature (K) and composition in mass fraction (can be omitted if the system is pure).
"""
ρTω_chemical_potential(::MEOSModel, args...) = throw(ErrorException("ρTω chemical potential not implemented for this model."))

"""
    ρTω_chemical_potential_res(model, ρ, t, [mass_fractions])

Get the residual chemical potential of the state in J/mol given a density (g/cm3), temperature (K) and composition in mass fraction (can be omitted if the system is pure).
"""
ρTω_chemical_potential_res(::MEOSModel, args...) = throw(ErrorException("Residual ρTω chemical potential not implemented for this model."))

"""
    activity(model, p, t, [mole_fractions])

Get the chemical activity of the state given a pressure (MPa), temperature (K) and composition in mole fraction (can be omitted if the system is pure).
"""
activity(::MEOSModel, args...) = throw(ErrorException("Activity not implemented for this model."))

"""
    ρTω_activity(model, ρ, t, [mass_fractions])

Get the chemical activity of the state given a density (g/cm3), temperature (K) and composition in mass fraction (can be omitted if the system is pure).
"""
ρTω_activity(::MEOSModel, args...) = throw(ErrorException("ρTω activity not implemented for this model."))

"""
    density_upper_bound(model, [mass_fractions])

Get the maximum possible mass density that the system can exhibit (e.g., for Sanchez Lacombe, this is its mixed characteristic density, for cubic EOSs, it is infinite), given its composition (omitted if system is pure).
"""
density_upper_bound(::MEOSModel, args...) = throw(ErrorException("Density upper bound not implemented for this model."))

# helpers
"""
    get_kij(modeltype, component_1::String, component_2::String; [default_value=0.0])

Attempt to look up the kij value between two components for a given EOS model type. 
- Returns `missing` if the value is not in the database.
- Equal pairings (i.e., "CO2" and "CO2" will return the `default_value`).
"""
get_kij(::MEOSModel, component_1::String, component_2::String; _=0.0) = throw(ErrorException("KIJ lookup not implemented for this model."))

"""
    get_kij_matrix(modeltype, components::AbstractVector{<:String})
Attempt to look up all pairings in a list of components and arrange them into a corresponding KIJ matrix. 
- When constructing an entire matrix at once, it is assumed you want to ignore missing values, so component pairings not in the database are assumed to interact ideally instead of being listed as `missing`.
e.g., 
```
get_kij_matrix(PR(), ["CO2", "CH4"])

> 2×2 Matrix{Union{Missing, Float64}}:
> 0.0   0.09
> 0.09  0.0
```

```
get_kij_matrix(PR(), ["CO2", "Something not in the database"])

> 2×2 Matrix{Union{Missing, Float64}}:
> 0.0  0.0
> 0.0  0.0
```
"""
get_kij_matrix(::MEOSModel, _=AbstractVector{<:String}) = throw(ErrorException("KIJ matrix initialization lookup not implemented for this model."))