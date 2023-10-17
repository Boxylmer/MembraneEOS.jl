# MembraneEOS.jl

Wrapper for various equation of state (EOS) implementations to be used in the membrane package ecosystem. 

This repository will hopefully soon wrap [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl) for EOS calculations, while providing convenience functions and data for membrane / polymer focused tasks. 


# Sections

## The Gist
To get the properties of a mixture of chemicals in a given composition and conditions (state), we use an equation of state (**EOS**).
This works by...
1. Constructing the parameters you will be needing for the EOS. (see [EOS Parameters](@ref))
2. (optional) Constructing a KIJ matrix to account for deviations from the EOS predictions. (see [KIJ Matrices](@ref))
3. Construction of the EOS model itself. (see [Models](@ref))
4. Querying properties from the model. (see [Exposed Functionality](@ref))

- the two ways to get parameters
    - by name
    - by values
- Constructing mixed models
    - KIJ initialization and setting
- unit limitations


## Quick Start
### Parameter generation
To construct a model, you first construct the parameters that go into the model. using [Peng Robinson](@ref) as an example, which uses "cubic parameters"...
```@example prindex
using MembraneEOS # hide
co2 = CubicParameters("CO2")
nothing # hide
```

They can also be constructed via values directly (see **[CubicParameters](@ref)**)
```@example prindex
ch4 = CubicParameters(190.8, 45.79, 0.012, 16.04) # Tc, Pc, ω, MW
nothing # hide
```
### Model construction
These parameters can then be combined to create an EOS and extract a property from it.
```@example prindex
components = [co2, ch4]
pr = PR(components)
pressure(pr, 22.414, 273.15, [0.5, 0.5]) # returns in mpa
```

### Unideal interactions
KIJ matrices, which encode deviations from the default EOS interactions, follow the same order that your parameters do in the input to the EOS model constructor.
We can use the helper function [`MembraneEOS.initmatrix`](@ref) to create an ideal KIJ matrix for some parameters. 


```@example prindex
kij = MembraneEOS.initmatrix(components) 
kij[1, 2] = 0.09
kij[2, 1] = 0.09
pr = PR(components, kij)
pressure(pr, 22.414, 273.15, [0.5, 0.5])
```

### Fast lookups for common components in gas separations.
Some default data is available for automatically looking up and using EOS parameters and kij values.

See "TestChemicalLookup.jl" for a bit more on how custom databases are used.

If looking up from a data base, you can skip the parameters and simply construct using strings (this will also try to initialize a KIJ matrix)
!!! note
    If you mix directly created parameters with tabulated ones, the KIJ search will not match any directly specified parameters.
    Eventually mixing these two methods of EOS model generation will throw an error to avoid ambiguity and unexpected, hard to debug behavior.

```@example
using MembraneEOS # hide
pr = PR(["CO2", "CH4"])
pressure(pr, 22.414, 273.15, [0.5, 0.5])
```



