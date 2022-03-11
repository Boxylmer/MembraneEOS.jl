#=
Cubic Equations of State
=#

struct CubicState{MT, CP<:ChemicalParameters, MOLE_T, KIJ_T, PT, VT, TT, BMT, AMT, ZT, BT, AT, ALP_T, PPT, FT, N, L} <: EquationOfState
    model::MT
    components::SVector{N, CP}
    mole_fractions::SVector{N, MOLE_T}
    kij::SMatrix{N, N, KIJ_T, L}
    pressure::PT
    volume::VT 
    temperature::TT

    b_mixed::BMT
    a_mixed::AMT
    z::ZT
    b::SVector{BT}  
    a::SVector{N, AT}
    alpha::SVector{N, ALP_T}
    partial_pressures::SVector{N, PPT}
    fugacities::FT 
end

# Peng Robinson 
struct PengRobinson end
function get_kij(::PengRobinson, component_1::String, component_2::String; default_value=missing)
    if component_1 == component_2 return 0.0 end
    return get_kij(PRKijLookup, component_1, component_2; default_value=default_value)  # maybe add a default value?
end
function get_kij_matrix(::PengRobinson, components::AbstractVector{<:String})
    return get_kij_matrix(PRKijLookup, components)
end 

EOS(::PengRobinson, args...; kwargs...) = CubicState(PengRobinson(), args...; kwargs...)  # forward the PREoS specific constructor to the eos ezpz constructor
"""
Peng-Robinson Equation of State
Wrapper around `CubicState` i.e., `CubicState(model=PengRobinson(), ...)`
"""
PengRobinsonState(args...; kwargs...) = CubicState(PengRobinson(), args...; kwargs...)
get_cubic_eos_constants(::PengRobinson) = 0.45724, 0.07780, 1-√2, 1+√2  # Ωa, Ωb, c1, c2
alpha(::PengRobinson, t, tc, acentric_factor) = (1 + m(PengRobinson(), acentric_factor) * (1 - sqrt(t/tc)))^2
m(::PengRobinson, acentric_factor) = 0.37464 + 1.54226*acentric_factor + 0.26992*acentric_factor^2


"""
    CubicState(model, components, mole_fractions, kijmatrix; p=nothing, v=nothing, t=nothing)


## Calculate a standard Cubic Equation of State

When building a cubic EOS manually, the following information is necessary:
1. A components vector (vector of `ChemicalParameter`s)
    * each `ChemicalParameters` object in the vector will need the following specified unless not required by the specific EOS
        - Critical pressure (**atm**)
        - Critical temperature (**K**)
        - Acentric factor (no units)
        - Molecular weight (**g/mol**)
2. A mole fractions vector (with elements satisfying `<:Number`)
3. A kij matrix (with elements satisfying `<:Number`)
4. Two of the three PVT parameters, with the remaining one calculated by the EOS. 
    - P: Pressure (**atm**)
    - V: Volume (**L/mol**)
    - T: Temperature (**kelvin**)

### Extra options
The following optional flags are available:
* `minimal_calculations`: if set to `true`, will only calculate the state values (pressure, temperature, and volume). 
        This is useful for creating many states without wasting computational effort on downstream calculations that won't be used. 
* `autodiff_mode`: if set to `true`, will remove any measurements from the input parameters, allowing for autodiff types (E.g., ForwardDiff's DualNumbers) to propagate without conflicts
        note that this only works for cases where V and T are known, since the function is only differentiable for P(V, T).
* `component_mode`: if set to `true`, will opt to use the `@uncertain` macro to propagate `measurement` uncertainty. 
        This ensures that `uncertainty_components` will be able to properly account for all sources of error. 
        Otherwise, any iterative result (e.g., when solving volume) will not have any calculation history attached to them, and the components used to calculate them will erroneously show up as zero.

### Convenience constructors
Additionally, 3 convenience constructors are available to simplify state definition for the following scenarios: 
* Parameters are available in the PolymerMembranes (e.g., the *chemical parameter database*, *PREoS KIJ database*, etc).
* Only one component is present.
* Both of the above.
These are documented below. Note that defining two of the three PVT values is still required and setting optional flags is still allowed.
"""
function CubicState(model, components::AbstractVector{<:ChemicalParameters}, mole_fractions, kijmatrix; p=nothing, v=nothing, t=nothing, minimal_calculations=false, autodiff_mode=false, component_mode=false)
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model)
    num_components = length(components)
    # If this state is being called by ForwardDiff, we can't have measurement types present as they don't play nice together
    if autodiff_mode
        p_value = strip_measurement_to_value(p)
        t_value = strip_measurement_to_value(t)
        mf_vals = strip_measurement_to_value(mole_fractions)
        kijm_vals = strip_measurement_to_value(kijmatrix)
        component_vals = strip_measurement_to_value(components)
    else
        p_value = p
        t_value = t
        mf_vals = mole_fractions
        kijm_vals = kijmatrix
        component_vals = components
    end
    
    # get a slight performance boost by making the vector components into static vectors
    _components = SVector{num_components}(component_vals)
    _mole_fractions = SVector{num_components}(mf_vals)
    matrix_size, _ = size(kijmatrix)
    _kijmatrix = SMatrix{matrix_size, matrix_size}(kijm_vals)

    # figure out which parameter wasn't specified and calculate it
    if isnothing(t)
        throw(ErrorException("Solving for temperature is not implemented yet"))
    end
    
    # now that we know we have temperature, lets do all the temperature related calculations we can
    b_values = cubic_b_parameters(omega_b, _components)
    alphas = cubic_alphas(model, t_value, components)
    a_values = cubic_a_parameters(omega_a, alphas, components)

    b_mixed = van_der_waals_mixing_b(b_values, _mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, _mole_fractions, _kijmatrix)

    if isnothing(p) 
        p = cubic_eos_pressure(R_ATM_L_K_MOL, t, v, a_mixed, b_mixed, c1, c2)
    end

    A_i = a_values .* p / (R_ATM_L_K_MOL * t) ^ 2
    B_i = b_values .* p / (R_ATM_L_K_MOL * t)
    A = a_mixed * p / (R_ATM_L_K_MOL * t) ^ 2
    B = b_mixed * p / (R_ATM_L_K_MOL * t)
    z = cubic_eos_compressibility(A, B, c1, c2)

    if isnothing(v) 
        v = z * R_ATM_L_K_MOL * t / p  # z = pv / rt -> v = rtz / p 
    end

    if minimal_calculations  # if we just want pressure, we can stop here.
        partial_pressures = @SVector zeros(num_components)
        fugacities = @SVector zeros(num_components)
    else
        partial_pressures = compute_ideal_partial_pressures(p, _mole_fractions)
        fugacities = cubic_eos_fugacities(_mole_fractions, partial_pressures, z, A, A_i, B, B_i, _kijmatrix, c1, c2)
    end
    CubicState(model, _components, _mole_fractions, _kijmatrix, p, v, t, b_mixed, a_mixed, z, b_values, a_values, alphas, partial_pressures, fugacities)

end

"""
    PengRobinsonState(component::ChemicalParameters)
Define a cubic EOS from a single ChemicalParameters component
"""
function CubicState(model, component::ChemicalParameters; kwargs...)
    return CubicState(model, [component], [1], [0][:, :]; kwargs...)
end

"""
    PengRobinsonState(components::AbstractVector{<:AbstractString}, mole_fractions)
Define a cubic EOS from a vector of chemical names and their corresponding mole fractions.
"""
function CubicState(model, components::AbstractVector{<:AbstractString}, mole_fractions; kwargs...)
    # call the base constructor 
    component_parameters = ChemicalParameters(components)
    kij_matrix = get_kij_matrix(PengRobinson(), components)
    return CubicState(model, component_parameters, mole_fractions, kij_matrix; kwargs...)
end

"""
    PengRobinsonState(component::String)
Define a cubic EOS from a single chemical name.
"""
function CubicState(model, component::String; kwargs...)
    params = ChemicalParameters(component)
    return CubicState(model, params; kwargs...)
end


# define all general computation methods

function cubic_eos_compressibility(A, B, c1, c2)
    d = 1
    c = (c1 + c2 - 1)*B - 1
    b = c1*c2*B^2 - (c1 + c2)*(B^2 + B) + A
    a = -B * (c1*c2*(B^2 + B) + A)

    poly = @SVector [a, b, c, d]
    roots = solve_cubic_eq(poly)
    roots = filterreal(roots)

    if length(roots) == 1
        return roots[1]

    elseif length(roots) == 3
        phis = [root > B ? cubic_eos_fugacity_coefficient(root, A, B, c1, c2) : Inf for root in roots] # if z < B, then that compressibility is impossible
        correct_index = argmin(phis)
        return roots[correct_index]
    else 
        throw(ErrorException("Expected 1 or 3 roots. Got: " * string(length(roots))))
    end
end 

function cubic_eos_pressure(r, t, v, a, b, c1, c2)
    return r*t/(v - b) - a/((v+c1*b)*(v+c2*b))  # note that a and b are the dimensional values, *not* the dimensionless A and B values.
end

function cubic_eos_fugacity_coefficient(mole_fractions, z, A, B, B_i, c1, c2, Aij, idx_i)
    term_1 = (z-1) * B_i[idx_i] / B

    term_2 = -log(z - B)
    
    Aij_summation = sum([Aij[idx_i, idx_j] * mole_fractions[idx_j] for idx_j in eachindex(mole_fractions)])

    if (c1==0 && c2==0) term_3_p1 = -A/Z 
    else term_3_p1 = A / ((c1 - c2) * B) end 
    term_3_p2 = 2 * Aij_summation / A - (B_i[idx_i] / B)
    if (c1==0 && c2==0) term_3_p3 = 1  
    else term_3_p3 = log((z + c2 * B) / (z + c1 * B)) end
    term_3 = term_3_p1 * term_3_p2 * term_3_p3

    fugacity_coeff = exp(term_1 + term_2 + term_3)
    return fugacity_coeff
end

function cubic_eos_fugacity_coefficient(z, A, B, c1, c2)
    return cubic_eos_fugacity_coefficient([1], z, A, B, B, c1, c2, [A], 1)
end

function cubic_eos_fugacities(mole_fractions, partial_pressures, z, A, A_i, B, B_i, kij_matrix, c1, c2) 
    Aij = sqrt.(A_i * A_i') .* (1 .- kij_matrix) # technically inefficient since the matrix is symmetric; unless the compiler is able to recognize this
    fugacities = partial_pressures .* [cubic_eos_fugacity_coefficient(mole_fractions, z, A, B, B_i, c1, c2, Aij, idx_i) for idx_i in eachindex(mole_fractions, partial_pressures)]
    return fugacities
end

van_der_waals_mixing_a(a_values, mole_fractions, kij_matrix) = dot(sqrt.(a_values * a_values') .* (1 .- kij_matrix), mole_fractions * mole_fractions')

van_der_waals_mixing_b(b_values, mole_fractions) = sum([b_values[idx] * mole_fractions[idx] for idx in eachindex(b_values, mole_fractions)])

function cubic_alphas(model, t, components::AbstractVector{<:ChemicalParameters})
    return  SVector{length(components)}([alpha(model, t, critical_temperature(component), acentric_factor(component)) for component in components])
end

function cubic_a_parameters(omega_a, alphas, components::AbstractVector{<:ChemicalParameters})
    return SVector{length(components)}(omega_a * (PolymerMembranes.R_ATM_L_K_MOL^2) * (critical_temperature.(components).^2) ./ critical_pressure.(components) .* alphas)
end

function cubic_b_parameters(omega_b, components::AbstractVector{<:ChemicalParameters})
    SVector{length(components)}(omega_b * R_ATM_L_K_MOL * critical_temperature.(components) ./ critical_pressure.(components))
end


######### stuff to eventuall replace the EOS system above, but for now will exist concurrently until everything is moved over
struct CubicModel{MT, CPT, KIJ_T, N, L}
    modeltype::MT  # eventually add types to PengRobinson, VanDerWaals, etc, and make this specific
    components::SVector{N, CPT}
    kij::SMatrix{N, N, KIJ_T, L}
end

function CubicModel(modeltype, components::Vector{<:ChemicalParameters}, kij::Matrix{<:Number})
    n = length(components)
    _kij = SMatrix{n, n}(kij)
    _components = SVector{n}(components)
    return CubicModel(modeltype, _components, _kij)
end
function CubicModel(modeltype, components::AbstractVector{<:ChemicalParameters})
    return CubicModel(modeltype, components, initmatrix(components))
end
function CubicModel(modeltype, components::AbstractVector{<:String})
    kij = get_kij_matrix(modeltype, components)
    component_parameters = ChemicalParameters(components)
    return CubicModel(modeltype, component_parameters, kij)
end
function CubicModel(modeltype, component::ChemicalParameters)
    return CubicModel(modeltype, [component])
end
function CubicModel(modeltype, component::String)
    return CubicModel(modeltype, [component])
end

# Need to define an easy constructor(s) for this and a bunch of functions like
function compressibility_factor(model::CubicModel, p, t, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    A = a_mixed * p / (R_ATM_L_K_MOL * t) ^ 2
    B = b_mixed * p / (R_ATM_L_K_MOL * t)
    return cubic_eos_compressibility(A, B, c1, c2)
end

function VT_compressibility_factor(model::CubicModel, v, t, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    p = cubic_eos_pressure(R_ATM_L_K_MOL, t, v, a_mixed, b_mixed, c1, c2)

    A = a_mixed * p / (R_ATM_L_K_MOL * t) ^ 2
    B = b_mixed * p / (R_ATM_L_K_MOL * t)
    return cubic_eos_compressibility(A, B, c1, c2)
end

function pressure(model::CubicModel, v, t, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    # now that we know we have temperature, lets do all the temperature related calculations we can
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    return cubic_eos_pressure(R_ATM_L_K_MOL, t, v, a_mixed, b_mixed, c1, c2)
end

function volume(model::CubicModel, p, t, mole_fractions=[1])
    z = compressibility_factor(model, p, t, mole_fractions)
    v = z * R_ATM_L_K_MOL * t / p
    return v
end

function fugacity(model::CubicModel, p, t, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)
    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    A_i = a_values .* p / (R_ATM_L_K_MOL * t) ^ 2
    B_i = b_values .* p / (R_ATM_L_K_MOL * t)
    A = a_mixed * p / (R_ATM_L_K_MOL * t) ^ 2
    B = b_mixed * p / (R_ATM_L_K_MOL * t)
    z = cubic_eos_compressibility(A, B, c1, c2)
    partial_pressures = compute_ideal_partial_pressures(p, mole_fractions)
    return cubic_eos_fugacities(mole_fractions, partial_pressures, z, A, A_i, B, B_i, model.kij, c1, c2) 
end
