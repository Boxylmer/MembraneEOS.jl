#=
Cubic Equations of State
=#
struct CubicParameters{TCT, PCT, ACT, MWT}
    critical_temperature_k::TCT
    critical_pressure_atm::PCT
    acentric_factor::ACT
    molecular_weight::MWT
end
CubicParameters(tc_k, pc_atm, ω) = CubicParameters(tc_k, pc_atm, ω, missing)
CubicParameters(cp::ChemicalParameters) = CubicParameters(critical_temperature(cp), critical_pressure(cp), acentric_factor(cp), molecular_weight(cp))
critical_temperature(cp::CubicParameters) = cp.critical_temperature_k
critical_pressure(cp::CubicParameters) = cp.critical_pressure_atm
acentric_factor(cp::CubicParameters) = cp.acentric_factor
molecular_weight(cp::CubicParameters) = cp.molecular_weight

struct CubicModel{MT, CPT, KIJ_T, N, L}
    modeltype::MT 
    components::SVector{N, CPT}
    kij::SMatrix{N, N, KIJ_T, L}
end
function CubicModel(modeltype, components::Vector{<:CubicParameters}, kij::Matrix{<:Union{<:Number, Missing}})
    n = length(components)
    _kij = SMatrix{n, n}(kij)
    _components = SVector{n}(components)
    return CubicModel(modeltype, _components, _kij)
end
function CubicModel(modeltype, components::AbstractVector{<:CubicParameters})
    return CubicModel(modeltype, components, initmatrix(components))
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
get_cubic_eos_constants(::PengRobinson) = 0.45724, 0.07780, 1-√2, 1+√2  # Ωa, Ωb, c1, c2
alpha(::PengRobinson, t, tc, acentric_factor) = (1 + m(PengRobinson(), acentric_factor) * (1 - sqrt(t/tc)))^2
m(::PengRobinson, acentric_factor) = 0.37464 + 1.54226*acentric_factor + 0.26992*acentric_factor^2

function PR(chemical::String)
    component_parameters = ChemicalParameters(chemical)
    cubic_parameters = CubicParameters(component_parameters)
    return PR(cubic_parameters)
end
function PR(chemicals::AbstractVector{String})
    KIJ_matrix = get_kij_matrix(PengRobinson(), chemicals)
    component_parameters = ChemicalParameters(chemicals)
    cubic_parameters = CubicParameters.(component_parameters)
    return PR(cubic_parameters, KIJ_matrix)
end
function PR(pc_atm::AbstractVector, tc_k::AbstractVector, omega::AbstractVector, KIJ_matrix=nothing)
	params = [CubicParameters(tc_k[i], pc_atm[i], omega[i]) for i in eachindex(tc_k, pc_atm, omega)]
    return PR(params, KIJ_matrix)
end
PR(pc_atm::Number, tc_k::Number, ω::Number) = PR([pc_atm], [tc_k], [ω])
PR(params::CubicParameters) = PR([params])
function PR(params::AbstractVector{<:CubicParameters}, KIJ_matrix=nothing)  # base method
    if isnothing(KIJ_matrix)
        KIJ_matrix = initmatrix(params)
    end
    return CubicModel(PengRobinson(), params, KIJ_matrix)
end

# define all general computation methods
function compute_ideal_partial_pressures(pressure, mole_fractions)
    N = length(mole_fractions)
    partial_pressures = SVector{N}(pressure .* mole_fractions)
    return partial_pressures
end

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

function cubic_alphas(model, t, components::AbstractVector{<:CubicParameters})
    return  SVector{length(components)}([alpha(model, t, critical_temperature(component), acentric_factor(component)) for component in components])
end

function cubic_a_parameters(omega_a, alphas, components::AbstractVector{<:CubicParameters})
    return SVector{length(components)}(omega_a * (R_ATM_L_K_MOL^2) * (critical_temperature.(components).^2) ./ critical_pressure.(components) .* alphas)
end

function cubic_b_parameters(omega_b, components::AbstractVector{<:CubicParameters})
    SVector{length(components)}(omega_b * R_ATM_L_K_MOL * critical_temperature.(components) ./ critical_pressure.(components))
end



# functionality

"Get pressure in atm"
function pressure(model::CubicModel, v_l_mol, t_k, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    # now that we know we have temperature, lets do all the temperature related calculations we can
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t_k, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    return cubic_eos_pressure(R_ATM_L_K_MOL, t_k, v_l_mol, a_mixed, b_mixed, c1, c2)
end

"Get volume in L/mol"
function volume(model::CubicModel, p_atm, t_k, mole_fractions=[1])
    z = compressibility_factor(model, p_atm, t_k, mole_fractions)
    v = z * R_ATM_L_K_MOL * t_k / p_atm
    return v
end

"Get a vector of fugacities in atm"
function fugacity(model::CubicModel, p_atm, t_k, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t_k, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)
    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    A_i = a_values .* p_atm / (R_ATM_L_K_MOL * t_k) ^ 2
    B_i = b_values .* p_atm / (R_ATM_L_K_MOL * t_k)
    A = a_mixed * p_atm / (R_ATM_L_K_MOL * t_k) ^ 2
    B = b_mixed * p_atm / (R_ATM_L_K_MOL * t_k)
    z = cubic_eos_compressibility(A, B, c1, c2)
    partial_pressures = compute_ideal_partial_pressures(p_atm, mole_fractions)
    return cubic_eos_fugacities(mole_fractions, partial_pressures, z, A, A_i, B, B_i, model.kij, c1, c2) 
end

function compressibility_factor(model::CubicModel, p_atm, t_k, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t_k, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    A = a_mixed * p_atm / (R_ATM_L_K_MOL * t_k) ^ 2
    B = b_mixed * p_atm / (R_ATM_L_K_MOL * t_k)
    return cubic_eos_compressibility(A, B, c1, c2)
end

function VT_compressibility_factor(model::CubicModel, v_l_mol, t_k, mole_fractions=[1])
    omega_a, omega_b, c1, c2 = get_cubic_eos_constants(model.modeltype)
    b_values = cubic_b_parameters(omega_b, model.components)
    alphas = cubic_alphas(model.modeltype, t_k, model.components)
    a_values = cubic_a_parameters(omega_a, alphas, model.components)

    b_mixed = van_der_waals_mixing_b(b_values, mole_fractions)
    a_mixed = van_der_waals_mixing_a(a_values, mole_fractions, model.kij)

    p = cubic_eos_pressure(R_ATM_L_K_MOL, t_k, v_l_mol, a_mixed, b_mixed, c1, c2)

    A = a_mixed * p / (R_ATM_L_K_MOL * t_k) ^ 2
    B = b_mixed * p / (R_ATM_L_K_MOL * t_k)
    return cubic_eos_compressibility(A, B, c1, c2)
end
