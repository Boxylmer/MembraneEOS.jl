struct SanchezLacombe end

struct SanchezLacombeParameters{CPT, CTT, CDT, MWT}
    characteristic_pressure_mpa::CPT
    characteristic_temperature_k::CTT
    characteristic_density_g_cm3::CDT
    molecular_weight::MWT
end
SanchezLacombeParameters(cp::ChemicalParameters) = SanchezLacombeParameters(characteristic_pressure(cp), characteristic_temperature(cp), characteristic_density(cp), molecular_weight(cp))
characteristic_pressure(slp::SanchezLacombeParameters) = slp.characteristic_pressure_mpa
characteristic_temperature(slp::SanchezLacombeParameters) = slp.characteristic_temperature_k
characteristic_density(slp::SanchezLacombeParameters) = slp.characteristic_density_g_cm3
molecular_weight(slp::SanchezLacombeParameters) = slp.molecular_weight


# core math functions

function sanchez_lacombe_reduced_pressure(reduced_temperature, reduced_density, φ, r_i)
    return -reduced_temperature * (log(1-reduced_density) + (1 - sum(φ ./ r_i)) * reduced_density) - reduced_density^2
end

function sanchez_lacombe_reduced_density(components, kij, mole_fractions, temperature, reduced_temperature, reduced_pressure, φ, r_i0, r_i, pure_characteristic_volume_i; tol=1e-12, maxiters=1000)

    target_function(reduced_density) = 1 - exp(-(reduced_density^2 + reduced_pressure)/reduced_temperature - reduced_density*(1 - sum(φ ./ r_i))) - reduced_density
    reduced_densities = find_zeros(target_function, 0, 1)
    if length(reduced_densities) > 1 
        # choose the root with the lowest free energy
        # get chemical potentials
        μ_i_at_density = [sanchez_lacombe_chemical_potentials(components, kij, φ, r_i0, r_i, pure_characteristic_volume_i, ρ̃ , temperature) for ρ̃ in reduced_densities]
        # get mole fractions
        # evaluate ∑(μᵢ * xᵢ) at each density 
        free_energy_at_density = [sum(μ_i .*  mole_fractions) for μ_i in μ_i_at_density]
        # pick the lowest value
        return reduced_densities[argmin(free_energy_at_density)]
    end
    
    # todo remove this and return a default, see if tests pass)
    if length(reduced_densities) == 0  # roots couldn't do it quickly, we will have to iterate (this actually probably means the situation is not realistic, 
        iter = 0
        reduced_density = 0
        const_sum = 1 - sum(φ ./ r_i)
        while iter <= maxiters
            new_reduced_density = 1 - exp(-(reduced_density^2 + reduced_pressure)/reduced_temperature - reduced_density*const_sum)
            if abs(new_reduced_density - reduced_density) < tol
                return new_reduced_density
            else
                reduced_density = new_reduced_density
            end
            iter+=1
        end
        throw(ErrorException("Max iters exceeded: " * string(maxiters)))
    
    end

    return reduced_densities[1]

end

function sanchez_lacombe_chemical_potentials(components::AbstractVector{<:SanchezLacombeParameters}, kij, mass_fractions, density, temperature)
    mixed_characteristic_density = sanchez_lacombe_mixed_characteristic_density(components, mass_fractions)
    reduced_density = density / mixed_characteristic_density
    φ = sanchez_lacombe_close_packed_volume_fractions(components, mass_fractions)
    r_i0 = sanchez_lacombe_ri0.(components)
    pure_characteristic_volume_i = sanchez_lacombe_pure_characteristic_volume.(components)
    p★ = sanchez_lacombe_mixed_characteristic_pressure(components, φ, kij)
    t★ = sanchez_lacombe_mixed_characteristic_temperature(components, φ, p★)
    v★ = sanchez_lacombe_mixed_characteristic_volume(t★, p★)
    r_i = sanchez_lacombe_ri(r_i0, pure_characteristic_volume_i, v★)
    return sanchez_lacombe_chemical_potentials(components, kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, temperature)
end

function sanchez_lacombe_chemical_potentials(components::AbstractVector{<:SanchezLacombeParameters}, kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, temperature)
    # J/mol
    p★i = characteristic_pressure.(components)
    Δpij = sanchez_lacombe_Δpij(components, kij)
    
    # if safeguards
    #     if reduced_density > 1; 
    #         println("ERROR 1"); 
    #         return zeros(length(components)); 
    #     end
    #     if (reduced_density) < 0; 
    #         println("ERROR 2"); 
    #         return zeros(length(components)); 
    #     end
    #     if minimum(φ) < 0; 
    #         println("ERROR 3"); 
    #         return zeros(length(components)); 
    #     end
    # end
    
    summation_terms = [sum([φ[jdx] * (p★i[jdx] - Δpij[idx, jdx]) for jdx in eachindex(φ, p★i)]) for idx in eachindex(φ, p★i)]
    term_1 = log.(reduced_density .* φ)
    term_2 = -log(1 - reduced_density) .* (r_i0 .+ (r_i .- r_i0) ./ reduced_density)
    term_3 = -r_i .+ 1
    term_4 = -reduced_density .* r_i0 .* (pure_characteristic_volume_i .* (p★i .+ summation_terms)) ./ (R_MPA_L_K_MOL * temperature) 
    result = (term_1 .+ term_2 .+ term_3 .+ term_4) * R_J_MOL_K * temperature
    # result = sanchez_lacombe_chemical_potentials(
    #     sanchez_lacombe_activities(components, kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, temperature),
    #     temperature)
    return result
end

function sanchez_lacombe_chemical_potentials(activities::AbstractVector{<:Number}, temperature::Number)
    return log.(activities) * R_J_MOL_K * temperature
end

function sanchez_lacombe_activities(chemical_potentials, temperature)
    return exp.(chemical_potentials / (R_J_MOL_K * temperature)) 
end

function sanchez_lacombe_activities(components::AbstractVector{<:SanchezLacombeParameters}, kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, temperature)
    # J/mol
    if (1 - reduced_density) < 0; return zeros(length(components)); end
    if (reduced_density) < 0; return zeros(length(components)); end
    if minimum(φ) < 0; return zeros(length(components)); end

    p★ = characteristic_pressure.(components)
    Δpij = sanchez_lacombe_Δpij(components, kij)
    
    summation_terms = [sum([φ[jdx] * (p★[jdx] - Δpij[idx, jdx]) for jdx in eachindex(φ, p★)]) for idx in eachindex(φ, p★)]
    term_1 = reduced_density * φ
    term_2 = (1 - reduced_density) .^ -(r_i0 .+ (r_i .- r_i0) ./ reduced_density)
    term_3 = (-r_i) .+ 1
    term_4 = -reduced_density .* r_i0 .* pure_characteristic_volume_i ./ (R_MPA_L_K_MOL * temperature) .* (p★ .+ summation_terms)
    return (term_1 .* term_2 .* exp.(term_3 .+ term_4)) 
    
    # return exp.(sanchez_lacombe_chemical_potentials(components, kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, temperature) / (R_J_MOL_K * temperature))
end

function sanchez_lacombe_ri0(chemical::SanchezLacombeParameters)
    # number of lattice cells occupied by a molecule of pure component i
    # MPa * g/mol / ((MPa * cm^3 /mol /K) * K * g/cm^33)
    ri0 = characteristic_pressure(chemical) * molecular_weight(chemical) / (R_MPA_CM3_K_MOL * characteristic_temperature(chemical) * characteristic_density(chemical))
    return ri0  
end

function sanchez_lacombe_pure_characteristic_volume(chemical::SanchezLacombeParameters)
    # MPa L / K*Mol * K / MPa = L/Mol
    return R_MPA_L_K_MOL * characteristic_temperature(chemical) / characteristic_pressure(chemical)
end

function sanchez_lacombe_mixed_characteristic_temperature(components::AbstractVector{<:SanchezLacombeParameters},  φ::AbstractVector{<:Number}, p_star)
    return p_star / sum(characteristic_pressure.(components) .* φ ./ characteristic_temperature.(components))
end

function sanchez_lacombe_Δpij(components::AbstractVector{<:SanchezLacombeParameters}, kij)
    p = characteristic_pressure.(components)
    return [p[idx] + p[jdx] - 2 * (1-kij[idx, jdx]) * sqrt(p[idx] * p[jdx]) for idx in eachindex(p), jdx in eachindex(p)]
end

function sanchez_lacombe_mixed_characteristic_pressure(components::AbstractVector{<:SanchezLacombeParameters}, φ::AbstractVector{<:Number}, kij)
    p = characteristic_pressure.(components)
    p★ = 0
    interaction_effects = 0
    
    Δpij = sanchez_lacombe_Δpij(components, kij)
    # for idx in eachindex(p, φ)
    #     @show p★ += φ[idx] * p[idx] - 0.5 * sum([φ[idx]*φ[jdx]*(p[idx] + p[jdx] - 2*(1-kij[idx, jdx])*sqrt(p[idx]*p[jdx])) for jdx in eachindex(p, φ)])
    # end  # delete this once confirming accuracy of model
    
    for idx in eachindex(p, φ)
        p★ += φ[idx] * p[idx]
    end

    for idx in eachindex(p, φ), jdx in eachindex(p, φ)
        interaction_effects += φ[idx] * φ[jdx] * Δpij[idx, jdx]
    end


    return p★ - 0.5 * interaction_effects
end

function sanchez_lacombe_mixed_characteristic_density(components::AbstractVector{<:SanchezLacombeParameters}, mass_fractions::AbstractVector{<:Number})
    inverse_rho_star = sum(mass_fractions ./ characteristic_density.(components)) # 1 / (g/cm^3)
    mixed_char_dens = 1 / inverse_rho_star
    return mixed_char_dens  
end

function sanchez_lacombe_close_packed_volume_fractions(components::AbstractVector{<:SanchezLacombeParameters}, mass_fractions::AbstractVector{<:Number})
    mass_frac_over_char_dens = mass_fractions ./ characteristic_density.(components)
    return mass_frac_over_char_dens ./ sum(mass_frac_over_char_dens)
end  

function sanchez_lacombe_mixed_characteristic_volume(t_star::Number, p_star::Number)
    # L/mol
    return t_star * R_MPA_L_K_MOL / p_star
end

function sanchez_lacombe_ri(ri0::Number, pure_characteristic_volume::Number, mixed_characteristic_volume::Number)
    return ri0 * pure_characteristic_volume / mixed_characteristic_volume
end


# stuff to eventually replace the EOS system above, but for now will exist concurrently until everything is moved over

struct SanchezLacombeModel{CPT, KIJ_T, N, L}
    components::SVector{N, CPT}
    kij::SMatrix{N, N, KIJ_T, L}
end
function SanchezLacombeModel(components::Vector{<:SanchezLacombeParameters}, kij::Matrix{<:Number})
    n = length(components)
    _kij = SMatrix{n, n}(kij)
    _components = SVector{n}(components)
    return SanchezLacombeModel(_components, _kij)
end
function SanchezLacombeModel(components::AbstractVector{<:SanchezLacombeParameters})
    return SanchezLacombeModel(components, initmatrix(components))
end
function SanchezLacombeModel(components::AbstractVector{<:String})
    kij = get_kij_matrix(modeltype, components)
    component_parameters = ChemicalParameters(components)
    return SanchezLacombeModel(component_parameters, kij)
end
function SanchezLacombeModel(component::ChemicalParameters)
    return SanchezLacombeModel([component])
end
function SanchezLacombeModel(component::String)
    return SanchezLacombeModel([component])
end
# function compressibility_factor(model::SanchezLacombeModel, p, t, mole_fractions=[1])
# end

# function VT_compressibility_factor(model::SanchezLacombeModel, v, t, mole_fractions=[1])
# end

function MembraneBase.pressure(model::SanchezLacombeModel, v, t, mole_fractions=[1])
    density = molar_volume_to_density(v, mole_fractions, molecular_weight.(model.components))
    mass_fractions = mole_fractions_to_mass_fractions(mole_fractions, molecular_weight.(model.components))
    mixed_characteristic_density = sanchez_lacombe_mixed_characteristic_density(model.components, mass_fractions)
    φ = sanchez_lacombe_close_packed_volume_fractions(model.components, mass_fractions)
    r_i0 = sanchez_lacombe_ri0.(model.components)
    pure_characteristic_volume_i = sanchez_lacombe_pure_characteristic_volume.(model.components)
    p★ = sanchez_lacombe_mixed_characteristic_pressure(model.components, φ, model.kij)
    t★ = sanchez_lacombe_mixed_characteristic_temperature(model.components, φ, p★)
    v★ = sanchez_lacombe_mixed_characteristic_volume(t★, p★)
    r_i = sanchez_lacombe_ri.(r_i0, pure_characteristic_volume_i, v★)
    reduced_density = density / mixed_characteristic_density
    pressure_mpa = sanchez_lacombe_reduced_pressure(
        t/t★, 
        reduced_density, 
        φ, r_i) * p★
        return pressure_mpa
end

function volume(model::SanchezLacombeModel, p, t, mole_fractions=[1])
    mass_fractions = mole_fractions_to_mass_fractions(mole_fractions, molecular_weight.(model.components))
    molecular_weights = molecular_weight.(model.components)
    mixed_characteristic_density = sanchez_lacombe_mixed_characteristic_density(model.components, mass_fractions)
    φ = sanchez_lacombe_close_packed_volume_fractions(model.components, mass_fractions)
    r_i0 = sanchez_lacombe_ri0.(model.components)
    pure_characteristic_volume_i = sanchez_lacombe_pure_characteristic_volume.(model.components)
    p★ = sanchez_lacombe_mixed_characteristic_pressure(model.components, φ, model.kij)
    t★ = sanchez_lacombe_mixed_characteristic_temperature(model.components, φ, p★)
    v★ = sanchez_lacombe_mixed_characteristic_volume(t★, p★)
    r_i = sanchez_lacombe_ri.(r_i0, pure_characteristic_volume_i, v★)
    reduced_temperature = t / t★
    reduced_pressure = p / p★
    reduced_density = sanchez_lacombe_reduced_density(model.components, model.kij, mole_fractions, t, reduced_temperature, reduced_pressure, φ, r_i0, r_i, pure_characteristic_volume_i)
    density_g_cm3 = reduced_density * mixed_characteristic_density
    volume = density_to_molar_volume(density_g_cm3, mole_fractions, molecular_weights)
    return volume
end

# function density(model::)

# function fugacity(model::SanchezLacombeModel, p, t, mole_fractions=[1])
# end

function chemical_potential(model::SanchezLacombeModel, p, t, mole_fractions=[1])
    v = volume(model, p, t, mole_fractions)
    return VT_chemical_potential(model, v, t, mole_fractions)
end

function VT_chemical_potential(model::SanchezLacombeModel, v, t, mole_fractions=[1])
    ρ = molar_volume_to_density(v, mole_fractions, molecular_weight.(model.components))
    ω = mole_fractions_to_mass_fractions(mole_fractions, molecular_weight.(model.components))
    return ρTω_chemical_potential(model, ρ, t, ω)
end

function ρTω_chemical_potential(model::SanchezLacombeModel, ρ, t, mass_fractions=[1]) 
    mixed_characteristic_density = sanchez_lacombe_mixed_characteristic_density(model.components, mass_fractions)
    reduced_density = ρ / mixed_characteristic_density

    φ = sanchez_lacombe_close_packed_volume_fractions(model.components, mass_fractions)
    
    p★ = sanchez_lacombe_mixed_characteristic_pressure(model.components, φ, model.kij)
    t★ = sanchez_lacombe_mixed_characteristic_temperature(model.components, φ, p★)
    v★ = sanchez_lacombe_mixed_characteristic_volume(t★, p★)
    
    pure_characteristic_volume_i = sanchez_lacombe_pure_characteristic_volume.(model.components)
    r_i0 = sanchez_lacombe_ri0.(model.components)
    r_i = sanchez_lacombe_ri.(r_i0, pure_characteristic_volume_i, v★)

    return sanchez_lacombe_chemical_potentials(model.components, model.kij, φ, r_i0, r_i, pure_characteristic_volume_i, reduced_density, t) 
end

function activity(model::SanchezLacombeModel, p, t, mole_fractions=[1])
    v = volume(model, p, t, mole_fractions)
    ρ = molar_volume_to_density(v, mole_fractions, molecular_weight.(model.components))
    ω = mole_fractions_to_mass_fractions(mole_fractions, molecular_weight.(model.components))
    return ρTω_activity(model, ρ, t, ω)
end

function ρTω_activity(model::SanchezLacombeModel, ρ, t, mass_fractions=[1])
    chemical_potentials = ρTω_chemical_potential(model, ρ, t, mass_fractions) 
    return sanchez_lacombe_activities(chemical_potentials, t)
end

function density_upper_bound(model::SanchezLacombeModel, mass_fractions::AbstractVector)
    return sanchez_lacombe_mixed_characteristic_density(model.components, mass_fractions)
end

