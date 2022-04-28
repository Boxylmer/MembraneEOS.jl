# Sanchez Lacombe
struct SanchezLacombe end
function get_kij(::SanchezLacombe, component_1::String, component_2::String)
    if component_1 == component_2 return 0.0 end
    return get_kij(SLKijLookup, component_1, component_2) 
end
function get_kij_matrix(::SanchezLacombe, components::AbstractVector{<:String})
    return get_kij_matrix(SLKijLookup, components)
end 

### things for new SL mixing rule

struct SLKRule2 <: Clapeyron.SLMixingRule
    components::Vector{String}
    k::PairParam{Float64}
end


function Clapeyron.sl_mix(unmixed_vol,unmixed_epsilon, mixmodel::SLKRule2)
    #dont mind the function names, it performs the correct mixing
    premixed_vol = Clapeyron.epsilon_LorentzBerthelot(unmixed_vol)
    premixed_epsilon = Clapeyron.sigma_LorentzBerthelot(unmixed_epsilon)
    return premixed_vol, premixed_epsilon
end

function Clapeyron.a_res(model::Clapeyron.SanchezLacombe,V,T,z=SA[1.0])
    Σz = sum(z)     
    r = model.params.segment.values
    mixing = model.mixing
    r̄ = dot(z,r)
    r̄ = r̄/Σz
    v_r,ε_r = Clapeyron.mix_vε(model,V,T,z,mixing,r̄,Σz)
    @show v_r, ε_r 
    v = V/Σz
    ρ̃ = r̄*v_r/v
    T̃ = Clapeyron.R̄*T/ε_r
    _1 = one(V+T+first(z))
    return r̄*(-ρ̃ /T̃ + (_1/ρ̃  - _1)*log1p(-ρ̃ )+_1)
end

function Clapeyron.mix_vε(model::Clapeyron.SL,V,T,z,mix::SLKRule2,r̄,Σz = sum(z)) 
    v = model.params.vol.diagvalues  # m^3/mol
    ε = model.params.epsilon.diagvalues  # J/mol 
    isone(length(z)) && return (only(v),only(ε))
    r =  model.params.segment.values  # unitless
    k = mix.k.values   # unitless
    r̄inv = one(r̄)/r̄
    ϕ = @. r* z* r̄inv/Σz   # unitless
    p = ε ./ v  # J/mol  /  m^3/mol  -> J/m^3 -> Pa*m^3 / m^3 -> Pa
    Δpij = [p[idx] + p[jdx] - 2 * (1-k[idx, jdx]) * sqrt(p[idx] * p[jdx]) for idx in eachindex(p), jdx in eachindex(p)]
    
    # p★_ideal = 0
    # for idx in eachindex(p, ϕ)
    #     p★_ideal += ϕ[idx] * p[idx]
    # end
    # interaction_effects = 0
    # for idx in eachindex(p, ϕ), jdx in eachindex(p, ϕ)
    #     interaction_effects += ϕ[idx] * ϕ[jdx] * Δpij[idx, jdx]
    # end

    # p★ = p★_ideal - 0.5 * interaction_effects   # square of the Hildebrand Solubility Parameter and cohesive energy density


    p★ = zero(eltype(z))
    for i in 1:length(z)
        ϕi = ϕ[i]
        p_i = p[i]
        p★ += ϕ[i]*p[i]
        for j in 1:i-1 # i != j
            p_j = p[j]
            ϕj = ϕ[j]
            Δpij = p_i + p_j - 2*(1 - k[i,j])*sqrt(p_i*p_j)
            p★ -= ϕi*ϕj*Δpij #0.5*2
        end
    end

    t = ε ./ Clapeyron.R̄   # J/mol / (m2 kg s-2 K-1 mol-1 = J/(molK))  --> K
    t★ = p★ / sum(p .* ϕ ./ t)  # Pa / (Pa/K) --> K
    v★ = t★ * Clapeyron.R̄ / p★  # K * (m2 kg s-2 K-1 mol-1) / Pa = m2 kg s-2 mol-1 * kg-1 m s2 = m3 / mol
    ε★ = t★ * Clapeyron.R̄  # K * J / (mol K) = J/mol
    return v★,ε★   # units: m^3/mol, J/mol 
end

### end of things for new SL mixing rule 


SL(component::String) = SL([component])
function SL(components::AbstractVector, KIJ_matrix = nothing)
    if isnothing(KIJ_matrix)
        KIJ_matrix = get_kij_matrix(SanchezLacombe(), components)
    end
    component_parameters = ChemicalParameters(components)
    t★ = characteristic_temperature.(component_parameters)
    p★ = characteristic_pressure.(component_parameters)
    ρ★ = characteristic_density.(component_parameters)
    mw = Vector{Float64}(molecular_weight.(component_parameters))
    return SL(p★, t★, ρ★, mw, KIJ_matrix)
end
SL(p★::Number, t★::Number, ρ★::Number, mw::Number) = SL([p★], [t★], [ρ★], [mw])
function SL(p★::AbstractVector, t★::AbstractVector, ρ★::AbstractVector, mw::AbstractVector, kij = zeros(length(mw),length(mw)))
    R = Clapeyron.R̄  # cm3*mpa / kmol
    icomponents = 1:length(p★)
    components = [string(i) for i in icomponents]
    _v★ = R .* t★./p★ .* 1e-6 #m3/mol
    _ε = p★ .* _v★ .* 1e6  # cm3 * MPa / mol -> J/mol
    _r = mw ./ (ρ★ .* _v★ .* 1e6)  # g/mol / (cm3/mol * g/cm3) -> no units
    v★ = SingleParam("vol", components, _v★)
    ε = SingleParam("epsilon", components, _ε)
    r = SingleParam("segment", components, _r)
    mwparam = SingleParam("Mw", components, mw)
    kij = PairParam("kij", components, kij .* 1.0)
    mixing = Clapeyron.SLKRule2(components, kij)
    ideal = Clapeyron.init_model(Clapeyron.BasicIdeal, components, String[], false)
    premixed_vol, premixed_epsilon = Clapeyron.sl_mix(v★, ε, mixing)
    # mixing_new = MembraneEOS.SLKRule2(components, kij)
    # @show premixed_vol_new, premixed_epsilon_new = Clapeyron.sl_mix(v★, ε, mixing_new)
    packagedparams = Clapeyron.SanchezLacombeParam(mwparam, r, premixed_epsilon, premixed_vol)
    return Clapeyron.SL(components, icomponents, mixing, packagedparams, ideal, String[])
end

# functionality
"Pressure in MPa"
MembraneBase.pressure(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.pressure(model, v_l_mol * 1e-3, t_k, z) * 1e-6

"Volume in L/mol"
volume(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.volume(model, p_mpa * 1e6, t_k, z) * 1000

"Fugacity in MPa"
fugacity(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.fugacity_coefficient(model, p_mpa * 1e6, t_k, z) .* p_mpa .* z

molecular_weight(model::Clapeyron.SL) = Clapeyron.mw.(model)

compressibility_factor(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.compressibility_factor(model, p_mpa * 1e6, t_k, z)

VT_compressibility_factor(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.VT_compressibility_factor(model, v_l_mol / 1000, t_k, z)

"Mass density in g/cm^3"
function mass_density(model::Clapeyron.SL, p_mpa, t_k, z=[1])
    return Clapeyron.mass_density(model, p_mpa * 1.0e6, t_k, z) * 0.001  # kg/m3 * 1000g/kg * 1/1e6 m3/cm3
end

"Mass density in g/cm^3"
VT_mass_density(model::Clapeyron.SL, v_l_mol, z=[1]) = Clapeyron.mw(model) ./ v_l_mol .* z  ./ 1000   # g/mol / l/mol * l/cm3

"Chemical potential in J/mol"
chemical_potential(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.chemical_potential(model, p_mpa * 1.0e6, t_k, z)

"Residual chemical potential in J/mol"
chemical_potential_res(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.chemical_potential_res(model, p_mpa * 1.0e6, t_k, z)

"Chemical potential in J/mol"
VT_chemical_potential(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.VT_chemical_potential(model, v_l_mol / 1000, t_k, z)

# todo 

"Chemical potential in J/mol"
function ρTω_chemical_potential(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1])
    mw = Clapeyron.mw(model)  # g/cm3
    z = mass_fractions_to_mole_fractions(ω, mw) 
    v = density_to_molar_volume(ρ_g_cm3, z, mw) ./ 1000  # l/mol to m^3/mol
    μ = Clapeyron.VT_chemical_potential(model, v, t_k, z) 
    return μ 
end

"Chemical potential in J/mol"
function ρTω_chemical_potential_res(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1])
    mw = Clapeyron.mw(model)  # g/cm3
    z = mass_fractions_to_mole_fractions(ω, mw) 
    v = density_to_molar_volume(ρ_g_cm3, z, mw) ./ 1000  # l/mol to m^3/mol
    μ = Clapeyron.VT_chemical_potential_res(model, v, t_k, z) 
    return μ 
end

activity(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = exp.(chemical_potential(model, p_mpa, t_k, z) ./ (MembraneBase.R_J_MOL_K * t_k))

activity_res(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = exp.(chemical_potential_res(model, p_mpa, t_k, z) ./ (MembraneBase.R_J_MOL_K * t_k))

ρTω_activity(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1]) = exp.(ρTω_chemical_potential(model, ρ_g_cm3, t_k, ω) ./ (MembraneBase.R_J_MOL_K * t_k))

ρTω_activity_res(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1]) = exp.(ρTω_chemical_potential_res(model, ρ_g_cm3, t_k, ω) ./ (MembraneBase.R_J_MOL_K * t_k))

"Density upper bound in g/cm^3."
function density_upper_bound(model::Clapeyron.SL, ω=[1]) 
    mole_fracs = mass_fractions_to_mole_fractions(ω, molecular_weight(model))
    return molar_volume_to_density(
        Clapeyron.lb_volume(model, mole_fracs) * 1000, # m^3/mol -> l/mol
        mole_fracs,
        molecular_weight(model)
    )

end