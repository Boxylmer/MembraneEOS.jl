# Sanchez Lacombe
struct SanchezLacombe end
function get_kij(::SanchezLacombe, component_1::String, component_2::String)
    if component_1 == component_2 return 0.0 end
    return get_kij(SLKijLookup, component_1, component_2) 
end
function get_kij_matrix(::SanchezLacombe, components::AbstractVector{<:String})
    return get_kij_matrix(SLKijLookup, components)
end 

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
    n = length(icomponents)
    components = [string(i) for i in icomponents]
    _v★ = R .* t★./p★ .* 1e-6 #m3/mol
    _ε = p★ .* _v★ .* 1e6  # cm3 * MPa / mol -> J/mol
    _r = mw ./ (ρ★ .* _v★ .* 1e6)  # g/mol / (cm3/mol * g/cm3) -> no units
    v★ = SingleParam("vol", components, _v★)
    ε = SingleParam("epsilon", components, _ε)
    r = SingleParam("segment", components, _r)
    mwparam = SingleParam("Mw", components, mw)
    kij = PairParam("kij", components, kij .* 1.0)
    k1ij =  PairParam("k1", components, zeros(Float64, n, n))
    lij =  PairParam("l", components, zeros(Float64, n, n))
    mixing = SLk0k1lMixingRule(components, kij, k1ij, lij)
    ideal = Clapeyron.init_model(Clapeyron.BasicIdeal, components, String[], false)
    premixed_vol, premixed_epsilon = Clapeyron.sl_mix(v★, ε, mixing)
    packagedparams = Clapeyron.SanchezLacombeParam(mwparam, r, premixed_epsilon, premixed_vol)
    return Clapeyron.SL(components, icomponents, mixing, packagedparams, ideal, String[])
end


# functionality
"Pressure in MPa"
pressure(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.pressure(model, v_l_mol * 1e-3, t_k, z) * 1e-6

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

ρTω_activity(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1]) = exp.(ρTω_chemical_potential(model, ρ_g_cm3, t_k, ω) ./ (MembraneBase.R_J_MOL_K * t_k))

ρTω_activity_res(model::Clapeyron.SL, ρ_g_cm3, t_k, ω=[1]) = exp.(ρTω_chemical_potential_res(model, ρ_g_cm3, t_k, ω) ./ (MembraneBase.R_J_MOL_K * t_k))