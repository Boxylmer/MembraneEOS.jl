
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
fugacity(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.fugacity_coefficient(model, p_mpa * 1e6, t_k, z) .* p_mpa

compressibility_factor(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.compressibility_factor(model, p_mpa * 1e6, t_k, z)

VT_compressibility_factor(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.VT_compressibility_factor(model, v_l_mol / 1000, t_k, z)

"Mass density in g/cm^3"
function mass_density(model::Clapeyron.SL, p_mpa, t_k, z=[1])
    return Clapeyron.mass_density(model, p_mpa * 1.0e6, t_k, z) * 0.001
end

"Chemical potential in J/mol"
chemical_potential(model::Clapeyron.SL, p_mpa, t_k, z=[1]) = Clapeyron.chemical_potential(model, p_mpa * 1.0e6, t_k, z)

"Chemical potential in J/mol"
VT_chemical_potential(model::Clapeyron.SL, v_l_mol, t_k, z=[1]) = Clapeyron.VT_chemical_potential(model, v_l_mol / 1000, t_k, z)
