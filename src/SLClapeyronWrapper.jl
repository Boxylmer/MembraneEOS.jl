
using Clapeyron


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
    kij = PairParam("kij", components, kij)
    k1ij =  PairParam("k1", components, zeros(n,n))
    lij =  PairParam("l", components, zeros(n,n))
    mixing = SLk0k1lMixingRule(components, kij ,k1ij ,lij)
    ideal = Clapeyron.init_model(Clapeyron.BasicIdeal, components, String[], false)
    premixed_vol, premixed_epsilon = Clapeyron.sl_mix(v★, ε, mixing)
    packagedparams = Clapeyron.SanchezLacombeParam(mwparam, r, premixed_epsilon, premixed_vol)
    return SL(components, icomponents, mixing, packagedparams, ideal, String[])
end

# polycarbonate T* = 755, P* = 534, r0 = 1.275, mw = 1.00E+09 
using Measurements
model_co2 = SL([630.0], [300.0], [1.515], [44.0])
model_pdms = SL2([302.0], [476.0], [1.104], [1.00E+30])
@show mass_density(model_pdms, 101325.0, 273.15) * 0.001
chemical_potential(model_pdms, 101325.0, 273.15)