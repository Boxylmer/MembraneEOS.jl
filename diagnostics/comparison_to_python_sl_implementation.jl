using MembraneEOS
using MembraneBase
using Plots

co2_model = SL([630], [300], [1.515], [44])
co2_tpbo_model = SL([474, 630], [900, 300], [1.6624, 1.515], [10e6, 44], [0 0; 0 0])  # a potentially optimal kij for tpbo-0.25/co2 is [0 -0.0356; -0.0356 0]

t = 308.15
dry_polymer_ρ = 1.393  # g/cm3

python_pressures_co2 = [0.105263158, 0.210526316, 0.315789474, 0.421052632, 0.526315789, 0.631578947, 0.736842105, 0.842105263, 0.947368421, 1.052631579, 1.157894737, 1.263157895, 1.368421053, 1.473684211, 1.578947368, 1.684210526, 1.789473684, 1.894736842, 2]
python_polphase_co2_ω = [0.027715355, 0.042083848, 0.051303895, 0.057906687, 0.062958964, 0.066999547, 0.070334648, 0.073153325, 0.075579719, 0.077699312, 0.07957324, 0.081246575, 0.082753451, 0.084120243, 0.085367731, 0.086512547, 0.08756819, 0.088545754, 0.089454455]
python_polphase_co2_μ = 1000 .* [-33.48841488, -31.72317255, -30.6950378, -29.968712, -29.40777782, -28.95147154, -28.56738029, -28.2361567, -27.94532156, -27.68635426, -27.45317714, -27.24130517, -27.04732909, -26.86860126, -26.70302237, -26.54890031, -26.40485116, -26.26972868, -26.14257288]
python_gasphase_co2_μ = 1000 .* [-33.4884142, -31.72317198, -30.69503721, -29.96871144, -29.40777724, -28.95147097, -28.56737971, -28.23615611, -27.94532096, -27.68635369, -27.4531772, -27.24130453, -27.04732844, -26.86860058, -26.7030217, -26.54889963, -26.40485047, -26.26972798, -26.14257218]
# unitless and j/mol respectively, chemical potential is configurational

# assuming no swelling, the polymer phase density is just dry_pol_ρ / ω_pol
polymer_phase_ρ = dry_polymer_ρ ./ (1 .- python_polphase_co2_ω)

julia_polphase_co2_μ = Vector{Float64}(undef, length(python_polphase_co2_ω))
julia_gasphase_co2_μ = Vector{Float64}(undef, length(python_polphase_co2_ω))


for i in eachindex(julia_polphase_co2_μ)
    mass_fracs_at_i = [1. - python_polphase_co2_ω[i], python_polphase_co2_ω[i]]
    julia_polphase_co2_μ[i] = ρTω_chemical_potential(co2_tpbo_model, polymer_phase_ρ[i], t, mass_fracs_at_i)[2]
    julia_gasphase_co2_μ[i] = chemical_potential(co2_model, python_pressures_co2[i], t, [1])[1]
end



comparison = plot(python_polphase_co2_ω, python_polphase_co2_μ, label="python", legend=:topleft)
plot!(comparison, python_polphase_co2_ω, julia_polphase_co2_μ, label="julia")
plot!(comparison, python_polphase_co2_ω, julia_gasphase_co2_μ, label="julia, gas phase")
savefig(comparison, joinpath(@__DIR__, "TPBO-0.25 and CO2 at 35C chemical potential comparison.png"))


