using MembraneEOS
using MembraneBase

co2_pure_model = SL(["CO2"])
p = 1; t = 273.15
z = [0.2, 0.8]

# pure system with only one gas component
@show v = volume(co2_pure_model, p, t)  # 2.136649093323808 L/mol
                      # Clapeyron returns 2.1366490933238094 L/mol --> normal a_res, Clapeyron mixing rule
                      # Clapeyron returns 2.1366490933238094 L/mol --> normal a_res, Membrane mixing rule
# mixed system with two regular gas components
co2_ch4_model = SL(["CO2", "CH4"], [0 0; 0 0.])
@show v = volume(co2_ch4_model, p, t, z)  # 2.216132223322739 L/mol
                        # Clapeyron returns 2.211124407739013 L/mol --> normal a_res, Clapeyron mixing rule
                        # Clapeyron returns 2.210750908462087 L/mol --> normal a_res, Membrane mixing rule

# mixed system with two regular gas components and nonideal interactions
co2_ch4_model_with_kij = SL(["CO2", "CH4"], [0 -0.1; 0 -0.1])
@show v = volume(co2_ch4_model_with_kij, p, t, z)  # 2.205729481036781 L/mol
                                 # Clapeyron returns 2.2013115585553185 L/mol --> normal a_res, Clapeyron mixing rule
                                 # Clapeyron returns 2.210750908462087 L/mol --> normal a_res, Membrane mixing rule

# pure system with only one polymer component
model_pdms = MembraneEOS.SL([302.0], [476.0], [1.104], [100000])
@show v = volume(model_pdms, p, t)  #  100.37221812334239 L/mol
                   # Clapeyron returns 100.37221812334239 L/mol --> normal a_res, Clapeyron mixing rule
                   # Clapeyron returns 100.37221812334239 L/mol --> normal a_res, Membrane mixing rule

# mixed system with a polymer and gas component
model_pdms_co2 = MembraneEOS.SL([302.0, 630.0], [476.0, 300.0], [1.104, 1.515], [100000, 44.0])
ω = [0.998, 0.002]
@show z = mass_fractions_to_mole_fractions(ω, molecular_weight(model_pdms_co2))  # [0.18003214274000456, 0.8199678572599954]
@show v = volume(model_pdms_co2, p, t, z)  # 18.108297158036713 L/mol
                         # Clapeyron returns 18.08932995502047  L/mol --> normal a_res, Clapeyron mixing rule
                         # Clapeyron returns 18.130532214647367 L/mol --> normal a_res, Membrane mixing rule
