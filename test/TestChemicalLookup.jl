# --- test_chemical_parameters()

# look up parameters from the database with a string
methane = MembraneEOS.ChemicalParameters("CH4")

# if we define any values manually, we assume that the user isn't wanting to do database lookups
manual_methane = MembraneEOS.ChemicalParameters("CH4", critical_pressure=100.0)
@test ismissing(manual_methane.critical_temperature)

# otherwise, we can ask the constructor to go get missing values from the database, 
# (also the name must come first, if it's going to be specified)
manual_methane = MembraneEOS.ChemicalParameters("CH4", critical_pressure=100.0, fill_missing_with_database=true)
@test !ismissing(manual_methane.critical_temperature)
# but keep anything manually specified
@test methane.critical_pressure != manual_methane.critical_pressure

# when we define manual values, everything else is of type Missing
manual_methane = MembraneEOS.ChemicalParameters(critical_pressure=100.0)
@test ismissing(manual_methane.name)



# --- test_kij_database_lookups()

# values should be able to be looked up agnostic of order
co2_ch4 = MembraneEOS.UnorderedChemicalPair("CO2", "CH4")
ch4_co2 = MembraneEOS.UnorderedChemicalPair("CH4", "CO2")
@test MembraneEOS.PRKijLookup[co2_ch4] == MembraneEOS.PRKijLookup[ch4_co2]

# we can look up values by name in any equation of state with a database
co2_ch4_kij = MembraneEOS.get_kij(MembraneEOS.PRKijLookup, "CO2", "CH4")
@test !ismissing(co2_ch4_kij)
# and it should also be order-agnostic
ch4_co2_kij = MembraneEOS.get_kij(MembraneEOS.PRKijLookup, "CH4", "CO2")
@test co2_ch4_kij == ch4_co2_kij

# if we want to make interaction tables manually, we can load the specific tables and some components
# this is automatically handled in the actual EOS functions (e.g., PengRobinson) though. 
co2_ch4_missing = MembraneEOS.get_kij_matrix(MembraneEOS.PRKijLookup, ["CH4", "CO2", "Something not in the database"]; missing_value=missing, ideal_value=10)
@test ismissing(co2_ch4_missing[1, 3])
@test co2_ch4_missing[1, 1] == 10
@test co2_ch4_missing[3, 3] == 10

