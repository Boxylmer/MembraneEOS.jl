# Make the Clapeyron.jl library compatible and extend its functionality a bit to the PolymerMembranes world

# struct Cubic{T,I,C} <: EoSModel
#     components::Vector{String}
#     a::PairParam{T}
#     b::PairParam{T}
#     kij::PairParam{T}
#     idealmodel::I
#     type::C
#   end


# Ωa,Ωb = 0.457235,0.077796 #Peng Robinson
# components = ["a","b","c","d"]
# k = PairParam("Kij",components,rand(5,5))
# pc = rand(5)
# Tc = rand(5)
# R = 0.0821 # ???
# ai = SingleParam("a",components, @. R^2*Tc^2/pc)
# bi = SingleParam("b",components, @. Ωb*R*Tc/pc)
# aij = epsilon_LorentzBerthelot(ai ,k)
# bij = sigma_LorentzBerthelot(bi)

# Clapeyron.split_model(bij, [[1,2,3]])


# model = eos_model(gss)
# v_fm1_s = Clapeyron.volume.(model, gss.sampling_chamber_final_pressures[step - 1],t=gss.temperature) 

