# Make the Clapeyron.jl library compatible and extend its functionality a bit to the PolymerMembranes world

# struct Cubic{T,I,C} <: EoSModel
#     components::Vector{String}
#     a::PairParam{T}
#     b::PairParam{T}
#     kij::PairParam{T}
#     idealmodel::I
#     type::C
#   end
using Clapeyron


# PRParam2 and PR2 are more general than PR in that these structs can store measurement types
struct PRParam2{AT, BT, TCT, PCT, MWT} <: EoSParam
    a::PairParam{AT}
    b::PairParam{BT}
    Tc::SingleParam{TCT}
    Pc::SingleParam{PCT}
    Mw::SingleParam{MWT}
	#PRParam2(a, b, Tc, Pc, Mw) = PRParam2(promote(a, b)..., promote(Tc, Pc, Mw)...)
end

struct PR2{T <: IdealModel,α,c,γ,K} <: Clapeyron.PRModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    translation::c
    params::PRParam2{K}
    idealmodel::T
    references::Array{String,1}
end


Base.length(model::PR2) = length(model.icomponents)
Base.show(io::IO,model::PR2) = Clapeyron.eosshow(io,model)
# Base.show(io::IO,::MIME"text/plain",model::PR2) = Clapeyron.eosshow(io,mime,model)
Base.show(io::IO, mime::MIME"text/plain", model::PR2) = Clapeyron.eosshow(io,mime,model)
Clapeyron.molecular_weight(model::PR2,z=1.0) = Clapeyron.comp_molecular_weight(Clapeyron.mw(model),z)
# Clapeyron.@registermodel PR2

vdW1fRule() = Clapeyron.vdW1fRule(Clapeyron.vdW1fRuleParam())
NoTranslation() = NoTranslation(Clapeyron.NoTranslationParam())

function cubic_maker(cubic,Tc_k,Pc_atm,w,Kij)
	icomponents = 1:length(Tc_k)
	components = [string(i) for i in icomponents]
	acentricfactor = SingleParam("acentric factor",components,w)
	if cubic === PR
		alpha = PRAlpha(Clapeyron.PRAlphaParam(acentricfactor))
		eosmodel = PR2
	elseif cubic == SRK
		alpha = SoaveAlpha(Clapeyron.SoaveAlphaParam(acentricfactor))
		eosmodel = SRK
	elseif cubic == RK
		alpha = RKAlpha()
		eosmodel = RK
	else
		error("input cubic not valid")
	end
	mixing = vdW1fRule()
	tc = SingleParam("Tc", components, Tc_k)
	pc = SingleParam("Pc", components, Pc_atm .* 101325)
	Mw = SingleParam("Mw", components, Pc_atm .* 0.0) #not provided
	kij = PairParam("kij", components, Kij)
	a,b = Clapeyron.ab_premixing(cubic,mixing,tc,pc,kij)
	if cubic === PR
		params = PRParam2(a,b,tc,pc,Mw)
	elseif cubic == SRK || cubic == RK
		params =  Clapeyron.RKParam(a,b,tc,pc,Mw)
	else
		error("input cubic not valid")
	end
	translation = Clapeyron.NoTranslation(components)
	eosmodel(components,icomponents,alpha,mixing,translation,params,BasicIdeal(),String[])
	# cubic(components,icomponents,alpha,mixing,translation,params,BasicIdeal(),0.0,String[])
	# eosmodel(components,icomponents,alpha,mixing,translation,params,BasicIdeal(),String[])  # For clapeyron 0.3.4
end


using LinearAlgebra
ident = Matrix{Float64}(I, 1, 1)

co2model = cubic_maker(PR, [304.13], [72.8], [0.239], ident)
volume(co2model, 101325, 308.15)

using Measurements
co2model = cubic_maker(PR, [304.13 ± 0.1], [72.8 ± 10], [0.239 ± 0.2], ident)
volume(co2model, 101325 ± 100, 308.15 ± 34)

# # wrap individual EOS that we intend to use.
# function PR(pc_atm::Number, tc_k::Number, ω::Number)
#     # return PR(...) from the clapeyron library
# end

# function PR(chemical::String)
# 	# look up using internal ChemicalParameters and return a Clapeyron struct
# end
# function PR(chemical::AbstractVector{String})
#     # look critical param CSV -> ChemicalParameters 
#     # look up from binary interaction CSV -> UnorderedChemicalPair
#     # use both of these to call the constructor below
# end
# function PR(pc_atm::AbstractVector, tc_k::AbstractVector, omega::AbstractVector, KIJ_matrix=nothing)
	
# end
