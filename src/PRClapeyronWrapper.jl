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


# struct PRParam2{T} <: EoSParam
#     a::PairParam{T}
#     b::PairParam{T}
#     Tc::SingleParam{T}
#     Pc::SingleParam{T}
#     Mw::SingleParam{T}
# end

# struct PR2{T <: IdealModel,α,c,γ,K} <:PRModel
#     components::Array{String,1}
#     icomponents::UnitRange{Int}
#     alpha::α
#     mixing::γ
#     translation::c
#     params::PRParam2{K}
#     idealmodel::T
#     references::Array{String,1}
# end

# Clapeyron.@registermodel PR2

vdW1fRule() = Clapeyron.vdW1fRule(Clapeyron.vdW1fRuleParam())
NoTranslation() = NoTranslation(Clapeyron.NoTranslationParam())

function cubic_maker(cubic,Tc,Pc,w,Kij)
	icomponents = 1:length(Tc)
	components = [string(i) for i in icomponents ]
	acentricfactor = SingleParam("acentric factor",components,w)
	if cubic === PR
		alpha = PRAlpha(Clapeyron.PRAlphaParam(acentricfactor))
	elseif cubic == SRK
		alpha = SoaveAlpha(Clapeyron.SoaveAlphaParam(acentricfactor))
	elseif cubic == RK
		alpha = RKAlpha()
	else
		error("input cubic not valid")
	end
	mixing = vdW1fRule()
	tc = SingleParam("Tc",components,Tc)
	pc = SingleParam("Pc",components,Pc .* 1e-5)
	Mw = SingleParam("Mw",components,Pc .* 0.0) #not provided
	kij = PairParam("kij",components,Kij)
	a,b = Clapeyron.ab_premixing(cubic,mixing,tc,pc,kij)
	if cubic === PR
		params = Clapeyron.PRParam(a,b,tc,pc,Mw)
	elseif cubic == SRK || cubic == RK
		params =  Clapeyron.RKParam(a,b,tc,pc,Mw)
	else
		error("input cubic not valid")
	end
	translation = Clapeyron.NoTranslation(components)
	cubic(components,icomponents,alpha,mixing,translation,params,BasicIdeal(),String[])
end


using LinearAlgebra
ident = Matrix{Float64}(I, 1, 1)

cubic_maker(PR, [120.], [34.], [0.012], ident)

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
