abstract type EquationOfState end

"""
    EOS(EquationOfState, args...; kwargs...)

Construct an equation of state generically. `args` and `kwargs` are forwarded to the specifed equation of state. It is only advised to use this constructor when using the chemical parameter database.
"""
function EOS end


# need to be moved to a KIJ function 
function initmatrix(components::AbstractVector; initial_value = 0.0)
    matrix_size = length(components)
    Matrix(initial_value*I,matrix_size,matrix_size)
end

# function symmetricify()
# end

# EOS Accessors 

"Get the state compressibility `z`."
compressibility(state::EquationOfState) = state.z

"Get the state molar volume in **liters per mol**"
volume(state::EquationOfState) = state.volume

"Get the state pressure in **atmospheres**"
pressure(state::EquationOfState) = state.pressure

"Get the state temperature in **Kelvin**"
temperature(state::EquationOfState) = state.temperature

"Get the mole fractions of the state in the same order as the `components`"
mole_fractions(state::EquationOfState) = state.mole_fractions

"Get a vector of ChemicalParameters in the state, in the order they were originally specified"
components(state::EquationOfState) = state.components

"Get a vector of the names of the chemicals in the state (if present), in the order they were originally specified"
component_names(state::EquationOfState) = [component.name for component in state.components]

"Get the partial pressures of the state."
partial_pressures(state::EquationOfState) = state.partial_pressures

"Get the fugacities (MPa) of the state."
fugacities(state::EquationOfState) = state.fugacities


# helper functions for general EOS use
function compute_ideal_partial_pressures(pressure, mole_fractions)
    N = length(mole_fractions)
    partial_pressures = SVector{N}(pressure .* mole_fractions)
    return partial_pressures
end
# partial properties
function compute_partial_pressures(state::EquationOfState)
    function pressure_func(n_i)
        return pressure(
            EOS(state.model, state.components, n_i ./ sum(n_i), state.kij; v=state.volume * sum(n_i), t=state.temperature)
        )
    end
    # function get_dp_dn(molfracs...)
    #     # @show FiniteDifferences.finite_difference_derivative(pressure_func, molfracs)        
    # end
    @show ForwardDiff.gradient(pressure_func, state.mole_fractions)
    unnormalized_dp_dn = @show grad(central_fdm(5, 1), pressure_func, state.mole_fractions)[1]
    # unnormalized_dp_dn = get_dp_dn(...)
    pis = -unnormalized_dp_dn .* state.mole_fractions
    return pis
end
# function subset_kij(kijmatrix::AbstractMatrix, start_index, end_index)

# end
# need to take in what is effectively a triangular matrix or set of variables in a particular order


# series of methods that will eventually replace the stateful method currently employed

# function pressure(model::EquationOfState, v, t, mole_fractions=missing, kij=missing)
# end