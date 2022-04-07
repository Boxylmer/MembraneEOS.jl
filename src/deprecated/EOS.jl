

# helper functions for general EOS use
function compute_ideal_partial_pressures(pressure, mole_fractions)
    N = length(mole_fractions)
    partial_pressures = SVector{N}(pressure .* mole_fractions)
    return partial_pressures
end
