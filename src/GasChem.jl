module GasChem

using EarthSciMLBase
using ModelingToolkit
using Catalyst
using Dates
using DynamicQuantities
using StaticArrays
using Interpolations
using ModelingToolkit:t

@register_unit molec 1
@register_unit mol_air 1u"mol"
@register_unit ppb 1u"mol/mol_air"
#TODO if @register_unit ppb 1e-9u"mol/mol_air" if @register_unit ppb 1e-9u"mol/mol_air", though it's physically correct, but this will result ModelingToolkit.ValidationError when coupling different models.

include("geoschem_ratelaws.jl")
include("geoschem_fullchem.jl")
include("SuperFast.jl")
include("Fast-JX.jl")
include("fastjx_couplings.jl")

end
