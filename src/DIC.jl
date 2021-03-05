module DIC

using Reexport
@reexport using Coordinates
using Colors
using FileIO, ImageMagick, ImageIO # io
using ImageView

using
    Statistics,
    MappedArrays,
    Interpolations,
# automatic differentiations
    StaticArrays,
    ForwardDiff,
    DiffResults,
# progress bar
    ProgressMeter

using Base: @_propagate_inbounds_meta

export
# io/visualize
    load,
    save,
    imshow,
# utils
    walkindices,
    neighborindices,
# searching
    coarse_search,
    fine_search,
# displacement
    strain,
    displacement_field

const PixelIndices{dim} = AbstractArray{<: CartesianIndex{dim}}

include("utils.jl")
include("searching.jl")
include("displacement.jl")

end # module
