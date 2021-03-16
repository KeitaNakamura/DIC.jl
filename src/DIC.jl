module DIC

using Reexport
@reexport using Coordinates
using Colors
using FileIO, ImageMagick, ImageIO # io
using ImageView
using LinearAlgebra: tr, dot

using
    MappedArrays,
    Interpolations,
# automatic differentiations
    StaticArrays,
    ForwardDiff,
    DiffResults,
# progress bar
    ProgressMeter,
# plot
    GLMakie

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
    displacement_field,
    plot

const PixelIndices{dim} = AbstractArray{<: CartesianIndex{dim}}

include("utils.jl")
include("searching.jl")
include("displacement.jl")
include("Camera.jl")

end # module
