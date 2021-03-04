module DIC

using Reexport
@reexport using Images
using FileIO, ImageMagick, ImageIO # io
using ImageView

using
    Statistics,
    MappedArrays,
    Interpolations,
# automatic differentiations
    StaticArrays,
    ForwardDiff,
    DiffResults

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
    fine_search

include("utils.jl")
include("searching.jl")

end # module
