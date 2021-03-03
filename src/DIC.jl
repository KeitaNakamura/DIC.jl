module DIC

using Reexport
@reexport using Images
using FileIO, ImageMagick, ImageIO # io
using ImageView

using Statistics, MappedArrays

export
# io/visualize
    load,
    save,
    imshow,
# utils
    walkindices,
    neighborindices,
# searching
    search

include("utils.jl")
include("searching.jl")

end # module
