"""
    walkindices(subset, image; region = CartesianIndices(image))

Return indices to walk `image` with size of `subset`.

```jldoctest
julia> image = rand(4,4);

julia> subset = rand(2,2);

julia> walkindices(subset, image)
3×3 Array{CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},2}:
 [CartesianIndex(1, 1) CartesianIndex(1, 2); CartesianIndex(2, 1) CartesianIndex(2, 2)]  …  [CartesianIndex(1, 3) CartesianIndex(1, 4); CartesianIndex(2, 3) CartesianIndex(2, 4)]
 [CartesianIndex(2, 1) CartesianIndex(2, 2); CartesianIndex(3, 1) CartesianIndex(3, 2)]     [CartesianIndex(2, 3) CartesianIndex(2, 4); CartesianIndex(3, 3) CartesianIndex(3, 4)]
 [CartesianIndex(3, 1) CartesianIndex(3, 2); CartesianIndex(4, 1) CartesianIndex(4, 2)]     [CartesianIndex(3, 3) CartesianIndex(3, 4); CartesianIndex(4, 3) CartesianIndex(4, 4)]
```
"""
function walkindices(subset::AbstractArray, image::AbstractArray; region::CartesianIndices = CartesianIndices(image))
    checkbounds(CartesianIndices(image), region)
    checkbounds(region, CartesianIndices(subset))
    origins = first(region):first(region)+CartesianIndex(size(region) .- size(subset))
    map(origins) do I
        I:I+CartesianIndex(size(subset).-1)
    end
end

"""
    neighborindices(subset::CartesianIndices, image, npixels::Int)

Return `npixels` outer indices around `subset`.
Violated indices in `image` are cut automatically.
This is useful to give `region` in [`coarse_search`](@ref).

```jldoctest
julia> image = rand(10,10);

julia> neighborindices(CartesianIndices((4:6, 3:6)), image, 2)
7×8 CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}}:
 CartesianIndex(2, 1)  CartesianIndex(2, 2)  …  CartesianIndex(2, 8)
 CartesianIndex(3, 1)  CartesianIndex(3, 2)     CartesianIndex(3, 8)
 CartesianIndex(4, 1)  CartesianIndex(4, 2)     CartesianIndex(4, 8)
 CartesianIndex(5, 1)  CartesianIndex(5, 2)     CartesianIndex(5, 8)
 CartesianIndex(6, 1)  CartesianIndex(6, 2)     CartesianIndex(6, 8)
 CartesianIndex(7, 1)  CartesianIndex(7, 2)  …  CartesianIndex(7, 8)
 CartesianIndex(8, 1)  CartesianIndex(8, 2)     CartesianIndex(8, 8)
```
"""
function neighborindices(subset::CartesianIndices, image::AbstractArray, npixels::Int)
    start = Tuple(first(subset)) .- npixels
    stop = Tuple(last(subset)) .+ npixels
    newstart = clamp.(start, 1, size(image))
    newstop = clamp.(stop, 1, size(image))
    CartesianIndex(newstart):CartesianIndex(newstop)
end

function neighborindices(point::CartesianIndex, image::AbstractArray, npixels::Int)
    neighborindices(point:point, image, npixels)
end

function testimage(name::String)
    if splitext(name)[1] == "buffalo"
        return load(joinpath(dirname(@__FILE__), "../images/buffalo.tif"))
    end
    throw(ArgumentError("test image $name is not exist"))
end
