"""
    zncc(image1, image2)

Perform zero-mean normalized cross-correlation.
"""
function zncc(A::AbstractArray{Float64}, B::AbstractArray{Float64})
    size(A) == size(B) || throw(DimensionMismatch("Dimensions must match."))
    Ā = mean(A)
    B̄ = mean(B)
    n = sum((Aᵢ - Ā) * (Bᵢ - B̄) for (Aᵢ, Bᵢ) in zip(A, B))
    d = sum((Aᵢ - Ā)^2 for Aᵢ in A) * sum((Bᵢ - B̄)^2 for Bᵢ in B)
    n / sqrt(d)
end

function zncc(A::AbstractArray{<: Gray}, B::AbstractArray{<: Gray})
    zncc(mappedarray(Float64, A), mappedarray(Float64, B))
end

function zncc(A::AbstractArray{<: RGB}, B::AbstractArray{<: RGB})
    zncc(mappedarray(Gray, A), mappedarray(Gray, B))
end

"""
    search(subset, image; region = CartesianIndices(image)) -> indices, R

Search `subset` in `image` using DIC.
Return the `indices` which has the highest correlation with `subset`.
Use `image[indices]` to get the found part of image.
The searching `region` (entire image by default) can also be specified
by `CartesianIndices` to reduce computations.

See also [`neighborindices`](@ref).

```jldoctest
julia> image = rand(10,10);

julia> subset = image[3:5, 2:3];

julia> search(subset, image)
(CartesianIndex{2}[CartesianIndex(3, 2) CartesianIndex(3, 3); CartesianIndex(4, 2) CartesianIndex(4, 3); CartesianIndex(5, 2) CartesianIndex(5, 3)], 1.0)
```
"""
function search(subset::AbstractArray, image::AbstractArray; region::CartesianIndices = CartesianIndices(image))
    Rs = map(walkindices(subset, image; region)) do inds
        zncc(view(image, inds), subset)
    end
    I = argmax(Rs)
    I:I+CartesianIndex(size(subset).-1), Rs[I]
end
