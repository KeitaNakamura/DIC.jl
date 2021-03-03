"""
    zncc(image1, image2)

Perform zero-mean normalized cross-correlation.
"""
function zncc(A::AbstractArray{Float64}, B::AbstractArray{Float64})
    size(A) == size(B) || throw(DimensionMismatch("Dimensions must match."))

    # mean values
    Ā = 0.0
    B̄ = 0.0
    @inbounds @simd for i in eachindex(A, B)
        Ā += A[i]
        B̄ += B[i]
    end
    Ā /= length(A)
    B̄ /= length(B)

    # numerator/denominator
    n = 0.0
    d_A = 0.0
    d_B = 0.0
    @inbounds @simd for i in eachindex(A, B)
        Aᵢ = A[i]
        Bᵢ = B[i]
        n += (Aᵢ - Ā) * (Bᵢ - B̄)
        d_A += (Aᵢ - Ā)^2
        d_B += (Bᵢ - B̄)^2
    end

    n / sqrt(d_A * d_B)
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
    inds = walkindices(subset, image; region)
    Rs = similar(inds, Float64)
    Threads.@threads for i in eachindex(inds, Rs)
        @inbounds Rs[i] = zncc(view(image, inds[i]), subset)
    end
    I = argmax(Rs)
    I:I+CartesianIndex(size(subset).-1), Rs[I]
end
