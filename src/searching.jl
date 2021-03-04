"""
    zncc(image1, image2)

Perform zero-mean normalized cross-correlation.
"""
function zncc(A::AbstractArray{T}, B::AbstractArray{U}) where {T <: Real, U <: Real}
    size(A) == size(B) || throw(DimensionMismatch("Dimensions must match."))

    # mean values
    Ā = zero(T)
    B̄ = zero(U)
    @inbounds @simd for i in eachindex(A, B)
        Ā += A[i]
        B̄ += B[i]
    end
    Ā /= length(A)
    B̄ /= length(B)

    # numerator/denominator
    n = zero(promote_type(T, U))
    d_A = zero(T)
    d_B = zero(U)
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
    coarse_search(subset, image; region = CartesianIndices(image)) -> indices, R

Perform coarse search `subset` in `image` using DIC.
Return the `indices` which has the highest correlation with `subset`.
Use `image[indices]` to get the found part of image.
The searching `region` (entire image by default) can also be specified
by `CartesianIndices` to reduce computations.

See also [`neighborindices`](@ref).

```jldoctest
julia> image = rand(10,10);

julia> subset = image[3:5, 2:3];

julia> coarse_search(subset, image)
(CartesianIndex{2}[CartesianIndex(3, 2) CartesianIndex(3, 3); CartesianIndex(4, 2) CartesianIndex(4, 3); CartesianIndex(5, 2) CartesianIndex(5, 3)], 1.0)
```
"""
function coarse_search(subset::AbstractArray, image::AbstractArray; region::CartesianIndices = CartesianIndices(image))
    inds = walkindices(subset, image; region)
    Cs = similar(inds, Float64)
    Threads.@threads for i in eachindex(inds, Cs)
        @inbounds Cs[i] = zncc(view(image, inds[i]), subset)
    end
    I = argmax(Cs)
    I:I+CartesianIndex(size(subset).-1), Cs[I]
end

# for 2D
solution_vector(::Val{2}) = zero(SVector{6, Float64})
function compute_correlation(subset::AbstractArray{<: Gray, 2}, image_itp::AbstractArray{<: Real, 2}, first_guess::CartesianIndices{2}, X::SVector{6})
    xc, yc = Tuple(first(first_guess) + last(first_guess)) ./ 2
    u, v, dudx, dudy, dvdx, dvdy = Tuple(X)
    sol = mappedarray(first_guess) do I
        x, y = Tuple(I)
        dx = x - xc
        dy = y - yc
        x′ = x + u + dudx*dx + dudy*dy
        y′ = y + v + dvdx*dx + dvdy*dy
        image_itp(x′, y′)
    end
    zncc(mappedarray(Float64, subset), sol)
end

# TODO: for 3D
# solution_vector
# compute_correlation

function fine_search(subset::AbstractArray{<: Gray, dim}, image::AbstractArray{<: Gray, dim}, first_guess::CartesianIndices{dim}) where {dim}
    @assert size(subset) == size(first_guess)
    image_itp = interpolate(mappedarray(Float64, image), BSpline(Linear())) # sub-pixel interpolation
    x = solution_vector(Val(dim))
    result = DiffResults.HessianResult(x)
    C = 0.0
    for i in 1:20
        result = ForwardDiff.hessian!(result, x -> compute_correlation(subset, image_itp, first_guess, x), x)
        C = DiffResults.value(result)
        C ≈ 1 && break
        ∇x = DiffResults.gradient(result)
        H = DiffResults.hessian(result)
        x += -H \ ∇x
    end
    center = Tuple(first(first_guess) + last(first_guess)) ./ 2
    ntuple(i -> center[i] + x[i], Val(dim)), C
end
