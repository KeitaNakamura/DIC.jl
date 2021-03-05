to_float(A::AbstractArray{<: Gray}) = mappedarray(Float64, A)
to_float(A::AbstractArray{<: RGB}) = to_float(Gray.(A))

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

function zncc(A::AbstractArray, B::AbstractArray)
    zncc(to_float(A), to_float(B))
end

"""
    coarse_search(subset, image; region = CartesianIndices(image)) -> indices, C

Perform coarse search `subset` in `image` using DIC.
Return the `indices` which has the highest correlation with `subset`.
Use `image[indices]` to get the found part of image.
The searching `region` (entire image by default) can also be specified
by `CartesianIndices` to reduce computations.

See also [`neighborindices`](@ref).

# Examples

```jldoctest
julia> image = rand(10,10);

julia> subset = image[3:5, 2:3];

julia> coarse_search(subset, image)
(CartesianIndex{2}[CartesianIndex(3, 2) CartesianIndex(3, 3); CartesianIndex(4, 2) CartesianIndex(4, 3); CartesianIndex(5, 2) CartesianIndex(5, 3)], 1.0)
```
"""
function coarse_search(subset::AbstractArray, image::AbstractArray; region::PixelIndices = CartesianIndices(image), parallel::Bool = true)
    inds = walkindices(subset, image; region)
    A = map(x -> x.indices, inds)
    Cs = similar(inds, Float64)
    if parallel
        Threads.@threads for i in eachindex(inds, Cs)
            @inbounds Cs[i] = zncc(view(image, inds[i]), subset)
        end
    else
        for i in eachindex(inds, Cs)
            @inbounds Cs[i] = zncc(view(image, inds[i]), subset)
        end
    end
    I = argmax(Cs)
    inds[I], Cs[I]
end

# for 2D
solution_vector(::Type{T}, ::Val{2}) where {T} = zero(SVector{6, T})
function compute_correlation(subset::AbstractArray{<: Real, 2}, image_itp::AbstractArray{T, 2}, first_guess::PixelIndices{2}, X::SVector{6}) where {T}
    xc, yc = Tuple(first(first_guess) + last(first_guess)) ./ 2
    u, v, dudx, dudy, dvdx, dvdy = Tuple(X)
    sol = mappedarray(first_guess) do I
        x, y = Tuple(I)
        dx = x - xc
        dy = y - yc
        x′ = x + u + dudx*dx + dudy*dy
        y′ = y + v + dvdx*dx + dvdy*dy
        # If calculted coordinates are out side of image, just return zero.
        # This means that out side of image are filled with black color.
        checkbounds(Bool, image_itp, x′, y′) ? image_itp(x′, y′) : zero(T)
    end
    zncc(subset, sol)
end

# TODO: for 3D
# solution_vector
# compute_correlation

"""
    fine_search(subset, image, first_guess::PixelIndices) -> center, C

Perform fine search `subset` in `image` based on the Newton-Raphson method.
The results by [`coarse_search`](@ref) can be used as `first_guess`.
Note that returned `center` is a center coordinates (not integer any more) of searched subset in `image`.

# Examples

```jldoctest
julia> image = DIC.testimage("buffalo");

julia> subset = image[100:300, 300:500];

julia> center, C = fine_search(subset, image, CartesianIndices((101:301, 301:501)))
([200.00000782067005, 400.00001094427904], 0.9999999999438116)
```
"""
function fine_search(subset::AbstractArray{T, dim}, image::AbstractArray{T, dim}, first_guess::PixelIndices{dim}) where {T <: Real, dim}
    @assert size(subset) == size(first_guess)
    image_itp = interpolate(image, BSpline(Linear())) # sub-pixel interpolation
    x = solution_vector(T, Val(dim))
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
    SVector(ntuple(i -> center[i] + x[i], Val(dim))), C
end

function fine_search(subset::AbstractArray, image::AbstractArray, first_guess::PixelIndices)
    fine_search(to_float(subset), to_float(image), first_guess)
end
