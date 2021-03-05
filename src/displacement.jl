struct DisplacementField{T <: SVector, N, ITP <: Interpolations.BSplineInterpolation{T, N}, V} <: AbstractArray{T, N}
    itp::ITP
    coords::Coordinate{N, Int, V}
    Cs::Array{Float64, N}
end

Base.size(x::DisplacementField) = size(x.itp)

@inline function Base.getindex(x::DisplacementField{<: Any, N}, I::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(x, I...)
    @inbounds x.itp[I...]
end

@inline function (x::DisplacementField{<: Any, N})(I::Vararg{Real, N}) where {N}
    @_propagate_inbounds_meta
    x.itp(I...)
end

@inline function strain(x::DisplacementField{<: Any, N}, I::Vararg{Real, N}) where {N}
    dudx = reduce(hcat, Interpolations.gradient(x, I...))
    (dudx + dudx') / 2
end

function displacement_field((before, after)::Pair{<: AbstractArray, <: AbstractArray}, sample_points::Coordinate{dim, Int}, npixels::Int, surrounding_npixels::Int = 5*npixels; thresh::Real = -Inf) where {dim}
    @assert size(before) == size(after)
    disp = similar(sample_points, SVector{dim, Float64})
    Cs = similar(sample_points, Float64)
    p = Progress(length(sample_points))
    Threads.@threads for i in eachindex(sample_points)
        indices = neighborindices(CartesianIndex(sample_points[i]), before, npixels)
        subset = before[indices]
        first_guess, C = coarse_search(subset, after; region = neighborindices(indices, after, surrounding_npixels), parallel = false)
        center, C = fine_search(subset, after, first_guess)
        disp[i] = center - SVector(Tuple(first(indices) + last(indices))) / 2
        Cs[i] = C
        next!(p)
    end
    DisplacementField(interpolate(disp, BSpline(Linear())), sample_points, Cs)
end
