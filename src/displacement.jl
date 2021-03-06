struct DisplacementField{T <: SVector, N, CT, ITP <: Interpolations.BSplineInterpolation{T, N}, V} <: AbstractArray{T, N}
    before::Array{CT}
    after::Array{CT}
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
    dudx = reduce(hcat, Interpolations.gradient(x.itp, I...))
    (dudx + dudx') / 2
end

function plot(disp::DisplacementField{<: Any, 2}; kwargs...)
    # NOTE: need to swich x and y since the coordinates of pixels is differen from general coordinate system

    m, n = size(disp)
    l = min(m, n)

    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98), resolution = round.(Int, (m, n) ./ l .* 1400))
    image = reverse(disp.before', dims = 1)

    # plot arrows
    ax_image = fig[1,1] = Axis(fig)
    image!(ax_image, image)

    x = coordinateaxes(disp.coords, 2)
    y = coordinateaxes(disp.coords, 1)
    u = [x[2] for x in disp']
    v = [x[1] for x in disp']
    arrows!(ax_image, x, y, u, v; kwargs...)
    ax_image.yreversed = true

    # plot correlation contour
    ax_cor, hm = contourf(fig[1,2], x, y, disp.Cs'; colormap = :jet)
    arrows!(ax_cor, x, y, u, v; kwargs...)
    ax_cor.yreversed = true

    linkaxes!(ax_image, ax_cor)
    ax_image.aspect = DataAspect()
    ax_cor.aspect = DataAspect()

    Colorbar(fig[1,3], hm, width = 30)

    # plot volumetric strains
    ϵv = mappedarray(i -> tr(strain(disp, Tuple(i)...)), CartesianIndices(disp))'
    ax_vol, hm = contourf(fig[2,1], x, y, ϵv; colormap = :jet)
    ax_vol.yreversed = true
    linkaxes!(ax_image, ax_vol)
    ax_vol.aspect = DataAspect()

    Colorbar(fig[3,1], hm, height = 30, vertical = false)

    # plot deviatoric strains
    ϵd = mappedarray(CartesianIndices(disp)) do i
        ϵ = strain(disp, Tuple(i)...)
        e = ϵ - tr(ϵ) * one(ϵ)
        sqrt(dot(ϵ, ϵ))
    end |> transpose
    ax_dev, hm = contourf(fig[2,2], x, y, ϵd; colormap = :jet)
    ax_dev.yreversed = true
    linkaxes!(ax_image, ax_dev)
    ax_dev.aspect = DataAspect()

    Colorbar(fig[3,2], hm, height = 30, vertical = false)

    fig
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
    DisplacementField(before, after, interpolate(disp, BSpline(Linear())), sample_points, Cs)
end
