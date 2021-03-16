mutable struct Camera{T}
    # internal parameters
    f_δx::MVector{2, T} # (focal length) * (pixel per distance)
    x₀::MVector{2, T} # offsets in image
    # external parameters
    t::MVector{3, T} # translation
    R::MMatrix{3, 3, T, 9} # rotation
    # lhs * P_vector = x
    P::MMatrix{3, 4, T} # projection matrix
    lhs::Vector{T} # lhs matrix (reshaped later)
    rhs::Vector{T} # rhs vector
end

function Camera{T}() where {T}
    f = zero(T)
    δx = x₀ = zero(MVector{2, T})
    t = zero(MVector{3, T})
    R = zero(MMatrix{3, 3, T, 9})
    P = zero(MMatrix{11, T})
    Camera{T}(f, δx, x₀, t, R, P)
end

Camera() = Camera{Float64}()

function Base.push!(camera::Camera, (x, X)::Pair{<: AbstractVector, <: AbstractVector})
    @assert length(x) == 2
    @assert length(X) == 3
    𝟘 = zeros(length(X) + 1)
    append!(camera.lhs, [X; 1;    𝟘; -x[1].*X])
    append!(camera.lhs, [𝟘;    X; 1; -x[2].*X])
    append!(camera.rhs, x)
    camera
end

function calibrate!(camera::Camera)
    length(camera.rhs) > 11 || @warn "number of sample points should be at least 6, but $(length(camera.rhs) ÷ 2)"

    # solve projection matrix
    sol = camera.lhs \ camera.rhs
    push!(sol, 1) # fill element at (3,4)
    camera.P .= reshape(sol, 4, 3)'

    R, A = qr(camera.P[1:3, 1:3])
    # internal parameters
    camera.f_δx .= [A[1,1], A[2,2]]
    camera.x₀ .= A[1:2, 3]
    # external parameters
    camera.R .= R
    camera.t .= inv(A) * camera.P[1:3, 4]

    camera
end
