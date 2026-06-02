#=
Utility functions for SimsFlanagan.jl
=#

export safe_norm

# =============================================================================
# AD-Safe Norm (Huber-style)
# =============================================================================

"""
    safe_norm(x, δ=1e-8)

Compute vector norm with well-defined gradients everywhere using Huber-style smoothing.

Standard `norm(x)` has gradient `x/norm(x)` which is NaN at x=0.
This version uses: `sqrt(|x|² + δ²) - δ`

Properties:
- safe_norm(0) = 0 (exact)
- safe_norm(x) ≈ norm(x) for |x| >> δ  
- Gradient at x=0 is [0, 0, 0] (well-defined)
- Gradient is continuous everywhere
- C∞ smooth (infinitely differentiable)
"""
@inline function safe_norm(x::AbstractVector{T}, δ::T=T(1e-8)) where {T}
    sq = dot(x, x)
    return sqrt(sq + δ^2) - δ
end

# Euclidean (2-) norm computed with an explicit reduction. This intentionally
# avoids `LinearAlgebra.norm`, whose `StridedVector{<:Union{...}}` methods route
# through `ReinterpretArray` padding code that triggers a JET false positive on
# Julia 1.12 when analyzed against generic `AbstractVector{<:Number}` signatures.
@inline function _euclidean_norm(v::AbstractVector{T}) where {T<:Number}
    s = abs2(zero(T))
    @inbounds for i in eachindex(v)
        s += abs2(v[i])
    end
    return sqrt(s)
end

