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
@inline function safe_norm(x::AbstractVector{T}, δ::T = T(1e-8)) where {T}
    sq = dot(x, x)
    return sqrt(sq + δ^2) - δ
end

