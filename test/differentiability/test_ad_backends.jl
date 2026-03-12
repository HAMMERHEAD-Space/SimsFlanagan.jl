# AD backend tests for SimsFlanagan.jl
#
# Validates that each AD backend produces Jacobians matching ForwardDiff
# (the reference) when differentiating through core SimsFlanagan functions.
#
# Expects:
#   _SF_BACKENDS — Tuple of (name::String, backend) pairs (defined in runtests.jl)

using DifferentiationInterface: jacobian, AutoForwardDiff

const _SF_REF_BACKEND = AutoForwardDiff()

# Test problem parameters
const _SF_MU = 398600.4418  # km³/s² (Earth)
const _SF_R0 = SVector{3}(7000.0, 0.0, 0.0)
const _SF_V0 = SVector{3}(0.0, 7.546, 0.0)
const _SF_RF = SVector{3}(30000.0, 30000.0, 1000.0)
const _SF_VF = SVector{3}(0.0, 3.0, 0.1)
const _SF_TOF = 86400.0 * 5  # 5 days
const _SF_SPACECRAFT = Spacecraft(200.0, 800.0, 0.5, 3000.0)

# Create test problem
const _SF_TEST_PROB = simsflanagan_problem(
    _SF_R0, _SF_V0, _SF_RF, _SF_VF, _SF_TOF, _SF_MU, _SF_SPACECRAFT;
    n_segments = 4
)

# Test throttle vectors
const _SF_THROTTLES = [
    SVector{3}(0.1, 0.2, 0.3),
    SVector{3}(0.2, 0.1, 0.0),
    SVector{3}(0.0, 0.1, 0.2),
    SVector{3}(0.3, 0.2, 0.1),
]

# Initial throttle guess (flattened)
const _SF_THROTTLE_X0 = Vector(vcat([0.1, 0.2, 0.3, 0.2, 0.1, 0.0, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1]...))

# Segment times for compute_total_Δv
const _SF_DT_SEGMENTS = fill(_SF_TOF / 4, 4)

"""
Map input vector x = [throttle1..., throttle2..., ..., throttleN...] to mismatch output.
"""
function _mismatch_map(x::AbstractVector{T}, prob) where {T}
    n_seg = prob.options.n_segments
    throttles = [SVector{3,T}(x[3*(i-1)+1], x[3*(i-1)+2], x[3*(i-1)+3]) for i in 1:n_seg]
    mismatch = compute_mismatch(prob, throttles)
    return Vector(mismatch)
end

"""
Map input vector x = [throttle1..., throttle2..., ..., throttleN...] to scaled mismatch output.
"""
function _scaled_mismatch_map(x::AbstractVector{T}, prob) where {T}
    n_seg = prob.options.n_segments
    throttles = [SVector{3,T}(x[3*(i-1)+1], x[3*(i-1)+2], x[3*(i-1)+3]) for i in 1:n_seg]
    mismatch = scaled_mismatch_constraints(prob, throttles)
    return Vector(mismatch)
end

"""
Map input vector x = [throttle1..., throttle2..., ..., throttleN...] to total ΔV.
"""
function _total_dv_map(x::AbstractVector{T}, prob, Δt_segments) where {T}
    n_seg = prob.options.n_segments
    throttles = [SVector{3,T}(x[3*(i-1)+1], x[3*(i-1)+2], x[3*(i-1)+3]) for i in 1:n_seg]
    Δv = compute_total_Δv(throttles, Δt_segments, prob.spacecraft)
    return [Δv]
end

"""
Map input vector x = [throttle1..., throttle2..., ..., throttleN...] to throttle magnitude constraints.
"""
function _throttle_mag_map(x::AbstractVector{T}) where {T}
    n_seg = length(x) ÷ 3
    throttles = [SVector{3,T}(x[3*(i-1)+1], x[3*(i-1)+2], x[3*(i-1)+3]) for i in 1:n_seg]
    return SimsFlanagan.throttle_magnitude_constraints(throttles)
end

"""
Map input state x = [r0..., v0..., Δt] to propagated state [rf..., vf...].
"""
function _kepler_map(x::AbstractVector{T}) where {T}
    r0 = SVector{3,T}(x[1], x[2], x[3])
    v0 = SVector{3,T}(x[4], x[5], x[6])
    Δt = x[7]
    rf, vf = SimsFlanagan.kepler_propagate(r0, v0, Δt, T(_SF_MU))
    return Vector(vcat(rf, vf))
end

"""
Map input x = [r..., v..., m, throttle..., Δt] to propagated segment output.
"""
function _segment_map(x::AbstractVector{T}) where {T}
    r = SVector{3,T}(x[1], x[2], x[3])
    v = SVector{3,T}(x[4], x[5], x[6])
    m = x[7]
    throttle = SVector{3,T}(x[8], x[9], x[10])
    Δt = x[11]
    
    rf, vf, mf, Δv = SimsFlanagan.propagate_segment(
        r, v, m, throttle, Δt, T(_SF_MU), _SF_SPACECRAFT
    )
    return Vector([rf..., vf..., mf, Δv])
end

# Test inputs
const _SF_KEPLER_X0 = Vector(vcat(_SF_R0, _SF_V0, 3600.0))  # 1 hour propagation

const _SF_SEGMENT_X0 = Vector([
    7000.0, 0.0, 0.0,      # r
    0.0, 7.546, 0.0,       # v
    1000.0,                # mass
    0.5, 0.3, 0.1,         # throttle
    86400.0                # Δt (1 day)
])

@testset "Kepler Propagation" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(_kepler_map, _SF_REF_BACKEND, _SF_KEPLER_X0)
            J_ad = jacobian(_kepler_map, backend, _SF_KEPLER_X0)
            @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-12
        end
    end
end

@testset "Segment Propagation" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(_segment_map, _SF_REF_BACKEND, _SF_SEGMENT_X0)
            J_ad = jacobian(_segment_map, backend, _SF_SEGMENT_X0)
            @test J_ad ≈ J_ref rtol = 1e-5 atol = 1e-10
        end
    end
end

@testset "Mismatch Computation" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(x -> _mismatch_map(x, _SF_TEST_PROB), _SF_REF_BACKEND, _SF_THROTTLE_X0)
            J_ad = jacobian(x -> _mismatch_map(x, _SF_TEST_PROB), backend, _SF_THROTTLE_X0)
            @test J_ad ≈ J_ref rtol = 1e-4 atol = 1e-8
        end
    end
end

@testset "Scaled Mismatch Constraints (solve objective)" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(x -> _scaled_mismatch_map(x, _SF_TEST_PROB), _SF_REF_BACKEND, _SF_THROTTLE_X0)
            J_ad = jacobian(x -> _scaled_mismatch_map(x, _SF_TEST_PROB), backend, _SF_THROTTLE_X0)
            @test J_ad ≈ J_ref rtol = 1e-4 atol = 1e-8
        end
    end
end

@testset "Total ΔV Computation (solve objective)" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(x -> _total_dv_map(x, _SF_TEST_PROB, _SF_DT_SEGMENTS), _SF_REF_BACKEND, _SF_THROTTLE_X0)
            J_ad = jacobian(x -> _total_dv_map(x, _SF_TEST_PROB, _SF_DT_SEGMENTS), backend, _SF_THROTTLE_X0)
            @test J_ad ≈ J_ref rtol = 1e-5 atol = 1e-10
        end
    end
end

@testset "Throttle Magnitude Constraints (solve inequality)" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            J_ref = jacobian(_throttle_mag_map, _SF_REF_BACKEND, _SF_THROTTLE_X0)
            J_ad = jacobian(_throttle_mag_map, backend, _SF_THROTTLE_X0)
            @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-14
        end
    end
end

@testset "Stumpff Functions" begin
    # Test Stumpff c2 differentiability
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname c2" begin
            for z in [-2.0, -0.5, -1e-8, 0.0, 1e-8, 0.5, 2.0]
                J_ref = jacobian(x -> [SimsFlanagan.stumpff_c2(x[1])], _SF_REF_BACKEND, [z])
                J_ad = jacobian(x -> [SimsFlanagan.stumpff_c2(x[1])], backend, [z])
                @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-14
            end
        end
    end
    
    # Test Stumpff c3 differentiability
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname c3" begin
            for z in [-2.0, -0.5, -1e-8, 0.0, 1e-8, 0.5, 2.0]
                J_ref = jacobian(x -> [SimsFlanagan.stumpff_c3(x[1])], _SF_REF_BACKEND, [z])
                J_ad = jacobian(x -> [SimsFlanagan.stumpff_c3(x[1])], backend, [z])
                @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-14
            end
        end
    end
end

@testset "safe_norm" begin
    for (bname, backend) in _SF_BACKENDS
        @testset "$bname" begin
            # Test at various points including near-zero
            for x in [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [0.001, 0.001, 0.001],
                [1e-10, 1e-10, 1e-10],
                [0.0, 0.0, 0.0],  # Critical: gradient at zero
            ]
                J_ref = jacobian(v -> [safe_norm(SVector{3}(v...))], _SF_REF_BACKEND, x)
                J_ad = jacobian(v -> [safe_norm(SVector{3}(v...))], backend, x)
                @test J_ad ≈ J_ref rtol = 1e-6 atol = 1e-14
                # Verify no NaN at zero
                @test all(isfinite, J_ad)
            end
        end
    end
end
