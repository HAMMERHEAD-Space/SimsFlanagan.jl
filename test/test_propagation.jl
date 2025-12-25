@testset "Propagation" begin

    @testset "Kepler propagation" begin
        # Simple circular orbit test
        μ = 398600.4418
        r0 = SVector{3}(7000.0, 0.0, 0.0)
        v0 = SVector{3}(0.0, sqrt(μ/7000.0), 0.0)  # Circular orbit velocity

        # Period of circular orbit
        T = 2π * sqrt(7000.0^3 / μ)

        # Propagate for one full period
        rf, vf = SimsFlanagan.kepler_propagate(r0, v0, T, μ)

        # Should return to initial state
        @test rf ≈ r0 rtol=1e-11  # Within 1e-11 km
        @test vf ≈ v0 rtol=1e-15  # Within 1e-15 km/s

        # Propagate for half period
        rf_half, vf_half = SimsFlanagan.kepler_propagate(r0, v0, T/2, μ)

        # Position should be on opposite side
        @test rf_half ≈ -r0 rtol=1e-15  # Within 1e-15 km
        @test vf_half ≈ -v0 rtol=1e-15  # Within 1e-15 km/s
    end

    @testset "Segment propagation" begin
        μ = 398600.4418
        r = SVector{3}(7000.0, 0.0, 0.0)
        v = SVector{3}(0.0, 7.546, 0.0)
        m = 1000.0
        Δt = 86400.0  # 1 day

        # Spacecraft(dry_mass, wet_mass, thrust, isp)
        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)

        # Zero throttle should just coast
        throttle_zero = SVector{3}(0.0, 0.0, 0.0)
        rf, vf, mf, Δv = SimsFlanagan.propagate_segment(r, v, m, throttle_zero, Δt, μ, sc)

        @test mf == m  # No mass change
        @test Δv == 0.0  # No ΔV

        # Non-zero throttle should consume mass and add ΔV
        throttle = SVector{3}(1.0, 0.0, 0.0)
        rf2, vf2, mf2, Δv2 = SimsFlanagan.propagate_segment(r, v, m, throttle, Δt, μ, sc)

        @test mf2 < m  # Mass decreased
        @test Δv2 > 0  # ΔV applied
    end

    @testset "Mismatch computation" begin
        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]  # Slight orbit change
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0 * 5

        # Spacecraft(dry_mass, wet_mass, thrust, isp)
        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4)

        # Zero throttle guess
        throttles = [SVector{3}(0.0, 0.0, 0.0) for _ = 1:4]

        mismatch = compute_mismatch(prob, throttles)

        @test length(mismatch) == 7  # 3 position + 3 velocity + 1 mass
        @test mismatch isa SVector{7}
    end

end

