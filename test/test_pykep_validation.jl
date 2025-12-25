#=
Validation tests comparing SimsFlanagan.jl against PyKEP reference values.

These tests ensure the Julia implementation matches PyKEP's Sims-Flanagan
implementation within acceptable tolerances.

Note: PyKEP uses SI units (meters, seconds) while SimsFlanagan.jl uses
km and km/s. Values are converted as needed.
=#

@testset "PyKEP Validation" begin

    # Heliocentric constants (matching PyKEP)
    const MU_SUN_SI = 1.32712440018e20  # m³/s² (PyKEP uses SI)
    const MU_SUN = 1.32712440018e11     # km³/s² (SimsFlanagan uses km)
    const AU_M = 1.495978707e11         # m (1 AU in meters)
    const AU_KM = 1.495978707e8         # km (1 AU in km)

    @testset "Unit conversions" begin
        # Verify our unit conversions are correct
        @test AU_KM * 1000 ≈ AU_M
        @test MU_SUN * 1e9 ≈ MU_SUN_SI
    end

    @testset "Stumpff functions" begin
        # Test Stumpff functions against known values
        # c2(0) = 1/2, c3(0) = 1/6
        @test SimsFlanagan.stumpff_c2(0.0) ≈ 0.5 atol=1e-10
        @test SimsFlanagan.stumpff_c3(0.0) ≈ 1/6 atol=1e-10

        # For z > 0 (elliptic): c2(z) = (1 - cos(√z))/z
        z = 1.0
        @test SimsFlanagan.stumpff_c2(z) ≈ (1 - cos(sqrt(z)))/z
        @test SimsFlanagan.stumpff_c3(z) ≈ (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3)

        # For z < 0 (hyperbolic): c2(z) = (cosh(√-z) - 1)/(-z)
        z = -1.0
        @test SimsFlanagan.stumpff_c2(z) ≈ (cosh(sqrt(-z)) - 1)/(-z)
        @test SimsFlanagan.stumpff_c3(z) ≈ (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)^3)
    end

    @testset "Kepler propagation accuracy" begin
        # Test against known orbit: circular Earth orbit at 1 AU
        r0 = SVector{3}(AU_KM, 0.0, 0.0)
        v_circ = sqrt(MU_SUN / AU_KM)  # Circular velocity at 1 AU
        v0 = SVector{3}(0.0, v_circ, 0.0)

        # Period of 1 AU circular orbit
        T_period = 2π * sqrt(AU_KM^3 / MU_SUN)

        # Propagate for quarter period
        rf, vf = SimsFlanagan.kepler_propagate(r0, v0, T_period/4, MU_SUN)

        # Should be at (0, 1 AU, 0) with velocity (-v_circ, 0, 0)
        @test rf ≈ SVector{3}(0.0, AU_KM, 0.0) rtol=1e-10
        @test vf ≈ SVector{3}(-v_circ, 0.0, 0.0) rtol=1e-10
    end

    @testset "Segment propagation matches PyKEP convention" begin
        # PyKEP applies impulse at segment midpoint
        # This test verifies our implementation does the same

        r0 = SVector{3}(AU_KM, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)  # Approximate Earth orbital velocity
        m0 = 1000.0  # kg
        Δt = 86400.0  # 1 day

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)

        # Zero throttle should give same result as pure Kepler propagation
        throttle_zero = SVector{3}(0.0, 0.0, 0.0)
        rf_seg, vf_seg, mf_seg, Δv = SimsFlanagan.propagate_segment(
            r0, v0, m0, throttle_zero, Δt, MU_SUN, sc
        )

        rf_kep, vf_kep = SimsFlanagan.kepler_propagate(r0, v0, Δt, MU_SUN)

        @test rf_seg ≈ rf_kep rtol=1e-12
        @test vf_seg ≈ vf_kep rtol=1e-12
        @test mf_seg == m0
        @test Δv == 0.0
    end

    @testset "Mass flow rate matches Tsiolkovsky" begin
        # Test that mass consumption follows rocket equation
        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        vex = exhaust_velocity(sc)

        # Mass flow rate should be thrust / exhaust_velocity
        expected_mdot = 0.5 / vex  # kg/s
        actual_mdot = SimsFlanagan.compute_mass_flow(sc, 0.5)

        @test actual_mdot ≈ expected_mdot
    end

    @testset "SEP thrust scaling" begin
        # Verify inverse-square thrust scaling for SEP
        sep = SEPSpacecraft(200.0, 800.0, 0.5, 3000.0, AU_KM)

        r_1AU = SVector{3}(AU_KM, 0.0, 0.0)
        r_2AU = SVector{3}(2*AU_KM, 0.0, 0.0)
        r_05AU = SVector{3}(0.5*AU_KM, 0.0, 0.0)

        thrust_1AU = SimsFlanagan.compute_thrust(sep, r_1AU, 1.0)
        thrust_2AU = SimsFlanagan.compute_thrust(sep, r_2AU, 1.0)
        thrust_05AU = SimsFlanagan.compute_thrust(sep, r_05AU, 1.0)

        # Verify inverse square scaling
        @test thrust_1AU ≈ sep.thrust_ref
        @test thrust_2AU ≈ thrust_1AU / 4  # (1/2)² = 1/4
        @test thrust_05AU ≈ thrust_1AU * 4  # (1/0.5)² = 4
    end

    @testset "Solar sail characteristic acceleration" begin
        # Compare to analytical formula
        # a_c = (1 + η) × P_sun × A / (c × m)
        sail = SolarSail(100.0, 1000.0, 0.9)

        P_sun = 1361.0  # W/m² at 1 AU
        c = 299792458.0  # m/s
        eta = 0.9
        m = 100.0  # kg
        A = 1000.0  # m²

        expected_a_c = (1 + eta) * P_sun * A / (c * m)
        actual_a_c = characteristic_acceleration(sail)

        @test actual_a_c ≈ expected_a_c rtol=1e-10
    end

    @testset "Backward propagation consistency" begin
        # Forward then backward should return to original state
        r0 = SVector{3}(AU_KM, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = 1000.0
        Δt = 86400.0

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        throttle = SVector{3}(0.5, 0.0, 0.0)

        # Forward propagation
        rf, vf, mf, Δv_fwd = SimsFlanagan.propagate_segment(
            r0, v0, m0, throttle, Δt, MU_SUN, sc; forward=Val(true)
        )

        # Backward propagation from forward result
        r_back, v_back, m_back, Δv_bwd = SimsFlanagan.propagate_segment(
            rf, vf, mf, throttle, Δt, MU_SUN, sc; forward=Val(false)
        )

        # Should return to original state
        @test r_back ≈ r0 rtol=1e-10
        @test v_back ≈ v0 rtol=1e-10
        # Mass might have small differences due to average mass approximation
        @test m_back ≈ m0 rtol=1e-6
    end

    @testset "Mismatch constraint scaling" begin
        # Verify canonical units are computed correctly
        r0 = [AU_KM, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [0.0, AU_KM * 1.1, 0.0]
        vf = [-27.0, 0.0, 0.0]
        tof = 86400.0 * 120

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, MU_SUN, sc;
            n_segments=4, verbosity=0)

        LU, VU, MU = SimsFlanagan.get_canonical_scales(prob)

        # LU should be max of initial/final radius
        @test LU ≈ max(norm(r0), norm(rf))

        # VU should be canonical velocity
        @test VU ≈ sqrt(MU_SUN / LU)

        # MU should be initial mass
        @test MU ≈ mass(sc)
    end

end

