@testset "Propulsion Types" begin

    @testset "SEP propagation" begin
        # Create SEP spacecraft
        r_1AU = 1.495978707e8  # km
        sep = SEPSpacecraft(1000.0, 0.5, 3000.0, r_1AU)

        # Sun-centric orbit starting at 1 AU
        μ_sun = 1.32712440018e11  # km³/s²
        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)  # ~Earth orbital velocity
        m0 = sep.mass

        # Simple throttle
        throttle = SVector{3}(0.5, 0.0, 0.0)
        Δt = 86400.0  # 1 day

        # Propagate segment
        rf, vf, mf, Δv = SimsFlanagan.propagate_segment(
            r0,
            v0,
            m0,
            throttle,
            Δt,
            μ_sun,
            sep;
            forward = Val(true),
        )

        # Check mass decreased
        @test mf < m0

        # Check position changed
        @test rf != r0

        # Check Δv is positive
        @test Δv > 0
    end

    @testset "SEP problem construction and solve" begin
        # Heliocentric transfer problem (non-collinear to avoid Lambert issues)
        r_1AU = 1.495978707e8  # km
        μ_sun = 1.32712440018e11  # km³/s²

        sep = SEPSpacecraft(1000.0, 0.5, 3000.0, r_1AU)

        # Non-collinear initial and final positions
        r0 = [r_1AU, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [0.0, r_1AU * 1.1, 0.0]  # 90 degree transfer
        vf = [-28.5, 0.0, 0.0]
        tof = 86400.0 * 90  # 90 days for quarter orbit

        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ_sun,
            sep;
            n_segments = 4,
            verbosity = 0,
        )

        @test prob isa SimsFlanaganProblem
        @test prob.spacecraft isa SEPSpacecraft

        # Solve with zero guess to avoid Lambert issues
        sol = simsflanagan_solve(prob; use_lambert_guess = Val(false))
        @test sol isa SimsFlanaganSolution
    end

    @testset "Solar sail propagation" begin
        # Create solar sail
        r_1AU = 1.495978707e8  # km
        sail = SolarSail(100.0, 1000.0, 0.9, r_1AU)

        # Sun-centric orbit starting at 1 AU
        μ_sun = 1.32712440018e11  # km³/s²
        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = sail.mass

        # Throttle pointing radially (sail perpendicular to sun)
        throttle = SVector{3}(1.0, 0.0, 0.0)
        Δt = 86400.0  # 1 day

        # Propagate segment
        rf, vf, mf, Δv = SimsFlanagan.propagate_segment(
            r0,
            v0,
            m0,
            throttle,
            Δt,
            μ_sun,
            sail;
            forward = Val(true),
        )

        # Check mass is CONSTANT for solar sail
        @test mf == m0

        # Check position changed
        @test rf != r0

        # Check Δv is positive (effective Δv)
        @test Δv > 0
    end

    @testset "Solar sail constant mass" begin
        r_1AU = 1.495978707e8
        sail = SolarSail(100.0, 1000.0, 0.9, r_1AU)
        μ_sun = 1.32712440018e11

        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = sail.mass

        # Propagate multiple segments
        throttles = [SVector{3}(1.0, 0.0, 0.0) for _ = 1:5]
        Δt_seg = 86400.0

        rf, vf, mf, total_Δv = SimsFlanagan.propagate_leg(
            r0,
            v0,
            m0,
            throttles,
            Δt_seg,
            μ_sun,
            sail;
            forward = Val(true),
        )

        # Mass should stay constant throughout
        @test mf == m0
    end

    @testset "Solar sail problem construction" begin
        r_1AU = 1.495978707e8
        μ_sun = 1.32712440018e11

        sail = SolarSail(100.0, 1000.0, 0.9, r_1AU)

        r0 = [r_1AU, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [r_1AU * 1.05, 0.0, 0.0]
        vf = [0.0, 29.0, 0.0]
        tof = 86400.0 * 60  # 60 days

        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ_sun,
            sail;
            n_segments = 4,
            verbosity = 0,
        )

        @test prob isa SimsFlanaganProblem
        @test prob.spacecraft isa SolarSail

        # Initial guess should fall back to zero for solar sails
        guess = SimsFlanagan.initial_guess_lambert(prob)
        @test all(g -> norm(g) == 0, guess)
    end

    @testset "compute_segment_masses for solar sail" begin
        r_1AU = 1.495978707e8
        μ_sun = 1.32712440018e11
        sail = SolarSail(100.0, 1000.0, 0.9, r_1AU)

        r0 = [r_1AU, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [r_1AU * 1.05, 0.0, 0.0]
        vf = [0.0, 29.0, 0.0]
        tof = 86400.0 * 30

        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ_sun,
            sail;
            n_segments = 4,
            verbosity = 0,
        )

        throttles = [SVector{3}(0.5, 0.0, 0.0) for _ = 1:4]
        masses = SimsFlanagan.compute_segment_masses(prob, throttles)

        # All masses should be equal for solar sail
        @test all(m -> m == sail.mass, masses)
    end

end

