@testset "Propulsion Types" begin

    # Heliocentric constants
    μ_sun = 1.32712440018e11  # km³/s²
    r_1AU = 1.495978707e8     # km

    @testset "SEP propagation" begin
        sep = SEPSpacecraft(200.0, 800.0, 0.5, 3000.0, r_1AU)

        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = mass(sep)

        throttle = SVector{3}(0.5, 0.0, 0.0)
        Δt = 86400.0

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

        @test mf < m0
        @test mf > 0
        @test rf != r0
        @test Δv > 0
    end

    @testset "SEP transfer" begin
        sep = SEPSpacecraft(200.0, 800.0, 0.5, 3000.0, r_1AU)

        r0 = [r_1AU, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [0.0, r_1AU * 1.1, 0.0]
        vf = [-27.0, 0.0, 0.0]
        tof = 86400.0 * 120

        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ_sun,
            sep;
            n_segments = 12,
            verbosity = 0,
            tol = 1e-4,
        )

        sol = solve(prob; max_iter = 500, initial_guess_strategy = RandomGuess(seed = 42))

        @test sol isa SimsFlanaganSolution
        @test sol.masses[end] <= mass(sep)
        @test sol.masses[end] > 0
    end

    @testset "Solar sail propagation" begin
        sail = SolarSail(100.0, 1000.0, 0.9)

        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = mass(sail)

        throttle = SVector{3}(1.0, 0.0, 0.0)
        Δt = 86400.0

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

        @test mf == m0  # Mass constant for solar sail
        @test rf != r0
        @test Δv > 0
    end

    @testset "Solar sail constant mass through leg" begin
        sail = SolarSail(100.0, 1000.0, 0.9)

        r0 = SVector{3}(r_1AU, 0.0, 0.0)
        v0 = SVector{3}(0.0, 29.78, 0.0)
        m0 = mass(sail)

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

        @test mf == m0
        @test total_Δv > 0
    end

    @testset "Solar sail problem" begin
        sail = SolarSail(100.0, 5000.0, 0.9)

        r0 = [r_1AU, 0.0, 0.0]
        v0 = [0.0, 29.78, 0.0]
        rf = [r_1AU * 1.02, 0.0, 0.0]
        vf = [0.0, 29.5, 0.0]
        tof = 86400.0 * 60

        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ_sun,
            sail;
            n_segments = 8,
            verbosity = 0,
            tol = 1e-4,
        )

        sol = solve(prob; initial_guess_strategy = RadialGuess())

        @test all(m -> m == mass(sail), sol.masses)
    end

    @testset "compute_segment_masses" begin
        sail = SolarSail(100.0, 1000.0, 0.9)
        prob_sail = simsflanagan_problem(
            [r_1AU, 0.0, 0.0],
            [0.0, 29.78, 0.0],
            [r_1AU * 1.05, 0.0, 0.0],
            [0.0, 29.0, 0.0],
            86400.0 * 30,
            μ_sun,
            sail;
            n_segments = 4,
            verbosity = 0,
        )

        throttles = [SVector{3}(0.5, 0.0, 0.0) for _ = 1:4]
        masses = SimsFlanagan.compute_segment_masses(prob_sail, throttles)
        @test all(m -> m == mass(sail), masses)

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob_sc = simsflanagan_problem(
            [r_1AU, 0.0, 0.0],
            [0.0, 29.78, 0.0],
            [r_1AU * 1.05, 0.0, 0.0],
            [0.0, 29.0, 0.0],
            86400.0 * 30,
            μ_sun,
            sc;
            n_segments = 4,
            verbosity = 0,
        )

        masses_sc = SimsFlanagan.compute_segment_masses(prob_sc, throttles)
        @test masses_sc[1] == mass(sc)
        @test all(i -> masses_sc[i+1] <= masses_sc[i], 1:4)
    end

    @testset "Spacecraft types produce valid solutions" begin
        # Simpler Earth orbit problem that converges reliably
        μ = 398600.4418  # Earth
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, sqrt(μ/7000.0), 0.0]
        rf = [0.0, 8000.0, 0.0]
        vf = [-sqrt(μ/8000.0), 0.0, 0.0]
        tof = 86400.0  # 1 day

        @testset "Low-Thrust" begin
            sc = Spacecraft(100.0, 400.0, 10.0, 3000.0)
            prob = simsflanagan_problem(
                r0,
                v0,
                rf,
                vf,
                tof,
                μ,
                sc;
                n_segments = 10,
                verbosity = 0,
                tol = 1e-4,
            )

            sol =
                solve(prob; max_iter = 500, initial_guess_strategy = RandomGuess(seed = 42))

            @test sol isa SimsFlanaganSolution
            @test sol.masses[end] <= mass(sc)
            @test sol.masses[end] > 0
            @test sol.Δv_total > 0
        end

        @testset "SEP" begin
            sep = SEPSpacecraft(100.0, 400.0, 10.0, 3000.0, 7000.0)
            prob = simsflanagan_problem(
                r0,
                v0,
                rf,
                vf,
                tof,
                μ,
                sep;
                n_segments = 10,
                verbosity = 0,
                tol = 1e-4,
            )

            sol =
                solve(prob; max_iter = 500, initial_guess_strategy = RandomGuess(seed = 42))

            @test sol isa SimsFlanaganSolution
            @test sol.masses[end] <= mass(sep)
            @test sol.masses[end] > 0
        end

        @testset "Solar Sail" begin
            # Solar sails are designed for heliocentric orbits, use sun-centered problem
            sail = SolarSail(100.0, 5000.0, 0.9)
            
            # Sun-centered initial and final conditions (similar to test at line 114)
            r0_sun = [r_1AU, 0.0, 0.0]
            v0_sun = [0.0, 29.78, 0.0]
            rf_sun = [r_1AU * 1.02, 0.0, 0.0]
            vf_sun = [0.0, 29.5, 0.0]
            tof_sun = 86400.0 * 370  # 60 days
            
            prob = simsflanagan_problem(
                r0_sun,
                v0_sun,
                rf_sun,
                vf_sun,
                tof_sun,
                μ_sun,
                sail;
                n_segments = 15,
                verbosity = 0,
                tol = 1e-4,
            )

            sol = solve(prob; initial_guess_strategy = RadialGuess())

            println(position_mismatch_norm(sol))
            println(velocity_mismatch_norm(sol))

            @test sol isa SimsFlanaganSolution
            @test sol.masses[1] == mass(sail)
            @test sol.masses[end] == mass(sail)
        end
    end
end
