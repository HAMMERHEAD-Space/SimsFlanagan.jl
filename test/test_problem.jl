@testset "Problem Construction" begin

    @testset "SimsFlanaganProblem construction" begin
        # Define a simple Earth orbit raising problem
        μ = 398600.4418  # km³/s²
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [42164.0, 0.0, 0.0]
        vf = [0.0, 3.075, 0.0]
        tof = 86400.0 * 30  # 30 days in seconds

        # Spacecraft(dry_mass, wet_mass, thrust, isp)
        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)

        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc)

        @test prob.r0 ≈ SVector{3}(r0)
        @test prob.v0 ≈ SVector{3}(v0)
        @test prob.rf ≈ SVector{3}(rf)
        @test prob.vf ≈ SVector{3}(vf)
        @test prob.tof == tof
        @test prob.μ == μ
        @test mass(prob.spacecraft) == mass(sc)  # Use mass() accessor
        @test prob.options.n_segments == 10
    end

    @testset "SciMLBase problem interface" begin
        using SciMLBase

        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4)

        # Test that problem inherits from AbstractSciMLProblem
        @test prob isa SciMLBase.AbstractSciMLProblem
        @test prob isa SimsFlanagan.AbstractSimsFlanaganProblem
    end

    @testset "SciMLBase remake" begin
        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4)

        # Test remake with single field
        prob_new_tof = remake(prob; tof = tof * 2)
        @test prob_new_tof.tof == tof * 2
        @test prob_new_tof.r0 == prob.r0
        @test prob_new_tof.v0 == prob.v0
        @test prob_new_tof.rf == prob.rf
        @test prob_new_tof.vf == prob.vf
        @test prob_new_tof.μ == prob.μ

        # Test remake with boundary conditions
        new_r0 = [8000.0, 0.0, 0.0]
        new_rf = [9000.0, 0.0, 0.0]
        prob_new_bc = remake(prob; r0 = new_r0, rf = new_rf)
        @test prob_new_bc.r0 == new_r0
        @test prob_new_bc.rf == new_rf
        @test prob_new_bc.tof == prob.tof

        # Test remake with new spacecraft
        sc2 = Spacecraft(100.0, 400.0, 1.0, 2000.0)
        prob_new_sc = remake(prob; spacecraft = sc2)
        @test prob_new_sc.spacecraft === sc2
        @test prob_new_sc.tof == prob.tof

        # Test remake with new options
        new_opts = SimsFlanaganOptions(n_segments = 8)
        prob_new_opts = remake(prob; options = new_opts)
        @test prob_new_opts.options.n_segments == 8
        @test prob_new_opts.tof == prob.tof

        # Test remake preserves type
        @test prob_new_tof isa SimsFlanaganProblem
        @test prob_new_bc isa SimsFlanaganProblem
    end

    @testset "Initial guess from Lambert" begin
        μ = 398600.4418
        # Use non-collinear positions (Earth to GEO with inclination change)
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [30000.0, 30000.0, 1000.0]  # Non-collinear final position
        vf = [0.0, 3.0, 0.1]
        tof = 86400.0 * 30

        # Spacecraft(dry_mass, wet_mass, thrust, isp)
        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 10)

        throttles = SimsFlanagan.initial_guess_lambert(prob)

        @test length(throttles) == 10
        @test all(t -> t isa SVector{3}, throttles)
    end

    @testset "Initial guess strategies" begin
        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0 * 5

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 6)

        # Test ZeroGuess
        throttles_zero = SimsFlanagan.generate_initial_guess(prob, ZeroGuess())
        @test length(throttles_zero) == 6
        @test all(t -> all(iszero, t), throttles_zero)

        # Test RandomGuess
        throttles_random = SimsFlanagan.generate_initial_guess(prob, RandomGuess())
        @test length(throttles_random) == 6
        @test all(t -> t isa SVector{3}, throttles_random)

        # Test ConstantGuess
        throttles_const = SimsFlanagan.generate_initial_guess(
            prob,
            ConstantGuess(direction = [1, 0, 0], magnitude = 0.5),
        )
        @test length(throttles_const) == 6
        @test all(t -> norm(t) ≈ 0.5, throttles_const)

        # Test RadialGuess
        throttles_radial =
            SimsFlanagan.generate_initial_guess(prob, RadialGuess(magnitude = 0.3))
        @test length(throttles_radial) == 6
    end
end

