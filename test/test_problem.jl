@testset "Problem Construction" begin

    @testset "SimsFlanaganProblem construction" begin
        # Define a simple Earth orbit raising problem
        μ = 398600.4418  # km³/s²
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [42164.0, 0.0, 0.0]
        vf = [0.0, 3.075, 0.0]
        tof = 86400.0 * 30  # 30 days in seconds

        sc = Spacecraft(1000.0, 0.5, 3000.0)

        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc)

        @test prob.r0 ≈ SVector{3}(r0)
        @test prob.v0 ≈ SVector{3}(v0)
        @test prob.rf ≈ SVector{3}(rf)
        @test prob.vf ≈ SVector{3}(vf)
        @test prob.tof == tof
        @test prob.μ == μ
        @test prob.spacecraft.mass == sc.mass
        @test prob.options.n_segments == 10
    end

    @testset "Initial guess from Lambert" begin
        μ = 398600.4418
        # Use non-collinear positions (Earth to GEO with inclination change)
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [30000.0, 30000.0, 1000.0]  # Non-collinear final position
        vf = [0.0, 3.0, 0.1]
        tof = 86400.0 * 30

        sc = Spacecraft(1000.0, 0.5, 3000.0)
        prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 10)

        throttles = SimsFlanagan.initial_guess_lambert(prob)

        @test length(throttles) == 10
        @test all(t -> t isa SVector{3}, throttles)
    end
end

