@testset "Solve" begin

    @testset "Basic solve" begin
        μ = 398600.4418
        # Use non-collinear positions
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [6000.0, 4000.0, 500.0]
        vf = [0.0, 7.8, 0.1]
        tof = 86400.0 * 5

        sc = Spacecraft(1000.0, 0.5, 3000.0)
        prob =
            simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 6, verbosity = 0)

        sol = simsflanagan_solve(prob)

        @test sol isa SimsFlanaganSolution
        @test length(sol.throttles) == 6
        @test length(sol.masses) == 7  # n_segments + 1
        @test sol.Δv_total >= 0
        @test sol.mismatch isa SVector{7}
    end

    @testset "Solve with zero throttle guess" begin
        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7000.0, 1000.0, 0.0]
        vf = [0.0, 7.5, 0.0]
        tof = 86400.0 * 2

        sc = Spacecraft(1000.0, 0.5, 3000.0)
        prob =
            simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4, verbosity = 0)

        sol = simsflanagan_solve(prob; use_lambert_guess = Val(false))

        @test sol isa SimsFlanaganSolution
    end

    @testset "Solve with NelderMead (derivative-free)" begin
        using OptimizationOptimJL: Optim

        μ = 398600.4418
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7000.0, 2000.0, 0.0]
        vf = [0.0, 7.4, 0.1]
        tof = 86400.0 * 3

        sc = Spacecraft(1000.0, 0.5, 3000.0)
        prob =
            simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4, verbosity = 0)

        sol = simsflanagan_solve(prob; alg = Optim.NelderMead())

        @test sol isa SimsFlanaganSolution
        @test sol.Δv_total >= 0
    end

end

