@testset "Solve" begin
    μ = 398600.4418  # Earth gravitational parameter

    @testset "LEO orbit raising" begin
        r1 = 7000.0
        r2 = 7500.0
        v1 = sqrt(μ / r1)
        v2 = sqrt(μ / r2)

        r0 = [r1, 0.0, 0.0]
        v0 = [0.0, v1, 0.0]
        rf = [r2, 0.0, 0.0]
        vf = [0.0, v2, 0.0]

        tof = 86400.0 * 2

        sc = Spacecraft(500.0, 500.0, 0.5, 3000.0)
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
            tol = 1e-6,
        )

        sol = solve(prob; initial_guess_strategy = RandomGuess(seed = 42))

        @test sol isa SimsFlanaganSolution
        @test length(sol.throttles) == 10

        @test position_mismatch(sol) isa SVector{3}
        @test velocity_mismatch(sol) isa SVector{3}
        @test mass_mismatch(sol) isa Number

        @test sol.converged == (sol.retcode == ReturnCode.Success)
    end

    @testset "Multiple segment counts" begin
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0 * 2

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)

        for n_seg in [4, 8]
            prob = simsflanagan_problem(
                r0,
                v0,
                rf,
                vf,
                tof,
                μ,
                sc;
                n_segments = n_seg,
                verbosity = 0,
            )

            sol = solve(prob)
            @test length(sol.throttles) == n_seg
        end
    end

    @testset "Infeasible problem" begin
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [42000.0, 0.0, 0.0]
        vf = [0.0, 3.07, 0.0]
        tof = 86400.0 * 0.1

        sc = Spacecraft(990.0, 10.0, 0.001, 3000.0)
        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ,
            sc;
            n_segments = 4,
            verbosity = 0,
            tol = 1e-6,
        )

        sol = solve(prob; max_iter = 100)

        @test !sol.converged
        @test sol.retcode != ReturnCode.Success
    end

    @testset "SciMLBase interface" begin
        r0 = [7000.0, 0.0, 0.0]
        v0 = [0.0, 7.546, 0.0]
        rf = [7500.0, 0.0, 0.0]
        vf = [0.0, 7.3, 0.0]
        tof = 86400.0

        sc = Spacecraft(200.0, 800.0, 0.5, 3000.0)
        prob =
            simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments = 4, verbosity = 0)

        sol1 = solve(prob)
        @test sol1 isa SimsFlanaganSolution

        sol2 = solve(prob; max_iter = 50, tol = 1e-4)
        @test sol2 isa SimsFlanaganSolution

        prob2 = remake(prob; tof = prob.tof * 1.5)
        @test prob2.tof == prob.tof * 1.5
        @test prob2.r0 == prob.r0

        sol_random = solve(prob; initial_guess_strategy = RandomGuess())
        @test sol_random isa SimsFlanaganSolution

        sol_zero = solve(prob; initial_guess_strategy = ZeroGuess())
        @test sol_zero isa SimsFlanaganSolution
    end
end
