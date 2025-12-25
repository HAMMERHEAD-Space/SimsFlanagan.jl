#=
Test Sundman transformation for adaptive mesh in trajectory optimization.
=#

using Test
using SimsFlanagan
using LinearAlgebra

# Constants
const MU_SUN = 1.32712440018e11  # km³/s²
const AU = 1.495978707e8  # km

@testset "Sundman Transformation" begin

    @testset "Equal segments when c=0" begin
        r0 = [AU, 0.0, 0.0]
        rf = [1.5*AU, 0.0, 0.0]
        tof = 200.0 * 86400.0

        Δt_fwd, Δt_bwd = SimsFlanagan.compute_sundman_leg_times(r0, rf, tof, 5, 5; c = 0.0)

        @test length(Δt_fwd) == 5
        @test length(Δt_bwd) == 5
        @test all(Δt_fwd .≈ tof / 10)
        @test all(Δt_bwd .≈ tof / 10)
    end

    @testset "Sundman c=1 gives varying segments" begin
        r0 = [AU, 0.0, 0.0]
        rf = [1.5*AU, 0.0, 0.0]
        tof = 200.0 * 86400.0

        Δt_fwd, Δt_bwd = SimsFlanagan.compute_sundman_leg_times(r0, rf, tof, 5, 5; c = 1.0)

        # Moving outward: later segments should be longer
        @test Δt_fwd[5] > Δt_fwd[1]
        @test sum(Δt_fwd) + sum(Δt_bwd) ≈ tof
    end

    @testset "Sundman c=1.5 for eccentric anomaly" begin
        r0 = [0.7*AU, 0.0, 0.0]
        rf = [1.5*AU, 0.0, 0.0]
        tof = 180.0 * 86400.0

        Δt_fwd, Δt_bwd = SimsFlanagan.compute_sundman_leg_times(r0, rf, tof, 5, 5; c = 1.5)

        # Larger c means more variation
        ratio = Δt_bwd[end] / Δt_fwd[1]
        @test ratio > 2.0
        @test sum(Δt_fwd) + sum(Δt_bwd) ≈ tof
    end
end

@testset "Transfer with Sundman" begin
    # Earth orbit transfer (7000 -> 8000 km)
    μ = 398600.4418
    r1, r2 = 7000.0, 8000.0

    r0 = [r1, 0.0, 0.0]
    v0 = [0.0, sqrt(μ/r1), 0.0]
    rf = [0.0, r2, 0.0]
    vf = [-sqrt(μ/r2), 0.0, 0.0]
    tof = 86400.0  # 1 day

    sc = Spacecraft(100.0, 400.0, 10.0, 3000.0)

    @testset "c=0 (equal segments)" begin
        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ,
            sc;
            n_segments = 10,
            n_fwd = 5,
            verbosity = 0,
            sundman_c = 0.0,
        )

        sol = solve(
            prob;
            initial_guess_strategy = RandomGuess(seed = 42),
            max_iter = 500,
            tol = 1e-4,
        )

        @test position_mismatch_norm(sol) < 10.0
    end

    @testset "c=1.0 (Sundman)" begin
        prob = simsflanagan_problem(
            r0,
            v0,
            rf,
            vf,
            tof,
            μ,
            sc;
            n_segments = 10,
            n_fwd = 5,
            verbosity = 0,
            sundman_c = 1.0,
        )

        sol = solve(
            prob;
            initial_guess_strategy = RandomGuess(seed = 42),
            max_iter = 500,
            tol = 1e-4,
        )

        @test position_mismatch_norm(sol) < 10.0
    end
end
