using AllocCheck

# AllocCheck reports spurious `jl_get_pgcstack_static` "allocating runtime
# call"s on macOS aarch64 with Julia 1.12+. These are not real heap allocations:
# the analyzed code is allocation-free on every other platform/version (Linux,
# Windows, and macOS on Julia 1.10/1.11). This is a known AllocCheck/Julia
# limitation, so the checks are skipped on the affected platform.
# Ref: https://github.com/SciML/SciMLStructures.jl/issues/59
const _SKIP_ALLOCCHECK = Sys.isapple() && Sys.ARCH === :aarch64 && VERSION >= v"1.12"

if _SKIP_ALLOCCHECK
    @info "Skipping AllocCheck allocation tests (spurious jl_get_pgcstack_static reports on macOS aarch64 + Julia 1.12+; see SciML/SciMLStructures.jl#59)."
end

function checked_allocs(f, types)
    _SKIP_ALLOCCHECK && return ()
    allocs = check_allocs(f, types)
    if !isempty(allocs)
        printstyled(stdout, "\n[ALLOC] "; color=:red, bold=true)
        println(stdout, f, " with ", types, " => ", length(allocs), " allocation(s)")
        for (i, a) in enumerate(allocs)
            println(stdout, "  ──────── allocation ", i, " ────────")
            show(stdout, MIME"text/plain"(), a)
            println(stdout)
        end
        flush(stdout)
    end
    return allocs
end

# Test parameters
const μ_perf = 398600.4418
const r0_perf = SVector{3}(7000.0, 0.0, 0.0)
const v0_perf = SVector{3}(0.0, 7.546, 0.0)
const rf_perf = SVector{3}(30000.0, 30000.0, 1000.0)
const vf_perf = SVector{3}(0.0, 3.0, 0.1)
const tof_perf = 86400.0 * 5
const Δt_perf = 86400.0
const throttle_perf = SVector{3}(0.5, 0.3, 0.1)
const mass_perf = 1000.0
const sc_perf = Spacecraft(200.0, 800.0, 0.5, 3000.0)
const sep_perf = SEPSpacecraft(200.0, 800.0, 0.5, 3000.0, 1.495978707e8)
const sail_perf = SolarSail(200.0, 100.0, 0.9, 1.495978707e8)

@testset "Non-Allocating Performance Tests" begin
    @testset "Stumpff Functions — Zero Allocation" begin
        @testset "stumpff_c2" begin
            @test length(checked_allocs(SimsFlanagan.stumpff_c2, (Float64,))) == 0
        end

        @testset "stumpff_c3" begin
            @test length(checked_allocs(SimsFlanagan.stumpff_c3, (Float64,))) == 0
        end
    end

    @testset "safe_norm — Zero Allocation" begin
        @test length(checked_allocs(safe_norm, (SVector{3,Float64},))) == 0
    end

    @testset "Spacecraft Type Construction — Zero Allocation" begin
        @testset "Spacecraft" begin
            @test length(
                checked_allocs(Spacecraft, (Float64, Float64, Float64, Float64))
            ) == 0
        end

        @testset "SEPSpacecraft" begin
            @test length(
                checked_allocs(SEPSpacecraft, (Float64, Float64, Float64, Float64, Float64))
            ) == 0
        end

        @testset "SolarSail" begin
            @test length(checked_allocs(SolarSail, (Float64, Float64, Float64, Float64))) ==
                0
        end
    end

    @testset "Exhaust Velocity — Zero Allocation" begin
        @testset "Spacecraft" begin
            @test length(checked_allocs(exhaust_velocity, (typeof(sc_perf),))) == 0
        end

        @testset "SEPSpacecraft" begin
            @test length(checked_allocs(exhaust_velocity, (typeof(sep_perf),))) == 0
        end

        @testset "SolarSail" begin
            @test length(checked_allocs(exhaust_velocity, (typeof(sail_perf),))) == 0
        end
    end

    @testset "Thrust Computation — Zero Allocation" begin
        @testset "Spacecraft" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.compute_thrust,
                    (typeof(sc_perf), SVector{3,Float64}, Float64),
                ),
            ) == 0
        end

        @testset "SEPSpacecraft" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.compute_thrust,
                    (typeof(sep_perf), SVector{3,Float64}, Float64),
                ),
            ) == 0
        end

        @testset "SolarSail" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.compute_thrust,
                    (typeof(sail_perf), SVector{3,Float64}, Float64),
                ),
            ) == 0
        end
    end

    @testset "Mass Flow Computation — Zero Allocation" begin
        @testset "Spacecraft" begin
            @test length(
                checked_allocs(SimsFlanagan.compute_mass_flow, (typeof(sc_perf), Float64))
            ) == 0
        end

        @testset "SEPSpacecraft" begin
            @test length(
                checked_allocs(SimsFlanagan.compute_mass_flow, (typeof(sep_perf), Float64))
            ) == 0
        end

        @testset "SolarSail" begin
            @test length(
                checked_allocs(SimsFlanagan.compute_mass_flow, (typeof(sail_perf), Float64))
            ) == 0
        end
    end

    @testset "Characteristic Acceleration — Zero Allocation" begin
        @test length(checked_allocs(characteristic_acceleration, (typeof(sail_perf),))) == 0
    end

    @testset "Kepler Propagation — Zero Allocation" begin
        @test length(
            checked_allocs(
                SimsFlanagan.kepler_propagate,
                (SVector{3,Float64}, SVector{3,Float64}, Float64, Float64),
            ),
        ) == 0
    end

    @testset "Segment Propagation — Zero Allocation" begin
        @testset "Spacecraft forward" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.propagate_segment,
                    (
                        SVector{3,Float64},
                        SVector{3,Float64},
                        Float64,
                        SVector{3,Float64},
                        Float64,
                        Float64,
                        typeof(sc_perf),
                    ),
                ),
            ) == 0
        end

        @testset "SEPSpacecraft forward" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.propagate_segment,
                    (
                        SVector{3,Float64},
                        SVector{3,Float64},
                        Float64,
                        SVector{3,Float64},
                        Float64,
                        Float64,
                        typeof(sep_perf),
                    ),
                ),
            ) == 0
        end

        @testset "SolarSail forward" begin
            @test length(
                checked_allocs(
                    SimsFlanagan.propagate_segment,
                    (
                        SVector{3,Float64},
                        SVector{3,Float64},
                        Float64,
                        SVector{3,Float64},
                        Float64,
                        Float64,
                        typeof(sail_perf),
                    ),
                ),
            ) == 0
        end
    end

    @testset "SimsFlanaganOptions — Zero Allocation" begin
        @testset "Default construction" begin
            function create_default_opts()
                return SimsFlanaganOptions()
            end
            @test length(checked_allocs(create_default_opts, ())) == 0
        end
    end
end
