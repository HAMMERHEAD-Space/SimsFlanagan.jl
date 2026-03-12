using AllocCheck

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
            allocs_vec = check_allocs(SimsFlanagan.stumpff_c2, (Float64,))
            @test length(allocs_vec) == 0
        end

        @testset "stumpff_c3" begin
            allocs_vec = check_allocs(SimsFlanagan.stumpff_c3, (Float64,))
            @test length(allocs_vec) == 0
        end
    end

    @testset "safe_norm — Zero Allocation" begin
        allocs_vec = check_allocs(safe_norm, (SVector{3,Float64},))
        @test length(allocs_vec) == 0
    end

    @testset "Spacecraft Type Construction — Zero Allocation" begin
        @testset "Spacecraft" begin
            allocs_vec = check_allocs(Spacecraft, (Float64, Float64, Float64, Float64))
            @test length(allocs_vec) == 0
        end

        @testset "SEPSpacecraft" begin
            allocs_vec = check_allocs(
                SEPSpacecraft,
                (Float64, Float64, Float64, Float64, Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "SolarSail" begin
            allocs_vec = check_allocs(
                SolarSail,
                (Float64, Float64, Float64, Float64),
            )
            @test length(allocs_vec) == 0
        end
    end

    @testset "Exhaust Velocity — Zero Allocation" begin
        @testset "Spacecraft" begin
            allocs_vec = check_allocs(exhaust_velocity, (typeof(sc_perf),))
            @test length(allocs_vec) == 0
        end

        @testset "SEPSpacecraft" begin
            allocs_vec = check_allocs(exhaust_velocity, (typeof(sep_perf),))
            @test length(allocs_vec) == 0
        end

        @testset "SolarSail" begin
            allocs_vec = check_allocs(exhaust_velocity, (typeof(sail_perf),))
            @test length(allocs_vec) == 0
        end
    end

    @testset "Thrust Computation — Zero Allocation" begin
        @testset "Spacecraft" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_thrust,
                (typeof(sc_perf), SVector{3,Float64}, Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "SEPSpacecraft" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_thrust,
                (typeof(sep_perf), SVector{3,Float64}, Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "SolarSail" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_thrust,
                (typeof(sail_perf), SVector{3,Float64}, Float64),
            )
            @test length(allocs_vec) == 0
        end
    end

    @testset "Mass Flow Computation — Zero Allocation" begin
        @testset "Spacecraft" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_mass_flow,
                (typeof(sc_perf), Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "SEPSpacecraft" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_mass_flow,
                (typeof(sep_perf), Float64),
            )
            @test length(allocs_vec) == 0
        end

        @testset "SolarSail" begin
            allocs_vec = check_allocs(
                SimsFlanagan.compute_mass_flow,
                (typeof(sail_perf), Float64),
            )
            @test length(allocs_vec) == 0
        end
    end

    @testset "Characteristic Acceleration — Zero Allocation" begin
        allocs_vec = check_allocs(characteristic_acceleration, (typeof(sail_perf),))
        @test length(allocs_vec) == 0
    end

    @testset "Kepler Propagation — Zero Allocation" begin
        allocs_vec = check_allocs(
            SimsFlanagan.kepler_propagate,
            (SVector{3,Float64}, SVector{3,Float64}, Float64, Float64),
        )
        @test length(allocs_vec) == 0
    end

    @testset "Segment Propagation — Zero Allocation" begin
        @testset "Spacecraft forward" begin
            allocs_vec = check_allocs(
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
            )
            @test length(allocs_vec) == 0
        end

        @testset "SEPSpacecraft forward" begin
            allocs_vec = check_allocs(
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
            )
            @test length(allocs_vec) == 0
        end

        @testset "SolarSail forward" begin
            allocs_vec = check_allocs(
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
            )
            @test length(allocs_vec) == 0
        end
    end

    @testset "SimsFlanaganOptions — Zero Allocation" begin
        @testset "Default construction" begin
            # Test that default option creation doesn't allocate
            function create_default_opts()
                return SimsFlanaganOptions()
            end
            allocs_vec = check_allocs(create_default_opts, ())
            @test length(allocs_vec) == 0
        end
    end
end
