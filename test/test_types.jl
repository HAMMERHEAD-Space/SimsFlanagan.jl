@testset "Types" begin

    @testset "Spacecraft" begin
        sc = Spacecraft(1000.0, 0.5, 3000.0)
        @test sc.mass == 1000.0
        @test sc.thrust == 0.5
        @test sc.isp == 3000.0

        # Test exhaust velocity
        vex = SimsFlanagan.exhaust_velocity(sc)
        @test vex ≈ 3000.0 * 9.80665  # Isp * g0
    end

    @testset "SimsFlanaganOptions" begin
        # Default options
        opts = SimsFlanaganOptions()
        @test opts.n_segments == 10
        @test opts.n_fwd == 5
        @test opts.tol == 1e-8
        @test opts.max_iter == 1000
        @test opts.verbosity == 1

        # Custom options
        opts2 = SimsFlanaganOptions(n_segments = 20, n_fwd = 8, tol = 1e-10)
        @test opts2.n_segments == 20
        @test opts2.n_fwd == 8
        @test opts2.tol == 1e-10

        # Validation
        @test_throws ArgumentError SimsFlanaganOptions(n_segments = 0)
        @test_throws ArgumentError SimsFlanaganOptions(n_segments = 10, n_fwd = 0)
        @test_throws ArgumentError SimsFlanaganOptions(n_segments = 10, n_fwd = 11)
        @test_throws ArgumentError SimsFlanaganOptions(tol = -1.0)
    end

    @testset "Vector conversion utilities" begin
        throttles = [SVector{3}(1.0, 2.0, 3.0), SVector{3}(4.0, 5.0, 6.0)]

        x = SimsFlanagan.throttles_to_vector(throttles)
        @test x == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

        throttles_back = SimsFlanagan.vector_to_throttles(x, 2)
        @test throttles_back == throttles
    end

end

