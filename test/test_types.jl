@testset "Types" begin

    @testset "Spacecraft (constant thrust)" begin
        sc = Spacecraft(1000.0, 0.5, 3000.0)
        @test sc.mass == 1000.0
        @test sc.thrust == 0.5
        @test sc.isp == 3000.0
        @test sc isa AbstractSpacecraft

        # Test exhaust velocity
        vex = SimsFlanagan.exhaust_velocity(sc)
        @test vex ≈ 3000.0 * 9.80665  # Isp * g0

        # Test compute_thrust (constant regardless of position)
        r1 = SVector{3}(1.0e8, 0.0, 0.0)
        r2 = SVector{3}(2.0e8, 0.0, 0.0)
        @test SimsFlanagan.compute_thrust(sc, r1, 1.0) == sc.thrust
        @test SimsFlanagan.compute_thrust(sc, r2, 1.0) == sc.thrust
        @test SimsFlanagan.compute_thrust(sc, r1, 0.5) == 0.5 * sc.thrust

        # Test propellant consumption
        @test SimsFlanagan.has_propellant_consumption(sc) == true
        @test SimsFlanagan.compute_mass_flow(sc, 0.5) ≈ 0.5 / vex
    end

    @testset "SEPSpacecraft" begin
        # SEP with 0.5 N thrust at 1 AU, 3000s Isp
        r_1AU = 1.495978707e8  # km
        sep = SEPSpacecraft(1000.0, 0.5, 3000.0, r_1AU)
        @test sep.mass == 1000.0
        @test sep.thrust_ref == 0.5
        @test sep.isp == 3000.0
        @test sep.r_ref == r_1AU
        @test sep isa AbstractSpacecraft

        # Test exhaust velocity
        vex = SimsFlanagan.exhaust_velocity(sep)
        @test vex ≈ 3000.0 * 9.80665

        # Test thrust at reference distance (1 AU)
        r_ref = SVector{3}(r_1AU, 0.0, 0.0)
        @test SimsFlanagan.compute_thrust(sep, r_ref, 1.0) ≈ sep.thrust_ref

        # Test thrust at 2 AU (should be 1/4 due to inverse square)
        r_2AU = SVector{3}(2 * r_1AU, 0.0, 0.0)
        @test SimsFlanagan.compute_thrust(sep, r_2AU, 1.0) ≈ sep.thrust_ref / 4

        # Test thrust at 0.5 AU (should be 4x)
        r_half_AU = SVector{3}(0.5 * r_1AU, 0.0, 0.0)
        @test SimsFlanagan.compute_thrust(sep, r_half_AU, 1.0) ≈ sep.thrust_ref * 4

        # Test throttle scaling
        @test SimsFlanagan.compute_thrust(sep, r_ref, 0.5) ≈ 0.5 * sep.thrust_ref

        # Test propellant consumption
        @test SimsFlanagan.has_propellant_consumption(sep) == true

        # Validation
        @test_throws ArgumentError SEPSpacecraft(0.0, 0.5, 3000.0)  # zero mass
        @test_throws ArgumentError SEPSpacecraft(1000.0, -0.5, 3000.0)  # negative thrust
    end

    @testset "SolarSail" begin
        # Solar sail: 100 kg, 1000 m² sail, 0.9 reflectivity
        sail = SolarSail(100.0, 1000.0, 0.9)
        @test sail.mass == 100.0
        @test sail.area == 1000.0
        @test sail.reflectivity == 0.9
        @test sail isa AbstractSpacecraft

        # Test characteristic acceleration
        # a_c = (1 + η) × P_sun × A / (c × m)
        # a_c = (1 + 0.9) × 1361 × 1000 / (299792458 × 100)
        a_c = SimsFlanagan.characteristic_acceleration(sail)
        expected_a_c = (1 + 0.9) * 1361.0 * 1000.0 / (299792458.0 * 100.0)
        @test a_c ≈ expected_a_c

        # Test thrust computation (F = m × a)
        r_1AU = SVector{3}(1.495978707e8, 0.0, 0.0)
        thrust_1AU = SimsFlanagan.compute_thrust(sail, r_1AU, 1.0)
        @test thrust_1AU ≈ sail.mass * a_c

        # Test thrust at 2 AU (1/4 due to inverse square)
        r_2AU = SVector{3}(2 * 1.495978707e8, 0.0, 0.0)
        thrust_2AU = SimsFlanagan.compute_thrust(sail, r_2AU, 1.0)
        @test thrust_2AU ≈ thrust_1AU / 4

        # Test NO propellant consumption
        @test SimsFlanagan.has_propellant_consumption(sail) == false
        @test SimsFlanagan.compute_mass_flow(sail, 100.0) == 0.0

        # Test exhaust velocity (returns Inf for sails)
        @test SimsFlanagan.exhaust_velocity(sail) == Inf

        # Validation
        @test_throws ArgumentError SolarSail(100.0, 1000.0, 1.5)  # reflectivity > 1
        @test_throws ArgumentError SolarSail(100.0, -1000.0, 0.9)  # negative area
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
