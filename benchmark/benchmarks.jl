using BenchmarkTools
using SimsFlanagan
using StaticArrays

const SUITE = BenchmarkGroup()

# Test parameters
const μ_bench = 398600.4418
const r0_bench = SVector{3}(7000.0, 0.0, 0.0)
const v0_bench = SVector{3}(0.0, 7.546, 0.0)
const rf_bench = SVector{3}(30000.0, 30000.0, 1000.0)
const vf_bench = SVector{3}(0.0, 3.0, 0.1)
const tof_bench = 86400.0 * 5
const Δt_bench = 86400.0
const throttle_bench = SVector{3}(0.5, 0.3, 0.1)
const mass_bench = 1000.0
const sc_bench = Spacecraft(200.0, 800.0, 0.5, 3000.0)
const sep_bench = SEPSpacecraft(200.0, 800.0, 0.5, 3000.0, 1.495978707e8)
const sail_bench = SolarSail(200.0, 100.0, 0.9, 1.495978707e8)

# Stumpff functions
SUITE["Stumpff"] = BenchmarkGroup()
SUITE["Stumpff"]["c2_positive"] = @benchmarkable SimsFlanagan.stumpff_c2(1.5)
SUITE["Stumpff"]["c2_negative"] = @benchmarkable SimsFlanagan.stumpff_c2(-1.5)
SUITE["Stumpff"]["c2_small"] = @benchmarkable SimsFlanagan.stumpff_c2(1e-8)
SUITE["Stumpff"]["c3_positive"] = @benchmarkable SimsFlanagan.stumpff_c3(1.5)
SUITE["Stumpff"]["c3_negative"] = @benchmarkable SimsFlanagan.stumpff_c3(-1.5)
SUITE["Stumpff"]["c3_small"] = @benchmarkable SimsFlanagan.stumpff_c3(1e-8)

# safe_norm
SUITE["Utils"] = BenchmarkGroup()
SUITE["Utils"]["safe_norm"] = @benchmarkable safe_norm($throttle_bench)
SUITE["Utils"]["safe_norm_zero"] = @benchmarkable safe_norm(SVector{3}(0.0, 0.0, 0.0))

# Kepler propagation
SUITE["Kepler"] = BenchmarkGroup()
SUITE["Kepler"]["propagate_circular"] =
    @benchmarkable SimsFlanagan.kepler_propagate($r0_bench, $v0_bench, 3600.0, $μ_bench)
SUITE["Kepler"]["propagate_long"] =
    @benchmarkable SimsFlanagan.kepler_propagate($r0_bench, $v0_bench, $Δt_bench, $μ_bench)

# Segment propagation
SUITE["Segment"] = BenchmarkGroup()
SUITE["Segment"]["Spacecraft"] = @benchmarkable SimsFlanagan.propagate_segment(
    $r0_bench,
    $v0_bench,
    $mass_bench,
    $throttle_bench,
    $Δt_bench,
    $μ_bench,
    $sc_bench,
)
SUITE["Segment"]["SEPSpacecraft"] = @benchmarkable SimsFlanagan.propagate_segment(
    $r0_bench,
    $v0_bench,
    $mass_bench,
    $throttle_bench,
    $Δt_bench,
    $μ_bench,
    $sep_bench,
)
SUITE["Segment"]["SolarSail"] = @benchmarkable SimsFlanagan.propagate_segment(
    $r0_bench,
    $v0_bench,
    $mass_bench,
    $throttle_bench,
    $Δt_bench,
    $μ_bench,
    $sail_bench,
)

# Leg propagation
const throttles_bench_4 = [SVector{3}(0.1*i, 0.05*i, 0.02*i) for i = 1:4]
const Δt_segments_bench = fill(Δt_bench/4, 4)
SUITE["Leg"] = BenchmarkGroup()
SUITE["Leg"]["Spacecraft_4seg"] = @benchmarkable SimsFlanagan.propagate_leg(
    $r0_bench,
    $v0_bench,
    $mass_bench,
    $throttles_bench_4,
    $Δt_segments_bench,
    $μ_bench,
    $sc_bench,
)

# Problem construction
SUITE["Problem"] = BenchmarkGroup()
SUITE["Problem"]["construction_4seg"] = @benchmarkable simsflanagan_problem(
    $r0_bench,
    $v0_bench,
    $rf_bench,
    $vf_bench,
    $tof_bench,
    $μ_bench,
    $sc_bench;
    n_segments = 4,
)
SUITE["Problem"]["construction_10seg"] = @benchmarkable simsflanagan_problem(
    $r0_bench,
    $v0_bench,
    $rf_bench,
    $vf_bench,
    $tof_bench,
    $μ_bench,
    $sc_bench;
    n_segments = 10,
)

# Mismatch computation
const prob_bench_4 = simsflanagan_problem(
    r0_bench,
    v0_bench,
    rf_bench,
    vf_bench,
    tof_bench,
    μ_bench,
    sc_bench;
    n_segments = 4,
    verbosity = 0,
)
const prob_bench_10 = simsflanagan_problem(
    r0_bench,
    v0_bench,
    rf_bench,
    vf_bench,
    tof_bench,
    μ_bench,
    sc_bench;
    n_segments = 10,
    verbosity = 0,
)
const throttles_bench_mismatch = [SVector{3}(0.1, 0.1, 0.0) for _ = 1:4]
const throttles_bench_10 = [SVector{3}(0.1*i/10, 0.05*i/10, 0.02*i/10) for i = 1:10]
SUITE["Mismatch"] = BenchmarkGroup()
SUITE["Mismatch"]["compute_4seg"] =
    @benchmarkable compute_mismatch($prob_bench_4, $throttles_bench_mismatch)
SUITE["Mismatch"]["compute_10seg"] =
    @benchmarkable compute_mismatch($prob_bench_10, $throttles_bench_10)

# Scaled mismatch constraints (used in solve)
SUITE["Constraints"] = BenchmarkGroup()
SUITE["Constraints"]["scaled_mismatch_4seg"] =
    @benchmarkable scaled_mismatch_constraints($prob_bench_4, $throttles_bench_mismatch)
SUITE["Constraints"]["scaled_mismatch_10seg"] =
    @benchmarkable scaled_mismatch_constraints($prob_bench_10, $throttles_bench_10)
SUITE["Constraints"]["throttle_magnitude_4seg"] =
    @benchmarkable SimsFlanagan.throttle_magnitude_constraints($throttles_bench_mismatch)
SUITE["Constraints"]["throttle_magnitude_10seg"] =
    @benchmarkable SimsFlanagan.throttle_magnitude_constraints($throttles_bench_10)

# Sundman transformation
SUITE["Sundman"] = BenchmarkGroup()
SUITE["Sundman"]["segments_uniform"] = @benchmarkable SimsFlanagan.compute_sundman_segments(
    $r0_bench,
    $rf_bench,
    $tof_bench,
    10;
    c = 0.0,
)
SUITE["Sundman"]["segments_c1"] = @benchmarkable SimsFlanagan.compute_sundman_segments(
    $r0_bench,
    $rf_bench,
    $tof_bench,
    10;
    c = 1.0,
)
SUITE["Sundman"]["leg_times"] = @benchmarkable SimsFlanagan.compute_sundman_leg_times(
    $r0_bench,
    $rf_bench,
    $tof_bench,
    5,
    5;
    c = 1.0,
)

# Spacecraft type operations
SUITE["Spacecraft"] = BenchmarkGroup()
SUITE["Spacecraft"]["mass"] = @benchmarkable mass($sc_bench)
SUITE["Spacecraft"]["exhaust_velocity"] = @benchmarkable exhaust_velocity($sc_bench)
SUITE["Spacecraft"]["compute_thrust"] =
    @benchmarkable SimsFlanagan.compute_thrust($sc_bench, $r0_bench, 0.5)

# SEP operations
SUITE["SEP"] = BenchmarkGroup()
SUITE["SEP"]["compute_thrust_1AU"] = @benchmarkable SimsFlanagan.compute_thrust(
    $sep_bench,
    SVector{3}(1.495978707e8, 0.0, 0.0),
    0.5,
)
SUITE["SEP"]["compute_thrust_2AU"] = @benchmarkable SimsFlanagan.compute_thrust(
    $sep_bench,
    SVector{3}(2.99e8, 0.0, 0.0),
    0.5,
)

# Solar sail operations
SUITE["SolarSail"] = BenchmarkGroup()
SUITE["SolarSail"]["characteristic_acceleration"] =
    @benchmarkable characteristic_acceleration($sail_bench)
SUITE["SolarSail"]["compute_thrust_1AU"] = @benchmarkable SimsFlanagan.compute_thrust(
    $sail_bench,
    SVector{3}(1.495978707e8, 0.0, 0.0),
    0.5,
)

# Total ΔV computation
SUITE["TotalDV"] = BenchmarkGroup()
SUITE["TotalDV"]["compute_4seg"] =
    @benchmarkable compute_total_Δv($throttles_bench_4, $Δt_segments_bench, $sc_bench)

# Initial guess generation
SUITE["InitialGuess"] = BenchmarkGroup()
SUITE["InitialGuess"]["zero"] =
    @benchmarkable SimsFlanagan.generate_initial_guess($prob_bench_4, ZeroGuess())
SUITE["InitialGuess"]["random"] =
    @benchmarkable SimsFlanagan.generate_initial_guess($prob_bench_4, RandomGuess())
SUITE["InitialGuess"]["constant"] = @benchmarkable SimsFlanagan.generate_initial_guess(
    $prob_bench_4,
    ConstantGuess(direction = [1, 0, 0], magnitude = 0.5),
)
SUITE["InitialGuess"]["lambert"] =
    @benchmarkable SimsFlanagan.generate_initial_guess($prob_bench_4, LambertGuess())

# Full solve benchmarks
SUITE["Solve"] = BenchmarkGroup()
SUITE["Solve"]["4seg_zero_guess"] = @benchmarkable solve(
    $prob_bench_4;
    initial_guess_strategy = ZeroGuess(),
    max_iter = 100,
)
SUITE["Solve"]["4seg_random_guess"] = @benchmarkable solve(
    $prob_bench_4;
    initial_guess_strategy = RandomGuess(),
    max_iter = 100,
)
SUITE["Solve"]["4seg_lambert_guess"] = @benchmarkable solve(
    $prob_bench_4;
    initial_guess_strategy = LambertGuess(),
    max_iter = 100,
)
SUITE["Solve"]["10seg_lambert_guess"] = @benchmarkable solve(
    $prob_bench_10;
    initial_guess_strategy = LambertGuess(),
    max_iter = 100,
)

# Helper functions (used in solve)
SUITE["Helpers"] = BenchmarkGroup()
const x_bench_flat = SimsFlanagan.throttles_to_vector(throttles_bench_mismatch)
SUITE["Helpers"]["throttles_to_vector_4"] =
    @benchmarkable SimsFlanagan.throttles_to_vector($throttles_bench_mismatch)
SUITE["Helpers"]["vector_to_throttles_4"] =
    @benchmarkable SimsFlanagan.vector_to_throttles($x_bench_flat, 4)
SUITE["Helpers"]["get_canonical_scales"] =
    @benchmarkable SimsFlanagan.get_canonical_scales($prob_bench_4)
SUITE["Helpers"]["compute_segment_masses_4"] =
    @benchmarkable SimsFlanagan.compute_segment_masses(
        $prob_bench_4,
        $throttles_bench_mismatch,
        $Δt_segments_bench,
    )
