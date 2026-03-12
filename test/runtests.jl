using Test
using SimsFlanagan
using LinearAlgebra
using Random
using StaticArrays

using AllocCheck
using Aqua
using JET

@testset "SimsFlanagan.jl Tests" begin
    include("test_types.jl")
    include("test_propagation.jl")
    include("test_problem.jl")
    include("test_solve.jl")
    include("test_propulsion_types.jl")
    include("test_sundman.jl")
    include("test_pykep_validation.jl")

    # Performance and allocation tests
    include("test_performance.jl")
end

# Differentiability tests — gated behind SIMSFLANAGAN_TEST_DIFF environment variable.
#   "true" / "all"       → all backends (ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote)
#   "ForwardDiff"        → ForwardDiff only (fast smoke-test)
#   comma-separated list → specific backends
#   unset / "false"      → skip
const _DIFF_ENV = get(ENV, "SIMSFLANAGAN_TEST_DIFF", "false")

if _DIFF_ENV ∉ ("false", "")
    using DifferentiationInterface
    using FiniteDiff
    using ForwardDiff

    _run_all = _DIFF_ENV ∈ ("true", "all")
    _requested = _run_all ? Set{String}() : Set(strip.(split(_DIFF_ENV, ",")))
    _need(name) = _run_all || name ∈ _requested

    _backend_list = Tuple{String,Any}[]

    if _need("ForwardDiff")
        push!(_backend_list, ("ForwardDiff", DifferentiationInterface.AutoForwardDiff()))
    end
    if _need("Enzyme")
        using Enzyme
        push!(
            _backend_list,
            (
                "Enzyme",
                DifferentiationInterface.AutoEnzyme(;
                    mode = Enzyme.set_runtime_activity(Enzyme.Forward),
                ),
            ),
        )
    end
    if _need("Mooncake")
        using Mooncake
        push!(
            _backend_list,
            ("Mooncake", DifferentiationInterface.AutoMooncake(; config = nothing)),
        )
    end
    if _need("PolyesterForwardDiff")
        using PolyesterForwardDiff
        push!(
            _backend_list,
            ("PolyesterForwardDiff", DifferentiationInterface.AutoPolyesterForwardDiff()),
        )
    end
    if _need("Zygote")
        using Zygote
        push!(_backend_list, ("Zygote", DifferentiationInterface.AutoZygote()))
    end

    if isempty(_backend_list)
        error(
            "SIMSFLANAGAN_TEST_DIFF=\"$_DIFF_ENV\" did not match any backend. " *
            "Valid names: ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote",
        )
    end

    const _SF_BACKENDS = Tuple(_backend_list)

    @info "Running AD backend tests with: $(join([b[1] for b in _SF_BACKENDS], ", "))"

    @testset "AD Backends" begin
        include("differentiability/test_ad_backends.jl")
    end
else
    @info "Skipping differentiability tests (set SIMSFLANAGAN_TEST_DIFF to enable)"
end

@testset "Aqua Tests" begin
    Aqua.test_all(
        SimsFlanagan;
        ambiguities = (recursive = false),
        deps_compat = (check_extras = false),
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        SimsFlanagan;
        toplevel_logger = nothing,
        target_modules = (@__MODULE__,),
    )
end
