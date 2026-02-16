include("../src/floquet.jl")
using JLD2

include("./parameters.jl")
nα = 1_000 # cyclic ratio
nρ = 1_000 # ρ = ratio T0/P where T0 = pseudo period and P the period
α = range(0.0, 1.0, nα)
ρ = 10.0 .^ range(-1.0, 1.0, nρ)

# WE ALWAYS ASSUME η = 1 for convenience (but without loss of generality)
P, S, S01 = floquet(Γ, γ, α, ρ; abstol=1e-10);
@save "./sims/floquet2d.jld2" Γ γ nα nρ α ρ S S01