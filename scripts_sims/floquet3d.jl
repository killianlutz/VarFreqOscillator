include("../src/floquet.jl")
include("./parameters.jl")
using JLD2

Γ1 = rmax_threshold(γ)

@load "./sims/optstepfun_natfreq.jld2" Γs zs

nα = 500
nρ = 500
α = range(0.0, 1.0, nα)
ρ = 10.0 .^ range(-1.0, 1.0, nρ)

floquet_diagrams = map(Γs) do Γ
    floquet(Γ, γ, α, ρ; abstol=1e-10)
end
xz = [[αi, ρj] for αi in α, ρj in ρ]
y = map(floquet_diagrams) do f
    _, _, S01 = f
    keep = S01 .> 1/2 # keep unstable region (bool == 1 <-> unstable)
    xz_keep = xz[keep] # coordinates of unstable region
    S01_keep = S01[keep]
    xz_keep, S01_keep
end

@save "./sims/floquet3d.jld2" α ρ γ Γ1 xz y
