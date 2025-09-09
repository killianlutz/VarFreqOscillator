using JLD2

@load "./sims/floquet2d.jld2" Γ γ nα nρ α ρ S S01
@load "./sims/optstepfun_amplification.jld2" Γ γ z Z r

# index of cyclic ratio maximizing spectral radius at constant frequency
iα_ρ = map(1:nρ) do j
    argmax(view(S, :, j))
end 
keep = iα_ρ .> 1 # only instability 
iα_ρ = iα_ρ[keep]
jρ = filter(k -> keep[k] == 1, 1:nρ) # keep frequency whose maximum is > 1

α_ρ = α[iα_ρ]
ρ_keep = ρ[keep]

# locus of maximizers
argmax_sr_cstfreq = [Point2(x, y) for (x, y) in zip(α_ρ, ρ_keep)]
@save "./sims/specrad_argmax.jld2" argmax_sr_cstfreq
