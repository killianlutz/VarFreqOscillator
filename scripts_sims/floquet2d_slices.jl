using JLD2

@load "./sims/floquet2d.jld2" Γ γ nα nρ α ρ S S01
@load "./sims/optstepfun_amplification.jld2" Γ γ z Z r

### r < ∞
k = 91 # index of amplification in vector r

logr = round(log10(r[k]), digits=2)
i = argmin(i -> abs(α[i] - Z[k][1]), eachindex(α))
j = argmin(j -> abs(ρ[j] - Z[k][2]), eachindex(ρ))

cst_freq_slice = S[:, j]
α_opt = α[i]
ρ_opt = ρ[j]
r_val = isinf(k) ? Inf : r[k]
z8 = Point2(α_opt, ρ_opt)

@save "./sims/floquet2d_slices_r_finite.jld2" k r_val α_opt ρ_opt α cst_freq_slice z8


### r = ∞
k = Inf # index of amplification in vector r

logr = Inf
i = argmin(i -> abs(α[i] - z[1]), eachindex(α))
j = argmin(j -> abs(ρ[j] - z[2]), eachindex(ρ))

cst_freq_slice = S[:, j]
α_opt = α[i]
ρ_opt = ρ[j]
r_val = isinf(k) ? Inf : r[k]

@save "./sims/floquet2d_slices_r_infinity.jld2" k r_val α_opt ρ_opt α cst_freq_slice
