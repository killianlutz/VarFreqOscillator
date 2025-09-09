include("../src/polarsolve.jl")
include("../src/polarsolve_refinements.jl")
using JLD2


include("./parameters.jl")
ω = η*Γ/2

# theoretical optimal solution for r = +∞
c1 = bissection(-5, -4, 1e-8) do c
    π + atan(c) - c + 0.5*c*log(1 + c^2) 
end
α0 = (π + atan(c1))/(π + atan(c1) - c1)
ρ0 = 2π/(π + atan(c1) - c1)
z = Point2(α0, ρ0) # pair α = cyclic ratio and ρ = (natural period)/period


# optimal solution for varying 1 <= r < ∞
r = 10.0 .^ range(0.0, 3.0, 300)
minimizers, p = solve_ocp(r, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
τ = tan.(minimizers[3, :])
P = map(t -> T1(t, p), τ)


λp = 2/(1 - γ) # "+" root of quadratic equation
λm = 2/(1 + γ) # "-" root
t = map(τ) do x # length of step at u_min
    num = λm * (λp - x)
    den = λp * (λm - x)
    -log(num/den)/γ ### assumes 0 ≤ γ < 1
end
α1 = (P .- t)./P
ρ1 = (2π ./ P)./ω
Z = [Point2(u, v) for (u, v) in zip(α1, ρ1)] # optimal pairs (α, ρ) as r varies

@save "./sims/optstepfun_amplification.jld2" Γ γ z Z r