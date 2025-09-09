include("../src/polarsolve.jl")
include("../src/polarsolve_refinements.jl")
using JLD2

include("./parameters.jl")

Γs = 4.0:3:20.0

zs = map(Γs) do Γ
    ω = η*Γ/2
    Vmax = η^2 * (1 + Γ^2)/4

    # theoretical solutions
    minimizers, p = solve_ocp(1.0, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
    tp, tm = rootsψ(p, 1e-12)
    P = T1(tm, p)
    
    λp = 2/(1 - γ) # "+" root of quadratic equation
    λm = 2/(1 + γ) # "-" root
    t = begin # length of step at u_min
        num = λm * (λp - tm)
        den = λp * (λm - tm)
        -log(num/den)/γ ### assumes 0 ≤ γ < 1
    end
    α0 = (P - t)/P
    ρ0 = (2π/P)/ω

    z = Point2(α0, ρ0)
end

z1s = map(Γs) do Γ
    ω = η*Γ/2
    Vmax = η^2 * (1 + Γ^2)/4

    # theoretical solutions
    minimizers, p = solve_ocp(1.0, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
    θ = minimizers[3, 1]
    k = minimizers[2, 1]
    T = minimizers[1, 1]
    τ = tan(θ)
    P = T1(τ, p)
    
    λp = 2/(1 - γ) # "+" root of quadratic equation
    λm = 2/(1 + γ) # "-" root
    t = begin # length of step at u_min
        num = λm * (λp - τ)
        den = λp * (λm - τ)
        -log(num/den)/γ ### assumes 0 ≤ γ < 1
    end
    α1 = (P - t)/P
    ρ1 = (2π/P)/ω

    z1 = Point2(α1, ρ1)
end

@save "./sims/optstepfun_natfreq.jld2" γ Γs zs z1s