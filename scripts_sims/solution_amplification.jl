include("../src/polarsolve.jl")
include("../src/polarsolve_refinements.jl")
using JLD2


include("./parameters_.jl")

# theoretical solution (T, k, θ) against amplification r
r = 10.0 .^ range(-4.0, 4.0, 200)
minimizers, p = solve_ocp(r, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
limit_ratio = limR(p) # maximal amplification per swing
critical_ratio = R(0.0, p) # amplification for constant control equal to max value

@save "./sims/solution_amplification.jld2" Γ γ minimizers p r limit_ratio critical_ratio