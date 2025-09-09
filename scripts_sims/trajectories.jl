include("../src/polarsolve.jl")
include("../src/polarsolve_refinements.jl")
include("../src/directsolve.jl")

include("./parameters_.jl")
r = [2.0, 15.0]

# theoretical solution computed numerically
minimizers, p = solve_ocp(r, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
fig = show_results(minimizers, p, r)

### r = 2
i = 1
r_target = r[i]
T, k, θ = minimizers[:, i]
τ = tan(θ)

# check feasible
t, x, dx, v = check_ode(r_target, T, τ, η, Vmin, Vmax);

@save "./sims/trajectory_r2.jld2" r_target t x dx v T  

### r = 15
i = 2
r_target = r[i]
T, k, θ = minimizers[:, i]
τ = tan(θ)

# check feasible
t, x, dx, v = check_ode(r_target, T, τ, η, Vmin, Vmax);

@save "./sims/trajectory_r15.jld2" r_target t x dx v T  