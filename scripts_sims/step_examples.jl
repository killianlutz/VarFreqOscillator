include("../src/stepfn.jl")
include("../scripts_sims/parameters.jl")
using JLD2

# one optimal control example
begin
    r = [2.0]
    minimizers, p = solve_ocp(r, η, Vmin, Vmax; solver=DE.RK4(), nk=20);
    T, k, θ = minimizers[:, 1]
    τ = tan(θ)

    # retrieve periodic feedback control over [0, T] with T = kP
    direct_check = check_ode(first(r), T, τ, η, Vmin, Vmax; dt=T/1e5);
    t = direct_check.t
    v = direct_check.v

    # compute cyclic ratio α and period P
    ds = control_to_steplengths(v, k)
    α = first(ds) + last(ds)
    P = T/k
    ω = 2π/P
end;

# select three examples of step functions, one unstable and two stable
begin
    ω0 = η*Γ/2
    Ω = zeros(3)
    A = zeros(3)
    is_stable = [0, 1, 1]

    Ω[1] = ω/ω0
    A[1] = α

    Ω[2] = 2.5
    A[2] = A[1]

    Ω[3] = Ω[1]
    A[3] = 0.15
end

# phase of initial points on upper half unit circle
ϕs = range(0.0, π, 10)
out = map(Ω, A, eachindex(Ω)) do ρ, α, i
    ω = ρ*ω0
    P = 2π/ω
    tspan = (0.0, 3P)
    u = t -> step_fn(α, ω, t)

    # several trajectories over THREE periods
    Y = sphere_orbit(ϕs, u, tspan; nt=500)
    # control over ONE period only
    t = range(0.0, P, 500)
    v = u.(t)/Vmax

    t, v, Y
end
ts = first.(out)
vs = map(x -> x[2], out)
Ys = last.(out)

@save "./sims/step_parameters.jld2" Ω A is_stable ts vs Ys ϕs ω0
