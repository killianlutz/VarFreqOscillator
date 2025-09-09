import DifferentialEquations as DE

function feedback(X, p)
    x, dx = X
    τ, Vmin, Vmax = p
    τ_min, τ_max = extrema((τ, 0))

    tangent = -x/dx
    (τ_min <= tangent < τ_max) ? Vmin : Vmax
end

function velocity(dx, x, p, _)
    dx[1] = x[2]
    dx[2] = - p.η * x[2] - feedback(x, p.Vfb) * x[1]
end

function check_ode(r, T, τ, η, Vmin, Vmax; dt=1.0/1000)
    Vfb = (τ, Vmin, Vmax)
    p = (; η, Vfb)

    target_angle = 3.0
    x0 = (r > 1) ? -target_angle/r : -1.0
    X0 = [x0, 0.0]
    tspan = (0.0, T)
    nt = ceil(Int, T/dt + 1)
    saveat = range(tspan..., nt)

    prob = DE.ODEProblem(velocity, X0, tspan)
    out = DE.solve(prob, DE.RK4(); saveat, p)
    
    x, dx = map((first, last)) do fun
        fun.(out.u)
    end
    v = map(out.u) do X
        feedback(X, Vfb)
    end

    return (; t=out.t, x, dx, v)
end
