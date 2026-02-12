include("../src/polarsolve.jl")
include("../src/polarsolve_refinements.jl")
include("../src/directsolve.jl")

function step_driven_oscillator(y::T, p, t) where T <: AbstractVector
    u = p 
    return [y[2], -u(t)*y[1] -η*y[2]]
end

function sdo(Y::T, p, t) where T <: AbstractArray
    mapslices(y -> step_driven_oscillator(y, p, t), Y; dims=1)
end

function step_fn(α, ω, t)
    # duty cycle α
    # period P
    P = 2π/ω
    m = mod(t, P)
    if m > α*P || isapprox(t, P)
        u = Vmin
    else
        u = Vmax
    end

    return u
end

function control_to_steplengths(vs, k)
    jumps = diff(vs)
    idx = filter(x -> jumps[x] != 0, eachindex(jumps))
    idx = [1; idx; length(vs)]

    # take into account the number of swings
    dts = k*diff(idx)/length(v)

    ds1 = first(dts)
    ds2 = dts[2]
    ds3 = last(dts)
    return ds1, ds2, ds3
end

function sphere_orbit(ϕs, p, tspan; nt=500)
    y0s = map(ϕ -> Point2(cos(ϕ), sin(ϕ)), ϕs)
    Y0s = reduce(hcat, y0s)

    saveat = range(tspan..., nt)
    ode = DE.ODEProblem(sdo, Y0s, tspan, p; saveat)
    Y = DE.solve(ode, DE.RK4()).u
    Y = cat(Y...; dims=3)
end