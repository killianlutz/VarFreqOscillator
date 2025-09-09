using GLMakie
import DifferentialEquations as DE
import ForwardDiff: derivative

function Gx(τ, p) 
    -p.a_up*π/p.α + p.offset -p.a_up/p.α*atan((τ - p.a_up)/p.α) + 1/2 * log((p.Vmax*τ^2 - p.damping_rate*τ + 1)/(p.Vmin*τ^2 - p.damping_rate*τ + 1)) + p.a_down/2/p.β * log((τ - p.zp)/(τ - p.zm))
end 

function Gp(τ, p)
    2*Gx(0.0, p) - Gx(τ, p)
end

function R(τ, p)
    τ <= 0 ? exp(Gx(τ, p)) : exp(Gp(τ, p))
end

function limR(γ_up, γ_down)
    limG = 1/2*log((1 + γ_up^2)/(1 - γ_down^2)) - 1/2/γ_down*log((1 + γ_down)/(1 - γ_down)) - atan(1/γ_up)/γ_up - π/2/γ_up
    exp(limG)
end

function limR(p)
    γ_up = sqrt(4*p.Vmax/p.damping_rate^2 - 1)
    γ_down = sqrt(1 - 4*p.Vmin/p.damping_rate^2)
    limR(γ_up, γ_down)
end

function dR(τ, p)
    derivative(t -> R(t, p), τ)
end

function Rinverse(r, p)
    idx = argmin(abs2.(r .- p.r_mesh))
    return p.Rinv[idx]
end

function φ(τ, p)
    atan((τ - p.a_up)/p.α)/(p.α*p.Vmax) - log((τ - p.zp)/(τ - p.zm))/(2*p.β*p.Vmin)
end

function φ(down, up, p)
    π/(p.α*p.Vmax) + φ(up, p) - φ(down, p)
end

function T1(τ, p)
    τ <= 0 ? φ(0.0, τ, p) : φ(τ, 0.0, p)
end

function Tk(r, k, p)
    τk = Rinverse(r^(1/k), p)
    return [k*T1(τk, p), τk]
end

function optimal_strategy(r, odd_k_values, p)
    r_max = limR(p)
    if r_max < 1 && r_max < r # unfeasible
        Tstar = Inf64
        kstar = Inf64
        θstar = Inf64
    else
        time_tangent_pairs = map(k -> Tk(r, k, p), odd_k_values)
        Tstar, idx = findmin(first, time_tangent_pairs)
        kstar = 2*(idx - 1) + 1
        τstar = last(time_tangent_pairs[idx])
        θstar = atan(τstar)
    end
    return [Tstar, kstar, θstar]
end

function show_results(minimizers, p, r)
    T, k, θ = eachrow(minimizers)
    critical_ratio = R(0, p)
    limit_ratio = limR(p)

    fig = Figure(size=(500, 500))
    axs = Axis(
        fig[1, 1], 
        aspect=AxisAspect(1), 
        xlabel=L"r", 
        xscale=log10
        )
    scatter!(axs, r, T, label=L"T^*", color=:purple, markersize=6, marker=:circle)
    scatter!(axs, r, θ, label=L"θ^*", color=:orange, markersize=6, marker=:cross)
    scatter!(axs, r, k, label=L"k^*", color=:teal, markersize=6, marker=:diamond)
    vlines!(axs, critical_ratio, color=:red, label=L"r_c", alpha=1.0, linestyle=:dot, linewidth=1.5)
    vlines!(axs, limit_ratio, color=:black, label=L"\bar{r}", alpha=1.0, linestyle=:dash, linewidth=1.5)
    axislegend(axs, position=:lt, orientation=:horizontal)
    
    # shade regions with counter-intuitive result
    # x = r[(1 .<= r .<= limit_ratio) .* (k .>= 3)]
    # band!(axs, x, -1, 4, color=:gray, alpha=0.1)

    return fig
end

function show_Rinv(p)
    critical_ratio = R(0, p)

    fig = Figure(size=(500, 500))
    axs = Axis(fig[1, 1], aspect=AxisAspect(1), xlabel=L"r", ylabel=L"τ=\tan(θ)")
    lines!(axs, p.r_mesh, p.Rinv, label=L"R^{-1}(r)", color=:blue, linewidth=2)
    vlines!(axs, critical_ratio, color=:red, label=L"r_c", alpha=0.3, linestyle=:dashdot)
    axislegend(axs, position=:lb)

    return fig
end

function ratios_to_tangents(n_mesh, distances, p, solver)
    initial_tangent, asymptotic_τ, dist_to_zm = distances

    # forward for r > r_c (negative tangents)
    τ0 = -abs(initial_tangent)
    r_min = R(τ0, p)
    r_max = R(asymptotic_τ, p)
    right_r, right_Rinv = ratios_to_tangents(n_mesh, (r_min, r_max), τ0, p, true, solver)

    # backward for r < r_c (positive tangents)
    τ0 = abs(initial_tangent)
    r_min = R(p.zm - dist_to_zm, p)
    r_max = R(τ0, p)
    left_r, left_Rinv = ratios_to_tangents(n_mesh, (r_min, r_max), τ0, p, false, solver)

    r_mesh = [left_r; right_r]
    Rinv = [left_Rinv; right_Rinv]
    return (r_mesh, Rinv)
end

function ratios_to_tangents(n_mesh, rspan, τ0, p, forward, solver)
    velocity_sign = forward ? 1.0 : -1.0
    vector_field = (τ, _, _) -> velocity_sign/dR(τ, p)

    saveat = range(rspan..., n_mesh)
    prob = DE.ODEProblem(vector_field, τ0, rspan)
    out = DE.solve(prob, solver; saveat)

    r_mesh = out.t
    Rinv = forward ? out.u : reverse(out.u)
    return (r_mesh, Rinv)
end

function compute_parameters(damping_rate, Vmin, Vmax)
    α = sqrt(1/Vmax - damping_rate^2/(4*Vmax^2))
    β = sqrt(damping_rate^2/(4*Vmin^2) - 1/Vmin)
    a_up = damping_rate/(2*Vmax)
    a_down = damping_rate/(2*Vmin)
    zp = a_down + β
    zm = a_down - β
    offset = - a_up/α * atan(a_up/α) - a_down/2/β * log(zp/zm)

    return (; α, β, a_up, a_down, zp, zm, offset, Vmin, Vmax, damping_rate)
end

function solve_ocp(r, damping_rate, Vmin, Vmax; nk=20, solver=DE.RK4(), n_mesh=5_000)
    p = compute_parameters(damping_rate, Vmin, Vmax)

    # numerically inverse the tangent to ratio mapping R(τ)
    distances = (; initial_tangent=1e-8, asymptotic_τ=-1e5, dist_to_zm=1e-3)
    r_mesh, Rinv = ratios_to_tangents(n_mesh, distances, p, solver)

    p = (; r_mesh, Rinv, p...)

    # solve control problem for all possible ratio values
    admissible_k = (2*i + 1 for i in 0:nk)
    minimizers = zeros(3, length(r)) # row-wise: Tstar, kstar, θstar

    for (i, r_target) in enumerate(r)
        minimizers[:, i] .= optimal_strategy(r_target, admissible_k, p)
    end

    return (minimizers, p)
end


function admissible_ranges(r, m)
    r_max = limR(p)
    threshold = log(r)/log(r_max)
    # because equality in constraint r <= r_max^k :
    thresh_low = 2*ceil(Int, 0.5*(threshold - 1)) + 1
    thresh_up = 2*floor(Int, 0.5*(threshold - 1)) + 1

    k_mesh = range(1.0, 2*m + 1, 1_000)
    k_odd = 2*(0:m) .+ 1
    k_even = 2*(1:m)

    if r_max > 1
        k_lb = max(1, thresh_low)
        k_adm = k_mesh[k_mesh .>= k_lb]
        k_odd = k_odd[k_odd .>= k_lb]
        k_even = k_even[k_even .>= k_lb]

        return (; k_adm, k_odd, k_even)
    else
        k_ub = min(2*m + 1, thresh_up) 

        if r <= r_max && k_ub >= 1
            k_adm = k_mesh[k_mesh .<= k_ub]
            k_odd = k_odd[k_odd .<= k_ub]
            k_even = k_even[k_even .<= k_ub]
            
            return (; k_adm, k_odd, k_even)
        else # unfeasible
            return nothing
        end
    end
end

function J(r, k, p)
    first(Tk(r, k, p))
end

function objective_function(r, m, p)    
    admissible = admissible_ranges(r, m)

    if isnothing(admissible)
        return nothing
    else
        transfer_times = map(k -> J(r, k, p), admissible.k_adm)
        oddint_transfer_times = map(k -> J(r, k, p), admissible.k_odd)
        evenint_transfer_times = map(k -> J(r, k, p), admissible.k_even)
        
        k_argmin = admissible.k_adm[argmin(transfer_times)]
        Tmin = first(Tk(r, k_argmin, p))
        optimum = Point2(k_argmin, Tmin)

        return (;
            admissible...,
            transfer_times,
            oddint_transfer_times,
            evenint_transfer_times,
            optimum,
        )
    end
end

# graph of J(k)
function graph_of_J(p; rs=[0.01, 0.1, 1.0, 10.0, 100.0], m=4)
    fig = Figure(size=(500, 500))
    axs = Axis(fig[1, 1], aspect=AxisAspect(1), xlabel=L"k", ylabel=L"T_k(r)")
    markersize = 14
    marker = :circle
    color = cgrad([:purple, :orange], length(rs))

    foreach(rs, color) do r, c
        data = objective_function(r, m, p)
        if !isnothing(data)
            scatter!(data.k_odd, data.oddint_transfer_times, color=:blue, label=L"\text{odd } k"; markersize, marker)
            # scatter!(data.k_even, data.evenint_transfer_times, color=:black, label=L"\text{even}"; markersize, marker)
            scatter!([data.optimum], color=:red, label=L"\text{argmin}"; markersize, marker=:cross)
            lines!(axs, data.k_adm, data.transfer_times, linewidth=2, label=L"r = %$(r)", alpha=0.5, color=c)
        end
    end

    axislegend(axs, position=:rb, unique=true)
    
    return fig
end