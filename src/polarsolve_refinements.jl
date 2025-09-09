### odd integers candidate for argmin J
function ψ(τ, p)
    return T1(τ, p) + τ * log(R(τ, p))
end

function bissection(f, a0, b0, abstol)
    n = 0
    a = a0
    b = b0
    dist = abs(b - a)

    while dist > abstol && n <= 500
        m = (a + b)/2
        if f(m)*f(a) <= 0
            b = m
        else
            a = m
        end
        dist = abs(b - a)
    end

    return (a + b)/2
end

function rootsψ(p, abstol)
    r_max = limR(p)

    if r_max > 1 # unimodal
        a0 = 0.0
        b0 = a0
        while ψ(b0, p) > 0
            b0 = b0 + (p.zm - b0)/2
        end

        d0 = 0.0
        c0 = d0
        while ψ(c0, p) > 0
            c0 -= 1.0
        end

        τp = bissection(τ -> ψ(τ, p), a0, b0, abstol)
        τm = bissection(τ -> ψ(τ, p), c0, d0, abstol)
        return (; τp, τm)

    else # monotonically decreasing
        a0 = 0.0
        b0 = 0.99*p.zm
        hasroot = ψ(b0, p) <= 0

        if hasroot
            τp = bissection(τ -> ψ(τ, p), a0, b0, abstol)
        else
            τp = nothing
        end
        τm = nothing

        return (; τp, τm)
    end
end

function ρ_bar(r, r_max)
    log(r)/log(r_max)
end

function ρ_pm(r, roots_psi, p; ispositive=true)
    if ispositive && !isnothing(roots_psi.τp)
        return log(r)/log(R(roots_psi.τp, p))

    elseif !ispositive && !isnothing(roots_psi.τm)
        return log(r)/log(R(roots_psi.τm, p))

    else
        return nothing
    end
end

function m_lb(r, r_max)
    if r >= 1 && r_max > 1
        ρ = ρ_bar(r, r_max)
        return ceil(0.5*(ρ - 1))
    else
        return 0
    end
end

function m_ub(r, r_max)
    if r < 1 && r_max < 1
        ρ = ρ_bar(r, r_max)
        return floor(0.5*(ρ - 1))
    else
        return Inf64
    end
end

function m_lr_pm(r, roots_psi, p; ispositive=true, left=true)
    ρ = ρ_pm(r, roots_psi, p; ispositive)
    if isnothing(ρ)
        return Inf64
    else
        m = 0.5*(ρ - 1)
        m = left ? floor(m) : ceil(m)
        return m
    end
end

function candidate_integers(rs, p; abstol=1e-8)
    r_max = limR(p)
    roots_psi = rootsψ(p, abstol)
    greater_than_one = rs .>= 1
    smaller_than_one = rs .< 1
    rs_gt1 = rs[greater_than_one]
    rs_st1 = rs[smaller_than_one]


    M = similar(rs, 4, length(rs)) # columns: m_min, m_max, m_left, m_right

    M[1, :] .= map(r -> m_lb(r, r_max), rs)
    M[2, :] .= map(r -> m_ub(r, r_max), rs)

    M[3, greater_than_one] .= map(rs_gt1) do r
        m_lr_pm(r, roots_psi, p; ispositive=false, left=true)
    end
    M[4, greater_than_one] .= map(rs_gt1) do r
        m_lr_pm(r, roots_psi, p; ispositive=false, left=false)
    end

    M[3, smaller_than_one] .= map(rs_st1) do r
        m_lr_pm(r, roots_psi, p; ispositive=true, left=true)
    end
    M[4, smaller_than_one] .= map(rs_st1) do r
        m_lr_pm(r, roots_psi, p; ispositive=true, left=false)
    end
    
    M[M .< 0] .= 0.0
    return M
end

# critical value of Vmax making 1 <= r_max
function F(x)
    x*log(1 + x^2)/2 + atan(x) - π
end

function dF(x)
    log(1 + x^2)/2 + 1
end

function bounds_critical_rmax(γ_lb, y)
    γ_ub = (γ_lb*dF(γ_lb) - F(γ_lb)) / (dF(γ_lb) - y)

    slope = (F(γ_lb) - F(γ_ub)) / (γ_lb - γ_ub)
    γ_lb = (γ_lb*slope - F(γ_lb)) / (slope - y)

    return (γ_lb, γ_ub)
end

function bounds_critical_rmax(y; n_repeat=2)
    γ_lb = 1.0
    γ_ub = Inf64
    for _ in 1:n_repeat
        γ_lb, γ_ub = bounds_critical_rmax(γ_lb, y)
    end

    return (γ_lb, γ_ub)
end


function critical_Vmax(η, Vmin)
    γ_down = sqrt(1 - 4Vmin/η^2)

    if isapprox(γ_down, 0.0)
        y = 1.0
    elseif isapprox(γ_down, 1.0)
        y = log(2)
    else
        y = log(1 - γ_down^2)/2 + log((1 + γ_down)/(1 - γ_down))/2/γ_down
    end

    bounds = bounds_critical_rmax(y)
    Vmax_bounds = map(bounds) do γ_up
        η^2 * (γ_up^2 + 1)/4
    end

    return Vmax_bounds
end