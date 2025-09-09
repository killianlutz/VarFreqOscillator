using LinearAlgebra
using Base.Threads

function maxflow(Γ, t)
    ω = Γ/2
    # C takes initial conditions to coefficients
    c11 = 1.0
    c12 = 0.0
    c21 = 1.0/(2ω)
    c22 = 1.0/ω
    C = [c11 c12; c21 c22]

    # B takes coefficients to solution at t
    c = cos(ω*t)
    s = sin(ω*t)
    E = exp(-t/2)
    b11 = E*c
    b12 = E*s
    b21 = E*(-c/2 - ω*s)
    b22 = E*(-s/2 + ω*c)
    B = [b11 b12; b21 b22]

    return B*C
end

function minflow(γ, t)
    λp = (-1 + γ)/2
    λm = (-1 - γ)/2

    # M takes initial conditions at time t to coefficients
    m11 = -λm/γ
    m12 = 1.0/γ
    m21 = λp/γ
    m22 = -1.0/γ
    M = [m11 m12; m21 m22]
    
    # A takes coefficients at time t to solution at time t + δ
    Ep = exp(λp*t)
    Em = exp(λm*t)

    a11 = Ep
    a12 = Em
    a21 = λp*Ep
    a22 = λm*Em
    A = [a11 a12; a21 a22]

    return A*M
end

function monodromy(α, P, Γ, γ)
    #### u = (1 + Γ^2)/4 maximal value
    t = α*P
    ϕ1 = maxflow(Γ, t)
    
    #### u = (1 - γ^2)/4 minimal value
    t = (1 - α)*P
    ϕ2 = minflow(γ, t)

    #### monodromy 
    # flow at time P = flow at (1-α)P times flow at αP
    return ϕ2 * ϕ1
end

function spectrum(ϕ)
    t = tr(ϕ)
    d = det(ϕ)
    Δ = t^2 - 4d

    δ = Δ ≥ 0 ? sqrt(Δ) : im*sqrt(-Δ)
    [(t + δ)/2, (t - δ)/2]
end

function isunstable(ϕ; abstol=0.0)
    vals = spectrum(ϕ)
    ρ = maximum(abs, vals)

    is_unstable = missing
    if ρ > 1.0 + abstol
        is_unstable = true 
    elseif ρ < 1.0 - abstol
        is_unstable = false
    else
        # eigenvalues of modulus = 1
        # since Real and Order 2, eigenvalues
        # are necessarily equal
        # -> check semi-simplicity
        λ = vals[1]
        is_unstable = !is_semisimple(ϕ, λ; abstol)
    end

    return (ρ, is_unstable)
end


function kernel_dimension(A; abstol=0.0)
    size(nullspace(A; atol=abstol), 2)
end

function is_semisimple(A, λ; abstol)
    is_ssimple = missing

    B = A .- λ.*I(2)
    geo_mult = kernel_dimension(B; abstol)

    if geo_mult == 2
        # 2 <= d_g <= d_a <= 2
        is_ssimple = true
    elseif geo_mult == 1
        # need only check B^2 cause A order 2
        alg_mult = kernel_dimension(B^2; abstol)
        
        if alg_mult == geo_mult
            # d_g = 1 = d_a 
            is_ssimple = true
        else
            # 1 = d_g < d_a
            is_ssimple = false
        end
    else
        is_ssimple = missing 
    end

    return is_ssimple
end

function floquet(Γ, γ, α, ρ; abstol=0.0)
    nα = length(α)
    nρ = length(ρ)

    ω = Γ/2  # natural pseudo pulsation ω0
    T = 2π/ω # natural pseudo period
    P = T ./ ρ

    x = vec([[αi, Pj] for αi in α, Pj in P])
    S = vec(zeros(nα, nρ))
    S01 = vec(zeros(Bool, nα, nρ))

    @threads for i in eachindex(S)
        αi, Pi = x[i]
        ϕ = monodromy(αi, Pi, Γ, γ)
        ρ, is_unstable = isunstable(ϕ; abstol)
        S[i] = ρ
        if ismissing(is_unstable)
            S01[i] = false
        else
            S01[i] = is_unstable
        end
    end

    S = reshape(S, nα, nρ)
    S01 = reshape(S01, nα, nρ)

    return (; P, S, S01)
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

function log_rmax(γ, Γ)
    if 0 < γ < 1 && Γ > 0
        return -0.5*log( (1 - γ^2)/(1 + Γ^2) ) - 0.5*log( (1 + γ)/(1 - γ) )/γ - atan(1/Γ)/Γ - π/(2Γ)
    else
        return Inf64
    end
end

function rmax_threshold(γ)
    bissection(x -> log_rmax(γ, x), 1e-3, 1e2, 1e-8)
end