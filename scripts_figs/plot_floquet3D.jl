using JLD2
using CairoMakie
CairoMakie.activate!()

function ωα_at_threshold(γ, Γ1)
    # ASSUMES Γ1 = rmax_threshold(γ), i.e. Γ1 is the smallest Γ such that sup(R) >= 1, given γ
    num = 2γ * (π - atan(Γ1))
    den = num + Γ1*(log(1 + γ) - log(1 - γ))
    α_thr = num/den # α^*

    num = 4π
    den = Γ1*(log(1 + Γ1^2) - log(1 - γ^2))
    ω_thr = num/den # ω^*/ω_0

    z_thr = Point2(α_thr, ω_thr)
    return z_thr
end

function fig_floquet3D(; fig=Figure())
    @load "./sims/optstepfun_natfreq.jld2" Γs zs z1s
    @load "./sims/floquet3d.jld2" α ρ γ Γ1 xz y
    ω0 = Γs/2 # natural frequency
    ω01 = Γ1/2 # natural frequency threshold

    fontsize = 20
    textsize = 24

    azimuth = 1.94π

    xlimits = (0.0, 1.0)
    ylimits = (ω01, maximum(ω0))
    zlimits = (1e-1, 3e0)

    xtickvalues = [0, 1/2, 1]
    xticklabels = [L"0", L"1/2", L"1"]
    ytickvalues = [2, 5, 8]
    yticklabels = [L"2", L"5", L"8"]
    ztickvalues = [-1, 0]
    zticklabels = [L"-1", L"0"]

    colors = cgrad([:teal, :orange], length(Γs); categorical=true);

    axs = Axis3(
        fig[1, 1], 
        aspect=:equal, 
        xlabel=L"α",
        xlabelsize=fontsize,
        ylabel=L"ω_0/η",
        ylabelsize=fontsize,
        ylabeloffset=25.0,
        zlabel=L"\log_{10}(ω/ω_0)",
        zlabelsize=fontsize,
        # title=L"\mathrm{Growing~instability~regions}",
        xticks=(xtickvalues, xticklabels),
        yticks=(ytickvalues, yticklabels),
        zticks=(ztickvalues, zticklabels),
        azimuth=azimuth
        )

    # Floquet stability for different values of ω0
    sc = map(ω0, y, colors, zs) do ω, Y, color, z
        xzz, _ = Y
        vertical_slice = map(xzz) do p
            Point3(first(p), ω, log10(last(p)))
        end |> vec
        p = Point3(first(z), ω, log10(last(z)))
        scatter!(axs, vertical_slice; color, markersize=5, transparency=true)
    end

    # Curve of optimal control for unit amplification: r = 1
    p = [Point3(first(z1), ω, log10(last(z1))) for (ω, z1) in zip(ω0, z1s)]
    scatter!(axs, p, color=:purple, marker=:circle, markersize=12)

    # Curve of optimal control for infinite amplification: r = +∞
    light_red = RGBAf(0.9, 0, 0.1, 0.9)
    p = [Point3(first(z), ω, log10(last(z))) for (ω, z) in zip(ω0, zs)]
    scatter!(axs, p, color=light_red, marker=:x, markersize=12)

    # plane highlighting critical ω0 value enabling amplification
    mz, Mz = log10.(zlimits)
    segment = [
        Point3(0.0, ω01, mz), 
        Point3(0.0, ω01, Mz), 
        Point3(1.0, ω01, Mz), 
        Point3(1.0, ω01, mz),
        Point3(0.0, ω01, mz),
        ]
    textposition = segment[4] + Point3(0, 0.3, 0)
    lines!(axs, segment, color=:blue, linestyle=:dash, linewidth=1.0)
    text!(axs, textposition; text=L"ω_{\mathrm{thr}}/\eta", color=:blue, fontsize=textsize)
    
    # limiting parameters as ω0 → ω_{thr}
    z_thr = ωα_at_threshold(γ, Γ1)
    p_thr = Point3(first(z_thr), ω01, log10(last(z_thr)))
    scatter!(axs, p_thr, color=:blue, marker=:circle, markersize=12)

    fig
end

### save 
fig = fig_floquet3D()
save("./plots/floquet3D.png", fig)