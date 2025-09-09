using JLD2
using CairoMakie
CairoMakie.activate!()

function fig_floquet2D(; fig=Figure())
    @load "./sims/floquet2d.jld2" Γ γ nα nρ α ρ S S01
    @load "./sims/optstepfun_amplification.jld2" Γ γ z Z r
    @load "./sims/specrad_argmax.jld2" argmax_sr_cstfreq

    fontsize = 20

    xlimits = (0.0, 1.0)
    ylimits = (1e-1, 3e0)
    Sthreshold = 1.0

    xtickvalues = [0, 1/4, 1/2, 3/4, 1]
    xticklabels = [L"0", L"1/4", L"1/2", L"3/4", L"1"]
    ytickvalues = [1e-1, 1]
    yticklabels = [L"10^{%$(n)}" for n in -1:0]

    hlinesfreq = reverse([1/2, 1, 2])
    hllabels = reverse([L"1/2", L"1", L"2"])
    hllinestyles = [:dot, :dash, :dashdot]

    axs = Axis(
        fig[1, 1], 
        aspect=AxisAspect(1), 
        xlabel=L"α",
        xlabelsize=fontsize,
        ylabel=L"ω/ω_0",
        ylabelsize=fontsize,
        # title=L"\mathrm{Floquet~stability}",
        xticks=(xtickvalues, xticklabels),
        xtickalign=1.0, # ticks pointing inside
        yscale=log10,
        yticks=(ytickvalues, yticklabels),
        ytickalign=1.0  # ticks pointing inside
        )
    
    hm = heatmap!(axs, α, ρ, S, colormap=:deep)
    contour!(axs, α, ρ, S, levels=[Sthreshold], linewidth=2/10, color=:black, linestyle=:dash)

    sc = scatter!(axs, Z, color=log10.(r), markersize=10, colormap=[:orange, :green])
    
    lines!(axs, argmax_sr_cstfreq, linewidth=1.2, color=:white)
    
    # control corresponding to r = 1
    z1 = first(Z)
    scatter!(axs, z1, color=:purple, markersize=15, marker=:circle)
    # hlines!(axs, last(z1); xmax=first(z1), linewidth=1.5, color=:red, alpha=1.0)
    vlines!(axs, first(z1); ymax=0.88, linewidth=1.5, color=:purple, alpha=1.0)

    # control corresponding to r = +∞
    light_red = RGBAf(0.9, 0, 0.1, 0.9)
    scatter!(axs, z, color=light_red, markersize=15, marker=:x)
    hlines!(axs, last(z); xmax=first(z), linewidth=1.5, color=light_red, alpha=1.0)
    vlines!(axs, first(z); ymax=0.70, linewidth=1.5, color=light_red, alpha=1.0)

    foreach(hlinesfreq, hllabels, hllinestyles) do v, l, s
        hlines!(axs, v, color=:blue, linestyle=s, linewidth=1.5, label=l)
    end 
    axislegend(axs, L"ω/ω_0"; position=:cb, titlesize=18)


    Colorbar(fig[1, 2][:, 1], hm, label=L"|λ|_{\mathrm{max}}", labelsize=fontsize)
    Colorbar(fig[1, 2][:, 2], sc, label=L"\log_{10}(r)", labelsize=fontsize)
    xlims!(axs, xlimits)
    ylims!(axs, ylimits)
    fig
end


### save 
fig = fig_floquet2D()
save("./plots/floquet2D.png", fig)
