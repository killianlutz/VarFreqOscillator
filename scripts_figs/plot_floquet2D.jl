using JLD2
using CairoMakie
CairoMakie.activate!()

function fig_floquet2D(; fig=Figure())
    @load "./sims/floquet2d.jld2" Γ γ nα nρ α ρ S S01
    @load "./sims/specrad_argmax.jld2" argmax_sr_cstfreq
    @load "./sims/step_examples.jld2" A Ω
    @load "./sims/optstepfun_amplification.jld2" Γ γ z Z r
    @load "./sims/floquet2d_slices_r_finite.jld2" z8

    ρ1 = last(first(Z))
    ρ1_round = round(ρ1, digits=3)
    ρinf = last(z)
    ρinf_round = round(ρinf, digits=3)
    ρ8 = last(z8)
    ρ8_round = round(ρ8, digits=3)

    fontsize = 20

    xlimits = (0.0, 1.0)
    ylimits = (1e-1, 3e0)
    Sthreshold = 1.0

    xtickvalues = [0, 1/4, 1/2, 3/4, 1]
    xticklabels = [L"0", L"1/4", L"1/2", L"3/4", L"1"]
    ytickvalues = [1e-1, 1]
    yticklabels = [L"10^{%$(n)}" for n in -1:0]

    hlinesfreq = reverse([ρ8, ρinf, ρ1])
    hllabels = reverse([L"%$(ρ8_round)", L"%$(ρinf_round)", L"%$(ρ1_round)"])
    hllinestyles = [:solid, :solid, :solid]
    hlcolors = reverse([:blue, :red, :purple])

    axs = Axis(
        fig[3, 1], 
        # aspect=AxisAspect(1), 
        xlabel=L"α",
        xlabelsize=1.3fontsize,
        ylabel=L"ω/ω_0",
        ylabelsize=1.5fontsize,
        # title=L"\mathrm{Floquet~stability}",
        xticks=(xtickvalues, xticklabels),
        xticklabelsize=fontsize,
        xtickalign=1.0, # ticks pointing inside
        yscale=log10,
        yticks=(ytickvalues, yticklabels),
        yticklabelsize=fontsize,
        ytickalign=1.0  # ticks pointing inside
        )
    
    hm = heatmap!(axs, α, ρ, S, colormap=:deep)
    contour!(axs, α, ρ, S, levels=[Sthreshold], linewidth=2/10, color=:black, linestyle=:dash)

    sc = scatter!(axs, Z, color=log10.(r), markersize=12, colormap=[:orange, :green])
    
    lines!(axs, argmax_sr_cstfreq, linewidth=2, color=:white)
    
    # control corresponding to r = 1
    z1 = first(Z)
    text = L"r = 1"
    align = (:right, :bottom)
    offset = (0, 10)
    scatter!(axs, z1, color=:purple, markersize=20, marker=:circle)
    text!(axs, z1; text, align, offset, fontsize=0.9fontsize, color=:purple)

    # control corresponding to r = 8
    text = L"r = 8"
    align = (:right, :bottom)
    offset = (0, 10)
    scatter!(axs, z8, color=:blue, markersize=25, marker=:star4)
    text!(axs, z8; text, align, offset, fontsize=0.9fontsize, color=:blue)

    # control corresponding to r = +∞
    text = L"r = ∞"
    align = (:right, :bottom)
    offset = (0, 10)
    scatter!(axs, z, color=:red, markersize=22, marker=:x)
    text!(axs, z; text, align, offset, fontsize=0.9fontsize, color=:red)

    foreach(hlinesfreq, hllabels, hllinestyles, hlcolors) do v, l, s, c
        hlines!(axs, v, color=c, linestyle=s, linewidth=1.5, label=l)
    end 
    axislegend(axs, L"ω/ω_0"; position=:cb, titlesize=22, labelsize=17)

    ### location of step function examples
    steps = [Point2(x, y) for (x, y) in zip(A, Ω)]
    text = [L"S_1", L"S_2", L"S_3"]
    align = [(:left, :top), (:right, :center), (:right, :center)]
    offset = [(3, -3), (-8, 0), (-8, 0)]
    scatter!(axs, steps; marker=:rect, markersize=16, color=:black)
    text!(axs, steps; text, align, offset, fontsize)

    ### colorbars
    Colorbar(fig[1, 1], hm, label=L"|λ|_{\mathrm{max}}", labelsize=fontsize, vertical=false, ticks=([1, 10, 30], [L"%$(x)" for x in [1, 10, 30]]), labelpadding=-10)
    Colorbar(fig[2, 1], sc, label=L"\log_{10}(r)", labelsize=fontsize, vertical=false, ticks=(0:3, [L"%$(x)" for x in 0:3]), labelpadding=-10)
    xlims!(axs, xlimits)
    ylims!(axs, ylimits)

    rowsize!(fig.layout, 3, Aspect(1, 1.0))    
    resize_to_layout!(fig)
    fig
end

### save 
fig = fig_floquet2D()
save("./plots/floquet2D.png", fig)
