using JLD2
using CairoMakie
CairoMakie.activate!()

function fig_cstfreqslice(; fig=Figure())
    fontsize = 20
    textsize = 21
    legendsize = 22

    axs = Axis(
        fig[1, 1], 
        aspect=AxisAspect(1), 
        xlabel=L"α",
        xlabelsize=fontsize,
        ylabel=L"|λ|_{\mathrm{max}}",
        ylabelsize=fontsize,
        # title=L"α ↦ |λ|_{\mathrm{max}}(ω, α)",
        # titlesize=fontsize,
        yticks=(1:2:7, [L"%$(k)" for k in 1:2:7]),
        xticks=([i/4 for i in 0:4], [L"0", L"1/4", L"1/2", L"3/4", L"1"])
        )
    hlines!(1.0, color=:black, alpha=0.5)
    
    ### finite amplification
    @load "./sims/floquet2d_slices_r_finite.jld2" k r_val α_opt α cst_freq_slice
    max_λ = cst_freq_slice[argmax(cst_freq_slice)]

    lines!(α, cst_freq_slice, linewidth=2, color=:black, label=L"ω^*_{r = 8}")
    vlines!(α_opt; ymax=0.95, linewidth=2.0, color=:blue)
    text!(α_opt, max_λ; text=L"α^*_{r = 8}", color=:blue, offset=(10, -12), fontsize=textsize)
    scatter!(α_opt, max_λ, color=:blue, markersize=12)

    ### infinite amplification
    @load "./sims/floquet2d_slices_r_infinity.jld2" k r_val α_opt α cst_freq_slice
    max_λ = cst_freq_slice[argmax(cst_freq_slice)]

    lines!(α, cst_freq_slice, linewidth=2, color=:black, linestyle=:dash, label=L"ω^*_{r = ∞}")
    vlines!(α_opt; ymax=0.44, linewidth=2.0, linestyle=:dash, color=:red)
    text!(α_opt, max_λ; text=L"α^*_{r = ∞}", color=:red, offset=(-13, 10), fontsize=textsize)
    scatter!(α_opt, max_λ, color=:red, markersize=12)
    
    axislegend(axs, position=:rt, labelsize=legendsize)
    hideydecorations!(axs, label=false, ticklabels=false, ticks=false)
    return fig
end

### save 
fig = fig_cstfreqslice()
save("./plots/cstfreqslice.png", fig)