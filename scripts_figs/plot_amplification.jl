using JLD2
using CairoMakie
CairoMakie.activate!()

function fig_amplification(; fig=Figure())
    legendsize = 22
    fontsize = 20

    @load "./sims/solution_amplification.jld2" minimizers p r limit_ratio critical_ratio
    
    T, k, θ = eachrow(minimizers)

    u_max = p.Vmax
    η = p.damping_rate
    Γ = sqrt(4u_max/η^2 - 1)
    ω0 = η*Γ/2
    T = ω0*T/15

    fig = Figure(size=(500, 500))
    axs = Axis(
        fig[1, 1], 
        aspect=AxisAspect(1), 
        xlabel=L"r",
        xlabelsize=fontsize,
        xscale=log10,
        yticks=([-1, 0, 1, 3, 5], [L"-1", L"0", L"1", L"3", L"5"]),
        xticks=([1e-4, 1e-2, 1, 1e2, 1e4], [L"10^{-4}", L"10^{-2}", L"10^{0}", L"10^{2}", L"10^{4}"])
        )
    scatter!(axs, r, T, label=L"ρω_0T^*" => (; markersize=13), color=:purple, markersize=6, marker=:circle)
    scatter!(axs, r, θ, label=L"θ^*" => (; markersize=13), color=:orange, markersize=6, marker=:cross)
    scatter!(axs, r, k, label=L"k^*" => (; markersize=13), color=:teal, markersize=6, marker=:diamond)
    # vlines!(axs, critical_ratio, color=:red, label=L"r_c", linestyle=:dot, linewidth=1.5)
    vlines!(axs, limit_ratio, color=:black, label=L"\bar{r}" => (; linewidth=2), linestyle=:dash, linewidth=1.0)
    axislegend(axs, position=:lt, orientation=:horizontal, labelsize=legendsize)
    
    # shade regions with counter-intuitive result
    x = r[(1 .<= r .<= limit_ratio) .* (k .>= 3)]
    band!(axs, x, -1, 4, color=:gray, alpha=0.1)

    return fig
end


### save 
fig = fig_amplification()
save("./plots/amplification.png", fig)