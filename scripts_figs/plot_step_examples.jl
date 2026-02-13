using JLD2
using CairoMakie
CairoMakie.activate!()

function fig_step_examples(; fig=Figure())
    @load "./sims/step_parameters.jld2" Ω A ts vs Ys ϕs ω0
    
    P1, P2 = map(x -> 2π/(x*ω0), Ω[1:2])
    Pmax = max(P1, P2)

    fontsize = 18
    linewidth = 1.5

    cmaps = [:viridis, :RdBu, :RdBu]
    colors = map(cmaps) do cmap; cgrad(cmap, length(ϕs); categorical=true); end
    linestyles = [i <= length(ϕs) ÷ 2 ? :solid : :solid for i in eachindex(ϕs)]

    colorrange = (minimum(first(vs)), 1.0)
    cmap = cgrad([:purple, :orange], 2)
    p = (; markersize=6, colormap=cmap, colorrange, marker=:cross)


    foreach(eachindex(Ω, A)) do i
        t = ts[i]
        v = vs[i]
        Y = Ys[i]

        ### phase space axis
        if i > 1
            xticks = (-1:1, [L"%$(x)" for x in -1:1])
            yticks = (-40:40:40, [L"%$(y)" for y in -40:40:40])

            axs = Axis(fig[i, 1]; 
                xlabel= i == 3 ? L"x\;\, \mathrm{[m]}" : "", 
                xlabelsize=fontsize,
                xticks,
                xtickalign=1.0, # ticks pointing inside
                ylabel=L"\dot{x}\;\, \mathrm{[m/s]}",
                yticks,
                ytickalign=1.0, # ticks pointing inside
                ylabelsize=fontsize, 
                # aspect=AxisAspect(1)
                )
            xlims!(-1.5, 1.5)
            ylims!(-55.0, 55.0)

        else
            xticks = (-7:7:7, [L"%$(x)" for x in -7:7:7])
            yticks = (-150:150:150, [L"%$(y)" for y in -150:150:150])

            axs = Axis(fig[i, 1]; 
                xlabel= i == 3 ? L"x\;\, \mathrm{[m]}" : "", 
                xlabelsize=fontsize,
                xticks,
                xtickalign=1.0, # ticks pointing inside
                ylabel=L"\dot{x}\;\, \mathrm{[m/s]}",
                ylabelsize=fontsize, 
                yticks,
                ytickalign=1.0, # ticks pointing inside
                # aspect=AxisAspect(1)
                )
        end

        # each trajectories and their starting point
        for (j, y) in enumerate(eachslice(Y; dims=2))
            Py = Point2.(eachcol(y))
            lines!(axs, Py; color=colors[i][j], linestyle=linestyles[j], linewidth)
            if i > 1
                # too cluttered if i == 1
                scatter!(axs, first(Py), color=colors[i][j], markersize=5)
            end
        end

        # same xlimits and common xticks for clarity
        if i == 2
            xticks = (ω0*[0, P2], [L"0", L"ω_0P_{%$(i)}"])
        else
            xticks = (ω0*[0, P1], [L"0", L"ω_0P_{%$(i)}"])
        end

        ### control axis
        axs = Axis(fig[i, 2]; 
            xlabel=i == 3 ? L"ω_0 t" : "", 
            xlabelsize=0.85fontsize,
            ylabel=L"S_{%$(i)}/u_{\max}", 
            ylabelsize=0.85fontsize,
            # aspect=AxisAspect(1),
            xticks,
            yticks=(0:1, [L"0", L"1"]),
            titlesize=15
            )
        xlims!(0.0, ω0*Pmax)

        # step function look
        φ = ω0*t
        stairs!(axs, φ, v; color=:black, alpha=0.2)
        scatter!(axs, φ, v; color=v, p...)
        hidespines!(axs)

        nothing
    end

    # colgap!(fig.layout, -200)
    # rowgap!(fig.layout, 10)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colsize!(fig.layout, 2, Aspect(1, 1.0))
    resize_to_layout!(fig)

    fig
end

### save 
fig = fig_step_examples()
save("./plots/step_examples.png", fig)
