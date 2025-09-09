using JLD2
using CairoMakie
CairoMakie.activate!()

function find_slope(x, dx, v)
    τ = -Inf64
    mean_extrema = 0.5*(maximum(v) + minimum(v))
    second_switch = false

    for i in eachindex(v)
        if v[i] < mean_extrema
            second_switch = true
        end

        if second_switch == true
            if v[i] > mean_extrema
                τ_before = -x[i-1]/dx[i-1]
                τ_after = -x[i]/dx[i]
                τ = 0.5*(τ_before + τ_after)
                break
            end
        end
    end
    slope = -1/τ

    return slope
end

function fig_trajectories(; fig=Figure())
    fontsize = 20
    
    ############### FIRST EXAMPLE
    @load "./sims/trajectory_r2.jld2" r_target t x dx v T   
 
    T = round(t[end], digits=2)
    r = round(Int, r_target)
    u = v/maximum(v)
    colorrange = (minimum(u), 1.0)
    cmap = cgrad([:purple, :orange], length(unique(u)))
    p = (; markersize=6, color=u, colormap=cmap, colorrange, marker=:cross)
    
    #### common axis for trajectories
    axtraj = Axis(fig[1, 1], 
        xlabel=L"x", 
        xlabelsize=fontsize,
        xticks=(-1:3, [L"%$(i)" for i in -1:3]),
        ylabel=L"\dot{x}",
        ylabelsize=fontsize, 
        yticks=(0:10:20, [L"%$(i)" for i in 0:10:20]),
        aspect=AxisAspect(1)
        )
    
    ##### commutation line
    θ_slope = find_slope(x, dx, v)
    lines!(axtraj, 0.0..3.0, z -> θ_slope*z, color=first(cmap), linestyle=:dot, linewidth=2.0)    

    ##### first trajectory
    scatter!(axtraj, x, dx; p...)

    ##### first control
    axctrl1 = Axis(fig[1, 2][1, 1:2], 
        xlabel=L"t/T", 
        xlabelsize=0.85fontsize,
        ylabel=L"u/\bar{u}", 
        ylabelsize=0.85fontsize,
        title=L"r = %$(r), \, T^* ≃ %$(T)", 
        aspect=AxisAspect(1),
        xticks=(0:1, [L"0", L"1"]),
        yticks=(0:1, [L"0", L"1"])
        )
    # step function look
    stairs!(axctrl1, t/T, u; color=:black, alpha=0.2)
    scatter!(axctrl1, t/T, u; p...)
    hidespines!(axctrl1)
    # subdivide into periods
    k = 1
    vlines!(axctrl1, [i/k for i in 0:k], color=:black, linestyle=:dash, alpha=1.0, linewidth=1.5)
    
    ############### SECOND EXAMPLE
    @load "./sims/trajectory_r15.jld2" r_target t x dx v T  

    T = round(t[end], digits=2)
    r = round(Int, r_target)
    u = v/maximum(v)
    colorrange = (minimum(u), 1.0)
    cmap = cgrad(:RdYlGn_9, length(unique(u)); rev=true, categorical=true)
    p = (; markersize=6, color=u, colormap=cmap, colorrange, marker=:cross)
    
    ##### commutation line
    θ_slope = find_slope(x, dx, v)
    lines!(axtraj, -1.5..3.0, z -> θ_slope*z, color=first(cmap), linestyle=:dot, linewidth=2.0)    

    ##### second trajectory
    scatter!(axtraj, x, dx; p...)

    ##### second control
    axctrl2 = Axis(fig[1, 2][2, 1:2], 
        xlabel=L"t/T", 
        xlabelsize=0.85fontsize,
        ylabel=L"u/\bar{u}", 
        ylabelsize=0.85fontsize, 
        title=L"r = %$(r), \, T^* ≃ %$(T)", 
        aspect=AxisAspect(1),
        xticks=(0:1, [L"0", L"1"]),
        yticks=(0:1, [L"0", L"1"])
        )
    # step function look
    stairs!(axctrl2, t/T, u; color=:black, alpha=0.2)
    scatter!(axctrl2, t/T, u; p...)
    hidespines!(axctrl2)
    # subdivide into periods
    k = 3
    vlines!(axctrl2, [i/k for i in 0:k], color=:black, linestyle=:dash, alpha=1.0, linewidth=1.5)
    
    return fig
end

### save 
fig = fig_trajectories()
save("./plots/trajectories.png", fig)