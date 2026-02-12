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
 
    u_max = maximum(v)
    η = 1.0 # by default in all simulations
    Γ = sqrt(4u_max/η^2 - 1)
    ω0 = η*Γ/2
    ω0T = round(Int, ω0*t[end])

    r = round(Int, r_target)
    u = v/maximum(v)
    colorrange = (minimum(u), 1.0)
    cmap = cgrad([:purple, :orange], length(unique(u)))
    p = (; markersize=6, color=u, colormap=cmap, colorrange, marker=:cross)

    #### axis for first trajectory
    axtraj = Axis(fig[1, 1], 
        xlabel=L"x\;\, \mathrm{[m]}", 
        xlabelsize=fontsize,
        xticks=([-1, 0, 2], [L"%$(i)" for i in [-1, 0, 2]]),
        ylabel=L"\dot{x}\;\, \mathrm{[m/s]}",
        ylabelsize=fontsize, 
        yticks=([0, 10], [L"%$(i)" for i in [0, 10]]),
        xtickalign=1.0, # ticks pointing inside
        ytickalign=1.0, # ticks pointing inside
        aspect=AxisAspect(1)
        )
    ylims!(axtraj, -1.0, 15.0)

    ##### commutation line
    θ_slope = find_slope(x, dx, v)
    lines!(axtraj, 0.0..r_target, z -> θ_slope*z, color=first(cmap), linestyle=:dot, linewidth=2.0)    

    ##### fill surface to indicate angle θ and add label
    xs = range(0.0, 15, 2)
    ylower = xs/θ_slope 
    yupper = 0*xs
    band!(axtraj, xs, ylower, yupper; direction=:y, color=first(cmap), alpha=0.1)

    position = Point2(1/2, 10.0)
    text = L"θ^*"
    offset = (-8, -35)
    text!(axtraj, position; text, fontsize, offset, color=:purple)

    ##### starting at (-1, 0) instead of (-3/r, 0).
    xx = r_target/3 * x
    dxx = r_target/3 * dx
    scatter!(axtraj, xx, dxx; p...)

    ##### end points of the trajectory
    positions = [Point2(-1.0, 0.0), Point2(r_target, 0.0)]
    text = [L"A", L"B"]
    align = [(:left, :bottom), (:right, :bottom)]
    offset = [(5, 0), (-7, 0)]
    scatter!(axtraj, positions, color=:black, marker=:circle, markersize=12)
    text!(axtraj, positions; text, align, fontsize=15, offset)

    ##### first control
    axctrl1 = Axis(fig[1, 2], 
        xlabel=L"t/T", 
        xlabelsize=0.85fontsize,
        ylabel=L"u/u_{\max}", 
        ylabelsize=0.85fontsize,
        title=L"r = %$(r), \quad \omega_0T^* ≃ %$(ω0T)", 
        aspect=AxisAspect(1),
        xticks=(0:1, [L"0", L"1"]),
        yticks=(0:1, [L"0", L"1"]),
        titlesize=15
        )
    # step function look
    stairs!(axctrl1, t/T, u; color=:black, alpha=0.2)
    scatter!(axctrl1, t/T, u; p...)
    hidespines!(axctrl1)
    # subdivide into periods
    k = 1
    vlines!(axctrl1, [i/k for i in 0:k], color=:black, linestyle=:dash, alpha=1.0, linewidth=1.5)

    ###############
    ############### SECOND EXAMPLE
    ###############
    @load "./sims/trajectory_r15.jld2" r_target t x dx v T  

    u_max = maximum(v)
    η = 1.0 # by default in all simulations
    Γ = sqrt(4u_max/η^2 - 1)
    ω0 = η*Γ/2
    ω0T = round(Int, ω0*t[end])

    r = round(Int, r_target)
    u = v/maximum(v)
    colorrange = (minimum(u), 1.0)
    cmap = cgrad(:RdYlGn_9, length(unique(u)); rev=true, categorical=true)
    p = (; markersize=6, color=u, colormap=cmap, colorrange, marker=:cross)
    
    #### axis for second trajectory
    axtraj = Axis(fig[2, 1], 
        xlabel=L"x\;\, \mathrm{[m]}", 
        xlabelsize=fontsize,
        xticks=([-8, 0, 15], [L"%$(i)" for i in [-8, 0, 15]]),
        ylabel=L"\dot{x}\;\, \mathrm{[m/s]}",
        ylabelsize=fontsize, 
        yticks=([0, 60], [L"%$(i)" for i in [0, 60]]),
        xtickalign=1.0, # ticks pointing inside
        ytickalign=1.0, # ticks pointing inside
        aspect=AxisAspect(1)
        )
    ylims!(axtraj, -40, 90)

    ##### commutation line
    θ_slope = find_slope(x, dx, v)
    lines!(axtraj, -7.0..r_target, z -> θ_slope*z, color=first(cmap), linestyle=:dot, linewidth=2.0)    

    ##### fill surface to indicate angle θ and add label
    # xs = range(-40.0, 0.0, 2)
    # ylower = 0*xs
    # yupper = xs/θ_slope
    # band!(axtraj, xs, ylower, yupper; direction=:y, color=first(cmap), alpha=0.1)

    xs = range(0.0, 90.0, 2)
    ylower = xs/θ_slope 
    yupper = 0*xs
    band!(axtraj, xs, ylower, yupper; direction=:y, color=first(cmap), alpha=0.1)

    position = Point2(2.5, 50.0)
    text = L"θ^*"
    offset = (-5, -25)
    text!(axtraj, position; text, fontsize, offset, color=first(cmap))

    ##### starting at (-1, 0) instead of (-3/r, 0).
    xx = r_target/3 * x
    dxx = r_target/3 * dx
    scatter!(axtraj, xx, dxx; p...)

    ##### end points of the trajectory
    positions = [Point2(-1.0, 0.0), Point2(r_target, 0.0)]
    text = [L"A", L"B"]
    align = [(:right, :bottom), (:right, :bottom)]
    offset = [(-5, 0), (-7, 0)]
    scatter!(axtraj, positions, color=:black, marker=:circle, markersize=12)
    text!(axtraj, positions; text, align, fontsize=15, offset)

    ##### second control
    axctrl2 = Axis(fig[2, 2], 
        xlabel=L"t/T", 
        xlabelsize=0.85fontsize,
        ylabel=L"u/u_{\max}", 
        ylabelsize=0.85fontsize,
        title=L"r = %$(r), \quad \omega_0T^* ≃ %$(ω0T)", 
        aspect=AxisAspect(1),
        xticks=(0:1, [L"0", L"1"]),
        yticks=(0:1, [L"0", L"1"]),
        titlesize=15
        )
    # step function look
    stairs!(axctrl2, t/T, u; color=:black, alpha=0.2)
    scatter!(axctrl2, t/T, u; p...)
    hidespines!(axctrl2)
    # subdivide into periods
    k = 3
    vlines!(axctrl2, [i/k for i in 0:k], color=:black, linestyle=:dash, alpha=1.0, linewidth=1.5)
    
    colgap!(fig.layout, -100)

    return fig
end


### save 
fig = fig_trajectories()
save("./plots/trajectories.png", fig)