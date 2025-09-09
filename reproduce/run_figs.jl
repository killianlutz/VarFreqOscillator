relpath = "../scripts_figs/"
files = [
    "plot_amplification.jl",
    "plot_cstfreqslice.jl",
    "plot_floquet2D.jl",
    "plot_floquet3D.jl",
    "plot_trajectories.jl"
]

foreach(files) do file
    include(relpath*file)
end