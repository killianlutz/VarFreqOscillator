relpath = "../scripts_sims/"
files = [
    "optstepfun_amplification.jl",
    "optstepfun_natfreq.jl",
    "floquet2d.jl",
    "floquet2d_slices.jl",
    "specrad_argmax.jl",
    "floquet3d.jl",
    "solution_amplification.jl",
    "trajectories.jl",
    "step_examples.jl"
]

foreach(files) do file
    include(relpath*file)
end