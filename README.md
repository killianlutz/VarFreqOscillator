## Purpose: reproducibility
This code, written in [Julia](https://julialang.org/), is provided to reproduce the figures in [this]() paper.

## Getting started 
Open the terminal and run
```
git clone https://github.com/killianlutz/VarFreqOscillator.git
```

## Reproduce all figures
Assuming Julia version `1.11.7` is correctly [installed](https://docs.julialang.org/en/v1/manual/installation/), there are three steps:

- Open a terminal and change directory to this project. Your directory now ends with `./VarFreqOscillator`.

- Download and install the dependencies missing on your machine by running
```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

- Perform the necessary calculations, stored in the `./sims` directory, by running
```
julia --project=. reproduce/run_sims.jl
```

- Generate the figures, stored in the `./plots` directory, by running
```
julia --project=. reproduce/run_figs.jl
```

## Files associated to figures
In the directory `./scripts_figs/`, each one of the 5 files corresponds to one of the 5 figures in the paper:
- Fig. 1: `plot_trajectories.jl`
- Fig. 2: `plot_step_examples.jl`
- Fig. 3: `plot_floquet2D.jl`
- Fig. 4: `plot_cstfreqslice.jl`
- Fig. 5: `plot_floquet3D.jl`

## Troubleshooting
Feel free to reach out to me: [Killian Lutz](https://killianlutz.github.io/).