## Purpose: reproducibility
This code, written in [`Julia`](https://julialang.org/), is provided to reproduce the figures shown in [this]() paper.

## Getting started 
Create a directory, open the terminal and run
```
git clone https://github.com/killianlutz/VarFreqOscillator.git
```

## Reproduce all figures
There are three steps:

- Open a terminal and change directory to the project VarfreqOscillator.

- Perform the necessary calculations, stored in the `./sims` directory, by running
```
julia --project=. reproduce/run_sims.jl
```

- Generate the figures, stored in the `./plots` directory, by running
```
julia --project=. reproduce/run_figs.jl
```

## Files associated to figures
Change directory to `./scripts_figs/`. There are five files, each of which is associated to one of the 5 figures of the paper:
- Fig. 1: `plot_trajectories.jl`
- Fig. 2: `plot_amplification.jl`
- Fig. 3: `plot_floquet2D.jl`
- Fig. 4: `plot_cstfreqslice.jl`
- Fig. 5: `plot_floquet3D.jl`

## Troubleshooting
Feel free to contact me by email: killian.lutz[at]inria.fr