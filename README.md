[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JanLeppert/IBS/HEAD?filepath=%2Fnotebooks%2FTGGCibs.ipynb)

# IBS

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> IBS

It is authored by Jan Leppert.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

## Project discription

This project is supplementary material to the paper "Simulation of ideal basic separations in gas chromatography with negative temperature gradients" by Jan Leppert, Leonid M. Blumberg, Matthias Wüst and Peter Boeker.

This github project includes a mybinder notebook with a version of the simulation.