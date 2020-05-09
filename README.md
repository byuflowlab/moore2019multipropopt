# Takeoff and Performance Tradeoffs of Retrofit Distributed Electric Propulsion for Urban Transport

_While vertical takeoff and landing aircraft have shown promise for urban air transport, distributed electric propulsion on existing aircraft may offer immediately implementable alter- natives. Distributed electric propulsion could potentially decrease takeoff distances enough to enable thousands of potential inter-city runways. This conceptual study explores the effects of a retrofit of open-bladed electric propulsion units. To model and explore the design space we use blade element momentum method, vortex lattice method, linear-beam finite element analysis, classical laminate theory, composite failure, empirically-based blade noise modeling, motor and motor-controller mass models, and gradient-based optimization. With liftoff time of seconds and the safe total field length for this aircraft type undefined, we focused on the minimum conceptual takeoff distance. We found that 16 propellers could reduce the takeoff distance by over 50% compared to the optimal 2 propeller case. This resulted in a conceptual minimum takeoff distance of 20.5 meters to clear a 50 ft (15.24 m) obstacle. We also found that when decreasing the allowable noise by approximately 10 dBa, the 8 propeller case performed the best with a 43% reduction in takeoff distance compared to the optimal 2 propeller case. This resulted in a noise-restricted conceptual minimum takeoff distance of 95 meters._

I don't claim this "glue code" as the cleanest or most finalized, however I hope that by making it open source, it will be of some value.  The majority of the sub-modules referenced in the paper are however of general high quality and maintained actively.  If I were a new researcher intending to build on this work, I would use the up-to-date sub modules and build my own glue code, referencing back to this occasionally as a supplement to better understand the theory and methods used in the published papers.

Contains data and source for paper and code.  The required Julia-0.6.0 packages are located in the deps folder in their state at the time of generation.  For a recreation of the environment, they would be placed in the .julia/v0.6 folder after julia-0.6.0 is installed.  SNOPT source code is NOT included.  At the time of writing, SNOPT was the only code required paid source.  Since the repositories of many of dependancies may have been massively updated, please refer to the official source for official requirements, etc.  Obviously not all of the required julia packages are included, specifically those which are widely used or are secondary or tertiary to the main ones used.  For those you would need to ensure you have installed the julia-0.6.0 compatible versions (June 11, 2018 is when most of mine were installed).

The license for my code does not apply for the dependent repositories.  They are simply included in their state to simplify environment recreation.

This repository was not intended to be used by anyone unfamiliar with Julia, MDAO, or the physics modeled in the study.  Therefore documentation is limited, however, many more details can be found in the thesis version of this work: https://scholarsarchive.byu.edu/etd/7537/  

Many engineering assumptions were made in this comparative study for the purposes of numerical stability and scope manageability, such as assuming airfoil rotational stall correction to be relatively constant for families of similar propeller designs and centrifugal forces to be second order to the intended comparison.

Within the src folder:
- airfoils: airfoil data and data corrections specific to the design case
- Optimization_Cases_Cruise: a .jl driver file for each cruise case
- Optimization_Cases_Noise: a .jl driver file for each noise case
- OptSummaries: Output data and processing scripts for figure generation, etc
- Remaining files are the optimization objective, sub files, an example driver, and several comparison/validation cases.
- Many of the .jld data files have been compressed to save space
