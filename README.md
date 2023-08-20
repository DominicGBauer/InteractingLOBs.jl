# Simulating Anomalous Diffusion Models and Calculating Price Impacts of The Lit And Latent Order-Book in Julia

## Authors
* Derick Diana
* Tim Gebbie

## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstmann for their assistance and advice with the formulation of a numerical solution of the SPDE necessary for our implementation. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Derick Diana would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.

## Capabilities

In its current form, the code can simulate a latent order book which obeys the equations shown in our paper. This includes a form which has anomalous diffusion. In addition, the code can simulate n different order books. It then allows the user to specify the source function and random kick function at any time via standard functions which are passed all information the user could possibly want to determine how the source and randomness functions should depend on the rest of the system and on time. In particular this allows the order books to be coupled by specifying the source or randomness terms of the any one order book to depend on the information contained in one or many others. Finally the code has the capibility to produce lucid plots showing the evolving PDEs over time, snap shots at various times and also has the ability to produce price impact functions via a standard outer function.

## Quick Start

* In terminal run `make`
*

1. Run `julia`
2. Run `]` and `activate .`
3. Run `include("src/Epps/epps.jl")`
4. Run `include("src/Epps/work.jl")`
5. Run `include("src/Epps/price_paths.jl")`
6. Run `include("src/Epps/price_surface.jl")`
7. Run `include("src/Epps/price_impact.jl")`
