# Technical Document

## Interacting Order Book

Used https://github.com/DerickDiana/InteractingLOBs.jl to generate simulations of a lit order book

Then took the Epps modelling in https://github.com/CHNPAT005/PCEPTG-MSC and used it to model the epps effect of the interacting order books

* Not really sure how to create an envelope

For calibration I adapted the NMTA in https://github.com/matthewdicks98/MDTG-MALABM/tree/Calibrated-ABM which uses method of moments to calibrate an optimal solution

* Issues I need to handle are the upper bounds

For L and M it seems that L/M needs to be < 0.67 otherwise it causes out of bounds error
I added a check for this in NMTA
