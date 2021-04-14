# To do  

- [X] solve wave equation
- [X] add taper
- [X] save time snaps to binary file
- [X] get a better taper
- [X] generate synthetic data
- [X] make reverse propagation
- [X] image condition
- [X] subtract direct wave from seismogram before reverse propagation
- [X] create a wave equation module for solving and related stuff
- [X] create specific code for 1D/2D signal propagation w/ and w/o saving
- [X] read velocity field from binary file
- [X] read signal from binary file
- [X] dx != dz
- [X] different laplacian orders
- [X] madagascar Julia API: sadly, it yet doesn't exist
- [X] Signal1d/Signal2d -> Vector{Signal1d}(n_sources)
- [X] code cleanup
- [X] receptors position --> local offsets (source position, n_recept., spacing, array)
- [X] generate synthetic seismogram (one shot)
- [X] generate synthetic seismograms (various shots)
- [X] make complete migration routine
- [X] adjust image condition overhead
- [X] read some articles (*busquem conhecimento...*  - ET Bilu)
  - [X] 1983_baysal_rtm
  - [X] 1984_levin_rtm_principle
  - [X] 1989_mcmechan_acoustic_rtm_imaging_review **(a 10 for this one's didactic!)**
  - [X] 1996_loewenthal_huygens_principle_x_exploding_reflector
  - [X] 2018_sortberg_wave_eq_solution_via_ann
  - [X] 2018_pestana_rem_wave_equation_time_evolution
  - [X] 2020_neural_ode_exploration_and_implementation
- [X] read REM articles
  - [X] 2010_pestana
  - [X] 2019_pestana
- [X] program REM
- [X] learn ANN
- [X] learn autograd with ANN exaple (AutoGrad --> JAX)
- [X] apply autograd to a different goal (to gain fluency)
- [X] code 1D FDM in Python
- [X] make a data feeder for a 1D AWM PINN
- [X] correct REM: the issue was on bessel function times coeficients array (cJ)
- [X] REM low frequency noise: solved by adding another Chebyshev expansion term
- [X] gain confidence with Fourier transform frequency manipulation
- [X] code 2D laplacian as a frequency filter (pseudo-spectral method)
- [X] join propagation functions to an only routine with different options
  - [X] time fdm + space fdm
  - [X] time rem + space fdm
  - [X] time rem + space pseudo
- [X] cosmetics: parameters file like Reynam's
- [ ] train our neural network


# BFP
Compare traditional modeling (maybe rtm) with ANN model approximation


- [ ] train neural network with 1 propagation cube
- [ ] study LSTM ANNs
- [ ] create velocity fields database
- [ ] create useful propagation snaps database
- [ ] create WE aproximator using ANNs
- [ ] U-Net for inversion
- [ ] U-Net for fault identification
