# reXcor — a model of the X-ray spectrum of AGNs that combines ionized reflection and a warm corona

## Description

<p>
reXcor is a new AGN X-ray (0.3–100 keV) spectral fitting model that self-consistently combines the effects of both ionized relativistic reflection and the emission from a warm corona. In this model, the accretion energy liberated in the inner disk is distributed between a warm corona, a lamppost X-ray source, and the accretion disk. The emission and ionized reflection spectrum from the inner 400 rg of the disk is computed, incorporating the effects of relativistic light-bending and blurring. The resulting spectra predict a variety of soft excess shapes and sizes that depend on the fraction of energy dissipated in the warm corona and lamppost. 
<p>
Full details of the reXcor models can be found in this paper: Xiang, X, Ballantyne, D.R., et al., 2022, MNRAS, in press (arXiv number to come).
  
## Parameters
  
The reXcor parameters are:
<table>  
  <tr><td>1.</td><td>fx —  the fraction of the accretion flux dissipated in the lamppost (the hot, power-law emitting corona). Range of values: 0.02 to 0.2 (in steps of 0.02).</td></tr>
  <tr><td>2.</td><td>Gamma — the photon-index of the power-law continuum. Range of values: 1.7 to 2.2 (in steps of 0.05)</td></tr>
  <tr><td>3.</td><td>hf — the fraction of the accretion flux dissipated in the warm corona. Range of values: 0.0 to 0.8 (in steps of 0.05)</td></tr>
  <tr><td>4.</td><td>tau — Thomson depth of the warm corona. Range of values: 10 to 30 (in steps of 2)</td></tr>
  <tr><td>5.</td><td>z — Redshift.</td></tr>
  <tr><td>6.</td><td>Normalization — determined by many quantities, including the distance to the source, the area of the disk emission region, the inclination angle of the disk, the black hole mass and spin, and the geometry of the X-ray source. We recommend using cflux to determine the normalization of the reXcor model.</td</tr>
</table>
<p>
  
## Table model files
  
<p>  
We provide 8 different reXcor grids that were computed for different values of the lamppost height (h), the black hole spin (a), and the AGN Eddington ratio (lambda). All models are computed assuming Solar abundances.
<table>
  <tr><td>reXcor_l001_a09_h5.fits</td><td> lambda =0.01, a=0.9, h=5 rg</td></tr>
  <tr><td>reXcor_l001_a09_h20.fits</td><td> lambda=0.01, a=0.9, h=20 rg</td></tr>
  <tr><td>reXcor_l001_a099_h5.fits</td><td> lambda=0.01, a=0.99, h=5 rg</td></tr>
  <tr><td>reXcor_l001_a099_h20.fits</td><td> lambda=0.01, a=0.99, h=20 rg</td></tr>
  <tr><td>reXcor_l01_a09_h5.fits</td><td> lambda=0.1, a=0.9, h=5 rg</td></tr>
  <tr><td>reXcor_l01_a09_h20.fits</td><td> lambda=0.1, a=0.9, h=20 rg</td></tr>
  <tr><td>reXcor_l01_a099_h5.fits</td><td> lambda=0.1, a=0.99, h=5 rg</td></tr>
  <tr><td>reXcor_l01_a099_h20.fits</td><td> lambda=0.1, a=0.99, h=20 rg</td></tr>
  </table>
  <p>
Questions about the use of these models can be send to david.ballantyne -at- physics.gatech.edu