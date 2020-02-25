# ThComp: Thermally Comptonized continuum

## Description

ThComp is a replacement of the Nthcomp model (Zdziarski et al. 1996, MNRAS, 283, 193). It agrees much better than Nthcomp with actual Monte Carlo spectra from Comptonization, see Zdziarski et al. (2019, arXiv:1910.04535) for details. See Nied&#x017A;wiecki et al. (2019, MNRAS, 485, 2942) for analogous comparison with Nthcomp, showing substantial discrepancies. ThComp describes spectra from Comptonization by thermal electrons emitted by a spherical source with the sinusoidal-like spatial distribution of the seed photons (as in compST, Sunyaev & Titarchuk 1980, A&A, 86, 121). It is a convolution model, and thus it can Comptonize any seed photon distribution, either hard or soft, and it describes both upscattering and downscattering (see Zdziarski et al. 2019, arXiv:1910.04535 for examples). In the case of upscattering of some seed photons (e.g. blackbody or disc blackbody), it is a much better description of the continuum shape from thermal Comptonization than an exponentially cutoff power law, but has similar corresponding free parameters, the spectral index, `Gamma`, and the high-energy cutoff, parameterized by the electron temperature (`kT_e`). That cutoff is much sharper than an exponential. The model also provides correct description of Comptonized spectra at energies comparable to those of the seed photons. Note that the model has no normalization parameter since its normalization follows from that of the seed photons.

Please reference Zdziarski et al. (2019, arXiv:1910.04535) if you use it. 

## Parameters for thcomp:

Par  | Name           | Description
---  | ----           | -------------
par1 | `Gamma`        | >0: the low-energy power-law photon index; <0: the Thomson optical depth (given by the absolute value).
par2 | `kT_e`         | electron temperature (high energy rollover)
par3 | `f_{sc}`       | the scattering fraction, 0 <= f_{sc} <= 1; if 1 all of the seed photons will be Comptonized; if set to 0, only the original seed photons will be seen.
par4 | `redshift`     | 
