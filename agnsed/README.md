
A model for the spectral energy distribution (SED) of an AGN developed
by Kubota & Done 2018 (KD18;
https://ui.adsabs.harvard.edu//#abs/2018MNRAS.480.1247K/abstract).
Following Done et al. 2012
(https://ui.adsabs.harvard.edu//#abs/2012MNRAS.420.1848D/abstract),
the SED model has three characteristic regions: the outer standard
disc region; the warm Comptonising region; and the inner hot
Comptonising region.

For the warm Comptonising region, this model adopts the passive disc
scenario tested by Petrucci et al. 2018
(https://ui.adsabs.harvard.edu//#abs/2018A&A...611A..59P/abstract). Here,
the flow is assumed to be completely radially stratified, emitting as
a standard disc blackbody from Rout to Rwarm, as warm Comptonisation
from Rwarm to Rhot and then makes a transition to the hard X-ray
emitting hot Comptonisation component from Rhot to RISCO. The warm
Comptonisation component is optically thick, so is associated with
material in the disc. Nonetheless, the energy does not thermalise to
even a modified blackbody, perhaps indicating that significant
dissipation takes place within the vertical structure of the disc,
rather than being predominantly released in the midplane.

At a radius below Rhot, the energy is emitted in the hot
Comptonisation component. This has much lower optical depth, so it is
not the disc itself. In the model, the albedo is fixed at a = 0.3, and
the seed photon temperature for the hot Comptonisation component is
calculated internally. In contrast to optxagnf, this model does not
take the color temperature correction into account.

There are two versions of the model, agnsed and qsosed. agnsed is the
full model, while qsosed is a simplified version of agnsed made by
fixing some parameters at their typical values and by including
reprocessing.  For qsosed the agnsed parameters are fixed at kTe\_hot =
100 keV, kTe\_warm = 0.2 keV, Gamma\_warm = 2.5, R\_warm = 2R\_hot, rout = rsg
and Htmax = 100. Also, Gamma\_hot is calculated via eq.(6) in KD18 and R\_hot
is calculated to satisfy Ldiss\_hot = 0.02LEdd.


Parameters for agnsed.

- par1     `mass`: black hole mass in solar masses
- par2     `dist`: comoving (proper) distance in Mpc
- par3     `logmdot`: mdot = Mdot/Mdot\_Edd where eta Mdot\_Edd c^2 = L\_Edd
- par4     `astar`: dimensionless black hole spin
- par5     `cosi`: cosine of the inclination angle i for the warm Comptonising component and the outer disc.
- par6     `kTe\_hot`: electron temperature for the hot Comptonisation component in keV. If this parameter is negative then only the
  	   hot `Comptonisation component is used.
- par7     `kTe\_warm`: electron temperature for the warm Comptonisation component in keV. If this parameter is negative then only the warm Comptonisation component is used.
- par8     `Gamma\_hot`: the spectral index of the hot Comptonisation component. If this parameter is negative, the code will use the value calculated via eq.(2) of KD18.
- par9     `Gamma\_warm`: the spectral index of the warm Comptonisation component. If this parameter is negative then only the outer disc component is used.
- par10    `R\_hot`: outer radius of the hot Comptonisation component in Rg
- par11    `R\_wwarm`: outer radius of the warm Comptonisation component in Rg
- par12    `logrout`: log of the outer radius of the disc in units of Rg. If this parameter is negative, the code will use the self gravity radius as calculated from Laor & Netzer 1989.
- par13    `Htmax`: the upper limit of the scaleheight for the hot Comptonisation component in Rg. If this parameter is smaller than parameter 10, the hot Comptonisation region is a sphere of radius Htmax by keeping Ldiss\_hot determined by R\_hot via eq.(2) of KD18.
- par14    `reprocess`: switching parameter for the reprocessing, 0 or 1. If this parameter is 0, reprocessing is not considered. If this parameter is 1, reprocessing is included.
- par15    `redshift`


Parameters for qsosed

- par1     `mass`: black hole mass in solar masses
- par2     `dist`: comoving (proper) distance in Mpc
- par3     `logmdot`: mdot = Mdot/Mdot\_Edd where eta Mdot\_Edd c^2 = L\_Edd
- par4     `astar`: dimensionless black hole spin
- par5     `cosi`: cosine of the inclination angle i for the warm Comptonising component and the outer disc.
- par6     `redshift`
