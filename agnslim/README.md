# agnslim, AGN super-Eddington accretion model

## Description

A broadband spectral model for a super-Eddington black hole accretion
disc developed by Kubota & Done (2019; KD19, MNRAS, 489, 524)}. This
is based on the slim disc emissivity (Abramowicz et al., 1988, ApJ,
332, 646; {Watarai et al., 2000, PASJ, 52, 133; Sadowki, 2011,
arXiv1108.0396), where radial advection keeps the surface luminosity
at the local Eddington limit, resulting in `L(r)` propto `1/r^2` rather
than the `1/r^3` expected from the Novikov-Thorne (standard,
sub-Eddington) disc emissivity. This is the only major change from the
sub-Eddington agnsed model (Kubota & Done 2018; KD18, MNRAS, 480,
1247), an updated version of optxagnf (Done et al. 2012, MNRAS, 420,
1848). The flow is radially stratified, with an outer standard disc
(from `R_out` to `R_warm`), an inner hot Comptonising region (`R_in` to
`R_hot`) and an intermediate warm Comptonising region to produce the
soft X-ray excess (`R_warm` to `R_hot`). A minor difference from agnsed is
that the disc is assumed to extend untruncated down to the inner
radius of the flow, `R_in`. This can be below the innermost stable
circular orbit as pressure forces are important. By default, the code
calculates its own expected value of `R_in` given the mass accretion
rate. However, we also allow this to be a free parameter e.g. for use
for the extreme super Eddington mass accretion rates probably truncate
at some radius from strong wind mass loss. Another minor difference
from agnsed is that we do not calculate the reprocessed emission as
the geometry of the inner disc is very uncertain but it probably
shields the outer flow.

The model calculates some useful quantities, such as the radius at
which the flux first goes above the local Eddington limit, and the
inner radius of the flow. These are not normally displayed but can be
seen by inputting the command chatter 20, and getting the model to
recalculate the fit e.g. by changing the normalisation to 1.0001. Set
this back to the default of chatter 10 to suppress all this
information if further fits are required.

## Parameters for agnslim:

Par  | Name           | Description
---  | ----           | -------------
par1 | `mass`         | black hole mass in solar masses
par2 | `dist`         | comoving (proper) distance in Mpc
par3 | `logmdot`      | `mdot = Mdot/Mdot_Edd` where `eta Mdot_Edd c^2 = L_Edd`
par4 | `astar`        | dimensionless black hole spin
par5 | `cosi`         | cosine of the inclination angle i for the warm Comptonising component and the outer disc.
par6 | `kTe_hot`      | electron temperature for the hot Comptonisation component in keV. If this parameter is negative then only the hot Comptonisation component is used.
par7 | `kTe_warm`     | electron temperature for the warm Comptonisation component in keV. If this parameter is negative then only the warm Comptonisation component is used.
par8 | `Gamma_hot`    | the spectral index of the hot Comptonisation component.
par9 | `Gamma_warm`   | the spectral index of the warm Comptonisation component. If this parameter is negative then only the outer disc component is used.
par10 | `R_hot`       | outer radius of the hot Comptonisation component in Rg
par11 | `R_warm`      | outer radius of the warm Comptonisation component in Rg
par12 | `logrout`     | log of the outer radius of the disc in units of Rg. If this parameter is negative, the code will use the self gravity radius as calculated from Laor & Netzer (1989, MNRAS, 238, 897L).
par13 | `R_in`        | the inner radius of the disc in Rg. If this parameter is -1 (the default), the model will use the radius calculated from KD19. This must be greater than `R_hot` for `mdot` greater than 6 and greater than `R_isco` for `mdot` less than 6.
par14 | `redshift`    |
par   | `norm`        | this must be fixed to 1.

