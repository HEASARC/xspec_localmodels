# Cluster Mass Mixing Models

NOTE: These models are now in the Heasoft release, as of Xspec version 12.11.

These models are for determining distributions of gravitating mass in
spherical, hydrostatic atmospheres.  A hot atmosphere is approximated
as a set of concentric spherical shells, each containing isothermal
gas.  The mixing models combine thermal spectra for the shells, with
weights determined by the gravitational potential for the model, to
produce projected spectra for annular regions centred on a cluster.
The models have much in common with the projct mixing model, but the
gas density distribution is determined by the gravitational potential
and the assumption of hydrostatic equilibrium.  Details of the models
are discussed in Nulsen, Powell & Vikhlinin (2010) and the clmass
model explained fully.  For the clmass model, the gravitating matter
density is assumed to be constant in each spherical shell.  The
monomass model is physically identical to clmass, but parametrized by
the differences in mass densities for adjacent shells to ensure that
the gravitating mass density is a non-increasing function of the
radius.  The nfwmass model treats an atmosphere as a nested set of
isothermal, spherical shells, but with the Navarro, Frenk & White form
for the gravitational potential.

Spectra should be extracted from concentric, circular annuli centered
on a cluster (elliptical annuli are rejected).  All spectra for one
annulus must belong to the same data group, while spectra for distinct
annuli must belong to separate data groups.  Spectra must be provided
for a complete set of annuli filling the range between the innermost
and outermost radius.  For each model, the inner radius of the
innermost annulus is specified as the first model parameter.

To handle cases where the observations do not cover the whole of a
cluster, X-ray emission from beyond the inner edge of the outermost
annulus may be modeled with an isothermal beta model.  This feature is
designed to deal with background X-ray emission from parts of a
cluster outside the region that has been observed.  It adds a
model-dependent element to the mass models.  When the beta model is
used, the pressure is assumed to be continuous between the two
outermost shells, but the gravitational potential is ignored for the
outermost shell.  If the beta model is disabled, the outermost shell
is treated like any other.


## Using the models

The models rely on the XFLT keywords in much the same way as projct.
Each spectrum must include, at least, the XFLT0001 keyword specifying
the outer radius of the corresponding annulus.  If present, XFLT0002
must equal XFLT0001 (annuli must be circular).  XFLT0003 is ignored.
If present, the following pairs of keywords (XFLT0004/5, XFLT0006/7,
etc) give ranges of angle that are summed and divided by 360 to
determine the fraction of the total annulus covered by a spectrum.
Note that the same effect can be achieved by specifying this fraction
in the AREASCAL keyword and leaving XFLT0004/5, etc, undefined (this
approach provides greater flexibility when there is more than one
spectrum in each data group).

It is essential for the models to link the temperature of the thermal
model for each shell to the corresponding temperature parameter of the
model.  All unused shell parameters must be frozen.  The norms of the
thermal model must be tied (equal) for all shells (which happens by
default).  The one free norm applies to gas occupying the intersection
between the innermost spherical shell and the 3-dimensional cylinder
corresponding to the innermost annulus.


## Model parameters

For clmass:
* rinner - inner radius of innermost shell (same units as XFLT0001)
* a      - core radius for beta model (same units as XFLT0001)
* beta   - beta model exponent
* switch - 1 enable, 0 disable beta model
* kTa... - shell temperatures - must be linked to the corresponding thermal component, or frozen
* dena... - gravitating matter densities for the shells - unused densities must be frozen, including the density for the outermost shell when the beta model is used

For monomass:
* rinner - inner radius of innermost shell (same units as XFLT0001)
* a      - core radius for beta model (same units as XFLT0001)
* beta   - beta model exponent
* switch - 1 enable, 0 disable beta model
* kTa... - shell temperatures - must be linked to the corresponding thermal component, or frozen
* dela... - gravitating matter density differences for the shells.  In terms of clmass the parameters, dela = dena - denb, delb = denb - denc, etc.  Unused parameters must be frozen, including that for the outermost shell when the beta model is used 

For nfwmass:
* rinner - inner radius of innermost shell (same units as XFLT0001)
* a      - core radius for beta model (same units as XFLT0001)
* beta   - beta model exponent
* switch - 1 enable, 0 disable beta model
* kTa... - shell temperatures - must be linked to the corresponding thermal component, or frozen
* nfwa   - NFW scale length (same units as XFLT0001)
* nfwpot - Normalization for NFW potential


## Units

Lengths can be specified in any unit, but must be consistent.
Internal units depend on the length unit.  If the physical length
corresponding to the unit in XFLT0001 is `u`, then the densities used
by the model are in units of `keV / (G \mu m_H u^2)`, where `keV` is
the energy of 1 keV, `G` is Newton's constant and `mu m_H` is the
mean mass per particle in the gas.

For nfwmass, the normalization constant, nfwpot, is `4 \pi G \rho_0 a^2 \mu m_H` in units of keV.  Here, `a` is the NFW scale length in physical units and `rho_0` is the normalizing density for the NFW potential (mass density is `rho_0 /[r/a (1 + r/a)^2]`).


## Notes

The number of shells available in these models is determined solely by
the number of entries in model.dat.  If you need more shells, simply
add more shell parameters to model.dat (being sure to add them in
pairs for clmass and monomass) and rebuild the models.

Code for the gravitational potentials has been separated from the
remainder of the code to make it relatively easy to add new potential
models.  See the notes with the model code.

## Installing the models

The model code is in model, useful support programs in support, tcl scripts in tcl, and test data and scripts in test. To build these models, go to the model sub-directory and follow the instructions in the file USING. To build the support programs, go to the support sub-directory and type make.
