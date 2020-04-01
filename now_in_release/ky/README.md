NOTE: These models are now in the Heasoft release, as of Xspec version 12.11.

A set of models from <a
href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2004ApJS..153..205D&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=4230b2429422869">Dovciak
et al. (2004)</a> for accretion disk spectra in the strong gravity
regime. The models are
<p>
<table>
<tr><td width=100>kyrline</td><td width=500>Relativistic line from accretion disc around Kerr black
hole, axisymmetric version.</td></tr>
<tr><td width=100>kyconv</td><td width=500>Axisymmetric convolution model for relativistic
smearing.</td></tr>
</table>

XSPEC12> initpackage ky lmodel.dat /path/to/KY/models

The KYDIR variable should be set so that the models find the FITS files with
the tables. e.g.

XSPEC12> xset KYDIR /path/to/KY/models

The models will also calculate and save the values of the black hole horizon, 
inner edge, and ISCO to KYRH, KYRIN, and KYRMS, respectively. To see these 
values use;
   
XSPEC12> xset


kyrline: black hole accretion disc line emission

Line emission from an accretion disc around a black hole. The broken power-law 
for the radial dependence and the limb darkening/brightening law for the 
emission directionality are used to define the local flux in the spectral line. 
All relativistic effects are taken into account, see 
<a
href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2004ApJS..153..205D&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=4230b2429422869">Dovciak M., Karas V. & Yaqoob T. (2004) ApJS, 153, 205-221.</a>

<table>
<tr><td width=100>par1</td><td width=500>the black hole angular momentum (0 <= a/M <= 1)</td></tr>
<tr><td width=100>par2</td><td width=500>the observer inclination in degrees (0 deg - pole, 90 deg - disc)</td></tr>
<tr><td width=100>par3</td><td width=500>the inner edge of an accretion disc in GM/c^2</td></tr>
<tr><td width=100>par4</td><td width=500>0 - means we always integrate from the disc inner edge, par3; 1 - if the disc inner edge, par3, is below the marginally stable orbit then we integrate emission from above the ISCO only</td></tr>
<tr><td width=100>par5</td><td width=500>the outer edge of an accretion disc in GM/c^2</td></tr>
<tr><td width=100>par6</td><td width=500>the rest energy of the intrinsically narrow spectral line (keV)</td></tr>
<tr><td width=100>par7</td><td width=500>the inner power-law index for the radial dependence of the emissivity that scales as r^(-par7) below the boundary radius, par9</td></tr>
<tr><td width=100>par8</td><td width=500>the outer power-law index for the radial dependence of the emissivity that scales as r^(-par8) above the boundary radius, par9</td></tr>
<tr><td width=100>par9</td><td width=500>the boundary radius (in units of GM/c^2)</td></tr>
<tr><td width=100>par10</td><td width=500>the overall Doppler shift</td></tr>
<tr><td width=100>par11</td><td width=500>defines the emission directionality:
        0 - isotropic emission (local flux ~ 1);
        1 - Laor's limb darkening (local flux ~ 1+2.06*mu_e);
        2 - Haardt's limb brightening (local flux ~ ln[1+1/mu_e])</td></tr>
<tr><td width=100>norm</td><td width=500>photons/cm^2/s in the spectral line</td></tr>
</table>

KYRH (the black hole horizon, r_h), KYRIN (the disc inner edge, r_in) and KYRMS 
(the marginally stable orbit, r_ms, ISCO) are added to the XSPEC internal 
switches. Use xset command to show their current values.


kyconv: black hole accretion disc emission

This convolution model takes an input flux as a definition of the
local flux across the accretion disc around a black hole. The broken power-law 
radial dependence and the limb darkening/brightening law for the emission 
directionality are used to define the local flux. The output
is the total spectrum of an accretion disc. All relativistic effects are taken
into account, see Dovciak M., Karas V. & Yaqoob T. (2004) ApJS, 153, 205-221.

<table>
<tr><td width=100>par1</td><td width=500>the black hole angular momentum (0 <= a/M <= 1)</td></tr>
<tr><td width=100>par2</td><td width=500>the observer inclination in degrees (0 deg - pole, 90 deg - disc)</td></tr>
<tr><td width=100>par3</td><td width=500>the inner edge of an accretion disc in GM/c^2</td></tr>
<tr><td width=100>par4</td><td width=500>0 - means we always integrate from the disc inner edge, par3; 1 - if the disc inner edge, par3, is below the marginally stable orbit then we integrate emission from above the ISCO only</td></tr>
<tr><td width=100>par5</td><td width=500>the outer edge of an accretion disc in GM/c^2</td></tr>
<tr><td width=100>par6</td><td width=500>the inner power-law index for the radial dependence of the emissivity that scales as r^(-par6) below the boundary radius, par8</td></tr>
<tr><td width=100>par7</td><td width=500> the outer power-law index for the radial dependence of the emissivity that scales as r^(-par7) above the boundary radius, par8</td></tr>
<tr><td width=100>par8</td><td width=500>the boundary radius (in units of GM/c^2)</td></tr>
<tr><td width=100>par9</td><td width=500>the overall Doppler shift</td></tr>
<tr><td width=100>par10</td><td width=500>defines the emission directionality: 0 - isotropic emission (local flux ~ 1); 1 - Laor's limb darkening (local flux ~ 1+2.06*mu_e); 2 - Haardt's limb brightening (local flux  ~ ln[1+1/mu_e])</td></tr>
<tr><td width=100>par11</td><td width=500>the number of grid points in the local energy (the energy resolution of the local flux)</td></tr>
<tr><td width=100>par12</td><td width=500>defines how to normalize the spectra (see the norm below)</td></tr>
<tr><td width=100>norm</td><td width=500>if par12 = 0 then norm means photons/cm^2/s in the line; if par12 > 0 then norm means photons/keV/cm^2/s at 1 keV; if par12 < 0 then the spectrum is not renormalized</td></tr>
</table>


KYRH (the black hole horizon, r_h), KYRIN (the disc inner edge, r_in) and KYRMS 
(the marginally stable orbit, r_ms, ISCO) are added to the XSPEC internal 
switches. Use xset command to show their current values.

Note: there are several restrictions that arise from the fact that existing 
XSPEC models are used for definition of the local flux:
- only the energy dependence of the photon flux can be defined by local XSPEC
  models,
- only a certain type of radial dependence of the local photon flux can be
  imposed, a broken power-law radial dependence was chosen,
- there is no intrinsic azimuthal dependence of the local photon flux, the
  only azimuthal dependence comes through limb darkening/brightening law
  (emission angle depends on azimuth),
- the local flux can highly depend on the energy resolution, i.e. on the energy
  binning used, if the energy resolution is not high enough. This is because the
  flux is defined in the centre of each bin. A large number of bins is needed
  for the highly varying local flux.


