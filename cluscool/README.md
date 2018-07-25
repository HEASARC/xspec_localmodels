This CpH (Cooling+Heating) model is a modification of mkcflow to
include heating as described in Zhoolideh Haghighi, Afshordi &
Khoshroshahi 2018 (https://arxiv.org/abs/1806.08822). The peak
temperature parameter is the peak of the emission measure distribution
and occurs where the cooling and heating timescales are equal.

The mkcflow model had parameters for low and high temperature but
these are not necessary for the CpH model since the emission measure
does not diverge at low temperatures. The model integrates the
emission measure distribution from 0.01 to 50 keV.


For Cooling+Heating flows model (CpH) the parameters are:

- par1		peakT: peak temperature (keV)
- par2		Abund: abundance relative to Solar
- par3		Redshift
- par4		switch (0 = calculate, 1 = interpolate, 2 = use AtomDB data)
- norm		Mass accretion rate (solar mass/yr)

For the version with variable abundances (vCpH) the parameters are:

- par1		peakT: peak temperature (keV)
- par2		He: He abundance relative to Solar
- par3		C: C abundance relative to Solar
- par4		N: N abundance relative to Solar
- par5		O: O abundance relative to Solar
- par6		Ne: Ne abundance relative to Solar
- par7		Na: Na abundance relative to Solar
- par8		Mg: Mg abundance relative to Solar
- par9		Al: Al abundance relative to Solar
- par10		Si: Si abundance relative to Solar
- par11		S: S abundance relative to Solar
- par12		Ar: Ar abundance relative to Solar
- par13		Ca: Ca abundance relative to Solar
- par14		Fe: Fe abundance relative to Solar
- par15		Ni: Ni abundance relative to Solar
- par16		Redshift
- par17		switch (0 = calculate, 1 = interpolate, 2 = use AtomDB data)
- norm		Mass accretion rate (solar mass/yr)





