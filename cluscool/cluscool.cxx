#include <functionMap.h>
#include <xsTypes.h>
#include <FunctionUtility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <iostream>
#include <sstream>


extern "C" void vcph(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString);


extern "C" void cph(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{

  //  XSPEC model subroutine to calculate modified cooling flow spectrum
  //  from sum of spectra.
  //  params array :
  //        0..................peak temperature
  //        1..................abundance
  //        2..................redshift
  //        3..................switch(0=calculate MEKAL model, 
  //                                  1=interpolate MEKAL model
  //                                  2=APEC model)

  //  Norm is mass accretion rate in units of Msun/yr


  RealArray vparams(17);
  vparams[0] = params[0];
  vparams[1] = 1.0;
  for (size_t i=2; i<15; i++) vparams[i] = params[1];
  vparams[15] = params[2];
  vparams[16] = params[3];
  vcph(energyArray, vparams, spectrumNumber, flux, fluxErr, initString);

  return;
}
