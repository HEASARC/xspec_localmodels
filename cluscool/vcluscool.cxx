//=======================This code is for a Cooling+Heating model with constant pressure.======================//

#include <functionMap.h>
#include <xsTypes.h>
#include <FunctionUtility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <iostream>
#include <sstream>

#define keVtoK 1.1604505e7

// function definitions in pfltmp.cxx, pflbol.cxx, and calcMultiTempPlasma.cxx
void pfltmp(int itype, float** tval, int* nmtval, int* status);
void pflbol(int itype, const float* abun, float* bolo, int* status);
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerArray& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);


extern "C" void vcph(const RealArray& energyArray, const RealArray& params,
		     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
		     const string& initString)
{

  //  XSPEC model subroutine to calculate modified cooling flow spectrum
  //  from sum of spectra. Variable abundances.
  //  params array :
  //        0..................peak temperature
  //        1..................He abundance
  //        2..................C     "
  //        3..................N     "
  //        4..................O     "
  //        5..................Ne    "
  //        6..................Na    "
  //        7..................Mg    "
  //        8..................Al    "
  //        9..................Si    "
  //       10..................S     "
  //       11..................Ar    "
  //       12..................Ca    "
  //       13..................Fe    "
  //       14..................Ni    "
  //       15..................redshift
  //       16..................switch(0=calculate MEKAL model, 
  //                                  1=interpolate MEKAL model
  //                                  2=APEC model)

  //  Norm is mass accretion rate in units of Msun/yr


  using namespace XSutility;

  const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
  IntegerArray Zarray(14);
  for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

  const size_t nE = energyArray.size()-1;
  flux.resize(nE);
  fluxErr.resize(0);


  // On initialization find the tabulated temperatures and calculate the
  // bolometric luminosities for each temperature.

  static bool isFirst = true;
  static float *ttab=0; // This memory is owned in ldpfil
  static int nmtval=0;
  static auto_array_ptr<float> apBolo(0);
  int status = 0;   
  int itype = 2;
  if (isFirst) {
    pfltmp(itype, &ttab, &nmtval, &status);
    if (status) return;      
    float *oldF = apBolo.reset(new float[nmtval]);
    delete [] oldF;
  }

  // Calculate the bolometric luminosities
  float *abunParam = new float[14];
  for (size_t i=0; i<14; i++) abunParam[i] = params[i+1];
  pflbol(itype, abunParam, apBolo.get(), &status);
  if (status) return;
  delete [] abunParam;

  // Get the number of temperature steps. 
  // Here is my changes in the nsteps
  Real tlow = 0.01;
  Real thigh = 50.0;
  Real altlow = log(tlow);
  Real althigh = log(thigh);
  float dlogT = 0.1;

  int nsteps = (althigh-altlow)/dlogT;
  

  
  Real slope = 0.0;
  RealArray abun(14);
  for (size_t i=0; i<14; i++) abun[i] = params[i+1];
  Real dens = 1.0;
  Real z = params[15];
  Real tpeak = params[0];
  Real sqrtTpeak5K = pow(keVtoK*tpeak,2.5);
  int switchPar = static_cast<int>(params[16]);

  if ( z <= 0.0 ) {
    FunctionUtility::xsWrite("\n VCLUSCOOL: Require z > 0 for cooling flow models",10);
    return;
  }

  int plasmaType;
  if ( switchPar == 0 || switchPar == 1) {
    plasmaType = switchPar + 3;
  } else if ( switchPar == 2 ) {
    plasmaType = switchPar + 4;
  } else {
    FunctionUtility::xsWrite("\n VCLUSCOOL: Invalid switch parameter value",2);
    FunctionUtility::xsWrite("            Must be 0, 1, or 2",2);
    return;
  }

  // Set up the temperature array


  RealArray tval(nsteps);
  for (size_t i=0; i<nsteps; i++) {
    tval[i] = exp(altlow + i*(althigh-altlow)/(nsteps-1));
  }

  // Get the cosmology parameters

  Real q0 = FunctionUtility::getq0();
  Real h0 = FunctionUtility::getH0();
  Real Lambda0 = FunctionUtility::getlambda0();

  // Fix up the norm. The numerical constant assumes distance linearly depends
  // on redshift with H0=50 hence cosmology factors correct this.

  Numerics::FZSQ fzsq;
  Real norm = 3.16e-15 * (h0/50.)*(h0/50.) / fzsq(z, q0, Lambda0);

  
  
  //  Now calculate the DEMs. The integral is an extended trapezium rule 
  //  over the temperatures. The integration is performed in log T space.

  RealArray dem(nsteps);
  for (size_t i=0; i<nsteps; i++) dem[i] = 0.0;

  // Loop over all the temperatures

  for ( size_t i=0; i<nsteps; i++) {

    Real tkeV = tval[i];

    Real factor = (althigh-altlow)/(nsteps-1);
    if ( i == 0 || i == nsteps-1 ) factor = factor/2.;

    //  Interpolate on the tabulated array of bolometric luminosities to
    //  get the value for this temperature

    Real boloi;
    if ( tkeV <= ttab[0] ) {
      boloi = apBolo.get()[0];
    } else if ( tkeV >= ttab[nmtval-1] ) {
      boloi = apBolo.get()[nmtval-1];
    } else {
      size_t j = 0;
      while ( tkeV > ttab[j+1] ) j++;
      boloi = ( apBolo.get()[j]*(ttab[j+1]-tkeV) + apBolo.get()[j+1]*(tkeV-ttab[j]) ) 
	/ (ttab[j+1]-ttab[j]);
    }

        
    //  calculate the emission-weighting for this temperature. This
    //  includes the division by the bolometric emissivity and an
    //  extra factor of T since we are integrating in log space.

    Real tstep = keVtoK*tval[i];
    dem[i] = pow((tstep/1.e7),slope) * tstep * factor * norm / boloi;

// Here tpeak is T* in Zhoolideh Haghighi et al. 2018.
// Interpolate on the tabulated array of bolometric luminosities to

    Real sqrtTstep5 = pow(tstep,2.5);
    Real tpeakBolo;
    if ( tpeak <= ttab[0] ) {
      tpeakBolo = apBolo.get()[0];
    } else if ( tpeak >= ttab[nmtval-1] ) {
      tpeakBolo = apBolo.get()[nmtval-1];
    } else {
      size_t j = 0;
      while ( tpeak > ttab[j+1] ) j++;
      tpeakBolo = ( apBolo.get()[j]*(ttab[j+1]-tpeak) + apBolo.get()[j+1]*(tpeak-ttab[j]) ) 
	/ (ttab[j+1]-ttab[j]);
    }

    dem[i] *= 1./fabs(1-sqrtTstep5/boloi/(sqrtTpeak5K/tpeakBolo));

  }

  calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z,
		      tval, dem, spectrumNumber, false, 0.0, flux, fluxErr);

   isFirst = false;
}
