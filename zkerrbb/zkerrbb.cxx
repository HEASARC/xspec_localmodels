#include <xsTypes.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <memory>

extern "C" void zrunkbb_(float* ear, int& nE, float& eta, float& astar, 
		float& theta, float& mbh, float& mdd, float& dbh, 
		float& fcol, int& rflag, int& lflag, float& zbh, 
                float* photar, float* fluxE, int* nex1, int* nex2);

extern "C" void zkerrbb(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Memory allocation wrapper function for Li-Xin Li's runkbb routine.

   int nE = static_cast<int>(energyArray.size()) - 1;

   // The memory for temporary arrays

   std::unique_ptr<float[]> apFluxE(new float[nE*2]);
   std::unique_ptr<int[]> apNex1(new int[nE*2]);
   std::unique_ptr<int[]> apNex2(new int[nE*2]);

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   std::unique_ptr<float[]> apEar(ear);
   std::unique_ptr<float[]> apPars(pars);
   std::unique_ptr<float[]> apPhotar(photar);
   std::unique_ptr<float[]> apPhoter(photer);

   // Set the physical quantities required by the main routine

   float eta = pars[0];
   float astar = pars[1];
   float theta = pars[2];
   float mbh = pars[3];
   // convert mdd from Msun/yr to 1e18 g/s units required by runkbb
   // (1.99e33/3.1457e7/1e18)
   float mdd = pars[4] * 6.309e7;
   float zbh = pars[5];
   float fcol = pars[6];
   int rflag = (int)round(pars[7]);
   int lflag = (int)round(pars[8]);

   // calculate distance in kpc 
   float q0 = FunctionUtility::getq0();
   float H0 = FunctionUtility::getH0();
   float Lambda0 = FunctionUtility::getlambda0();
   float clight = 2.99792458e5;
   Numerics::FZSQ fzsq;

   float dbh = (clight/H0)*sqrt(fzsq(zbh, q0, Lambda0)) * 1000.0;

   // call main routine

   zrunkbb_(ear, nE, eta, astar, theta, mbh, mdd, dbh, fcol, rflag, lflag,
	    zbh, photar, apFluxE.get(), apNex1.get(), apNex2.get());
   XSFunctions::floatFluxToStl<float>(photar, photer, nE, false, flux, fluxErr);

   // no flux errors associated with this model
   fluxErr = 0.0;
}
