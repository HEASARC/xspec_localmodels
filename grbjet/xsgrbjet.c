#include "cubature.h"
#include "grbjet.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define VERBOSE 0
#define NR_END 1
#define FREE_ARG char*

typedef double Real;

#if defined(PCUBATURE)
#define cubature pcubature
#else
#define cubature hcubature
#endif

/*===========================================================================================*/
double* xsgrbjet(Real* energy, int Nflux, Real* local_parameter, int spectrum, Real* flux, Real* fluxError, char* init)
{
    int ie, ii;
    double gamma, thobs, thjet, redshift, dl, ktbb, r12;
    double model, delta_ene, center_ene, square_law, spec_norm, k_norm;

    double val[2];
    double err[2];
    static double xmin[2];
    static double xmax[2];

    double function_params[15];

    double val_dl;
    double err_dl;

    double bound_left;
    double bound_right;

    double clight = 3 * pow(10, 10);
    double H0 = 70 * pow(10, 5);
    double erg2kev = 6.2415 * pow(10, 8);
    double t_spa;

    uint16_t maxEval = 100;

    /*=====================================================================*/
    /*XSPEC model parameters*/
    /*=====================================================================*/

    thobs = ((double*)local_parameter)[0];
    thjet = ((double*)local_parameter)[1];
    gamma = ((double*)local_parameter)[2];
    r12 = ((double*)local_parameter)[3];

    ktbb = ((double*)local_parameter)[10];

    model = ((double*)local_parameter)[11];

    redshift = ((double*)local_parameter)[12];

    /*=====================================================================*/
    /*Compute the luminosity distance*/
    /*=====================================================================*/

    bound_left = 0;
    bound_right = redshift;

    hcubature(1, luminosity_distance, &redshift, 1, &bound_left, &bound_right,
              maxEval, 0, 1e-10, ERROR_INDIVIDUAL, &val_dl, &err_dl);
    dl = clight / H0 * (1 + redshift) * val_dl;

    function_params[13] = dl;

    /*=====================================================================*/
    /*Parameters to be passed to the integrating function*/
    /*=====================================================================*/

    function_params[0] = thobs / 180 * PI;
    function_params[1] = thjet / 180 * PI;

    for (ii = 2; ii <= 12; ii++) {
        function_params[ii] = local_parameter[ii];
    }

    /*=====================================================================*/
    /*Define the boundaries of integration*/
    /*=====================================================================*/

    xmin[0] = cos(thjet / 180 * PI);
    xmax[0] = 1;

    xmin[1] = 0;
    xmax[1] = 2 * PI;

    /*=====================================================================*/
    /*Compute the duration of the single pulse*/
    /*=====================================================================*/

    t_spa = time_interval(r12, thobs, thjet, redshift);

    /*=====================================================================*/
    /*Jet radius is in units of 10^(12) cm, while 1 Mpc=3.08*10^(24 cm)*/
    /*Ratio (R_0/Dcm)^2= 1.05*10^(-25) (r12/DL_Mpc)^2 */
    /*=====================================================================*/

    square_law = 1.05 * pow(10, -25);

    /*=====================================================================*/
    /*Normalization constant for the comoving frame emissivity*/
    /*=====================================================================*/

    if (model == 1 || model == 2) {
        spec_norm = pow(10, 20);
    }
    else {
        spec_norm = 5.05 * pow(10, 25);
    }

    k_norm = erg2kev * spec_norm * (1 + redshift) * square_law / t_spa;

    /*=====================================================================*/

    for (ie = 0; ie <= Nflux; ie++) {
        if (ie > 0) {
            ii = ie - 1;

            delta_ene = (energy[ie] - energy[ie - 1]);

            center_ene = 0.5 * (energy[ie] + energy[ie - 1]);

            function_params[14] = center_ene;

            /*=====================================================================
              To avoid underflow errors and issues in the numerical integration
              use the steepest descent method for the BB spectrum and higly
              relativistic flow and thobs < thjet
              The factor (rjetcm/dlcm)^2 with physical units is already
              considered
              in the routine steepest_bb this is why there is 1/square_law
              for the BB flux computation
              =====================================================================*/

            if (model == 1 || model == 2 || (model == 3 && gamma < 50) ||
                (model == 3 && thobs > thjet)) {
                hcubature(1, surface_flux, function_params, 2, xmin, xmax, 2500,
                          0, 1e-40, ERROR_INDIVIDUAL, val, err);

                flux[ii] = k_norm * val[0] / center_ene * delta_ene;
            }

            else {
                flux[ii] = k_norm * steepest_bb(center_ene, gamma, thobs, r12,
                                                dl, redshift, ktbb) /
                           center_ene * delta_ene;
            }

            /*=====================================================================*/
        }
    }

    return flux;
}

/*=========================================================================================*/

int surface_flux(unsigned dim, const double* x, void* params, unsigned fdim, double* retval)
{
    double beta, cosalpha, doppler_factor, angular_term, thobs, redshift;
    double gamma_max, gamma, r12, dl, thjet;
    double gmin, p_index, theta;
    double energy;
    double E0, p1, p2, delta;
    double index_pl, ecut;
    double ktbb;
    double spec;
    double value;
    /*    double argument; */
    int model;

    double costheta = x[0];
    double phi = x[1];

    thobs = ((double*)params)[0];
    thjet = ((double*)params)[1];

    gamma_max = ((double*)params)[2];
    r12 = ((double*)params)[3];
    p1 = ((double*)params)[4];
    p2 = ((double*)params)[5];
    E0 = ((double*)params)[6];
    delta = ((double*)params)[7];

    index_pl = ((double*)params)[8];
    ecut = ((double*)params)[9];

    ktbb = ((double*)params)[10];

    model = ((double*)params)[11];

    redshift = ((double*)params)[12];

    dl = ((double*)params)[13];

    energy = ((double*)params)[14];

    /*=====================================================================*/

    cosalpha = compute_cosalpha(dl, r12, thobs, costheta, phi);

    if (cosalpha <= 0) {
        return 0;
    }

    /*=====================================================================*/

    gmin = ((double*)params)[15];
    p_index = ((double*)params)[16];

    if (gmin > 1) {
        theta = acos(costheta);

        gamma = gamma_vs_theta(gamma_max, gmin, p_index, thjet, theta);
    }
    else {
        gamma = gamma_max;
    }

    /*=====================================================================*/
    /*Compute the Doppler factor*/
    /*=====================================================================*/

    beta = sqrt(gamma * gamma - 1) / gamma;
    doppler_factor = doppler(gamma, beta, cosalpha);

    /*=====================================================================*/
    angular_term = asymp_angterm(dl, r12, thobs, costheta, phi);
    /*=====================================================================*/
    if (model == 1) {
        spec = bknpl(energy / doppler_factor * (1 + redshift), E0, p1, p2, delta);
    }
    else if (model == 2) {
        spec = cpl(energy / doppler_factor * (1 + redshift), index_pl, ecut);
    }
    else {
        spec = blackbody(energy / doppler_factor * (1 + redshift), ktbb);

	/*        argument = 3 * log(energy / doppler_factor * (1 + redshift)) -
		  energy / (doppler_factor * ktbb) * (1 + redshift); */
    }

    value = spec * angular_term * pow(doppler_factor, 3);

    *retval = value;

    return EXIT_SUCCESS;
}
