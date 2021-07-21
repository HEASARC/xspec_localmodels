#include "cubature.h"
#include "grbjet.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================================*/
int luminosity_distance(unsigned dim, const double* x, void* params, unsigned fdim, double* retval)
{
    double omega = 0.29;
    double result;

    result = 1 / sqrt(1 - omega + pow(1 + x[0], 3) * omega);

    *retval = result;

    return 0;
}

/*=======================================================================================*/
double gamma_vs_theta(double gamma_max, double gamma_min, double p_index, double thjet, double theta)
{
    double gamma_new;

    gamma_new = gamma_min +
                (gamma_max - gamma_min) / sqrt(1 + pow(theta / thjet, 2 * p_index));

    return gamma_new;
}

/*=======================================================================================*/
double time_interval(double r12, double thobs, double thjet, double redshift)
{
    double thobs_rad, thjet_rad, rjetcm, clight, res;

    clight = 3 * pow(10, 10);

    thobs_rad = thobs / 180 * PI;
    thjet_rad = thjet / 180 * PI;

    rjetcm = r12 * pow(10, 12);

    if (thobs <= thjet) {
        res = rjetcm * (1 - cos(thobs_rad + thjet_rad)) / clight * (1 + redshift);
    }
    else {
        res = 2 * rjetcm * sin(thobs_rad) * sin(thjet_rad) / clight * (1 + redshift);
    }

    return res;
}


/*====================================================================================*/
double steepest_bb(double energy, double gamma, double thobs, double r12,
		double dl, double redshift, double ktbb)
{
    double a, beta, result, dlcm, rjetcm, E;

    dlcm = dl * 3.086 * pow(10, 24);
    rjetcm = r12 * pow(10, 12);

    a = rjetcm / dlcm;

    beta = sqrt(gamma * gamma - 1) / gamma;
    E = 2.73;

    thobs = thobs / 180 * PI;

    result =
        (2 * pow(a, 2) * pow(energy, 3) * PI * pow(1 + redshift, 3) * sin(thobs)) /
        ((-1 + pow(E, (((-1 + a) * beta * dlcm + sqrt(pow(-1 + a, 2) * pow(dlcm, 2))) *
                       energy * gamma * (1 + redshift)) /
                          (sqrt(pow(-1 + a, 2) * pow(dlcm, 2)) * ktbb))) *
         sqrt((2 * a * beta * energy * gamma * (1 + redshift) *
                   (-2 * (-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb)) * ktbb +
                    beta * energy * gamma * (1 + redshift) +
                    ((-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb)) * ktbb -
                     beta * energy * gamma * (1 + redshift)) *
                        cos(2 * thobs)) +
               ((-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb)) * ktbb -
                beta * energy * gamma * (1 + redshift)) *
                   ((-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb)) *
                        ktbb * pow(cos(thobs), 2) +
                    (2 * (-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb)) * ktbb -
                     beta * energy * gamma * (1 + redshift)) *
                        pow(sin(thobs), 2))) /
              (pow(-1 + pow(E, ((-1 + beta) * energy * gamma * (1 + redshift)) / ktbb), 2) *
               pow(ktbb, 2))));

    return result;
}


/*======================================================================== */
double doppler(double gamma, double beta, double cosalpha)
{
    double value = 1 / (gamma * (1 - beta * cosalpha));

    return value;
}


/* ==========================================================================================*/
double compute_cosalpha(double dl, double rjet, double thobs, double u, double phi)
{
    double value;

    rjet = rjet * pow(10, 12);
    dl = dl * 3.086 * pow(10, 24);

    value = -((rjet - dl * u * cos(thobs) + dl * sqrt(1 - pow(u, 2)) * sin(phi) * sin(thobs)) /
              sqrt(pow(dl, 2) + pow(rjet, 2) - 2 * dl * u * cos(thobs) +
                   2 * dl * rjet * sqrt(1 - pow(u, 2)) * sin(phi) * sin(thobs)));

    return value;
}

/* ==========================================================================================
 */
double asymp_angterm(double dl, double rjet, double thobs, double costheta, double phi)
{
    /* The luminosity distance dl is in Mpc, the jet radius in units of 10^(12) cm*/

    double value;

    value = pow(rjet / dl, 2) * (costheta * cos(thobs) -
                                 sqrt(1 - costheta * costheta) * sin(thobs) * sin(phi));

    return value;
}

/* ==========================================================================================*/
double blackbody(double energy, double ktbb)
{
    double value;
    double argument;

    if (energy / ktbb <= 10) {
        value = pow(energy, 3) / (exp(energy / ktbb) - 1);
    }
    else {
        argument = 3 * log(energy) - energy / ktbb;
        value = exp(argument);
    }

    return value;
}

/* ==================================================================== */
double cpl(double energy, double index_pl, double ecut)
{
    double value;

    value = pow(energy, -index_pl) * exp(-energy / ecut);

    return value;
}

/* ==================================================================== */
double bknpl(double energy, double E0, double p1, double p2, double delta)
{
    double value;

    value = pow(2, -delta * (p1 - p2)) *
            pow(1 + pow(energy / E0, 1 / delta), delta * (p1 - p2)) *
            pow(energy / E0, -p1);

    return value;
}

/* ==================================================================== */
