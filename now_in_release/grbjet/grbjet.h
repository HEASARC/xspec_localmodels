/*===========================================================================================*/
/*Function prototypes*/
/*===========================================================================================*/

#define PI 3.14159265358979323846

int luminosity_distance(unsigned dim, const double *x, void* params, unsigned fdim, double *retval);
int surface_flux(unsigned dim, const double *x, void* params, unsigned fdim, double *retval);

double doppler(double gamma, double beta, double cosalpha);
double compute_cosalpha (double dl, double rjet, double thobs, double u, double phi);

double bknpl(double energy, double E0, double p1, double p2, double delta);
double cpl(double energy, double index_pl, double ecut);
double blackbody(double energy, double ktbb);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

double speepest_bknpl (double energy, double gamma, double thobs, double rjetcm, double dlcm, double redshift, double
		model, double E0, double p1, double p2, double delta, double index_pl, double ecut, double ktbb);


double asymp_angterm(double dl, double rjet, double thobs, double costheta, double phi);

void steepest_params(double* array_energy, int nbins, double* params, int npars, double*spectrum);

double steepest_bb (double energy, double gamma, double thobs,
		    double rjetcm, double dlcm, double redshift, double ktbb);

double time_interval(double rjet, double thobs, double thjet, double redshift);

double gamma_vs_theta (double gamma_max, double gamma_min,
		       double p_index, double thjet, double theta);



