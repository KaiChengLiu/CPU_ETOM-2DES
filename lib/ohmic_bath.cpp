#include "ohmic_bath.h"

double k = h_bar * h_bar * h_bar * h_bar;

const double ohmic_params[3][3] =
/*p_k			Gamma_k		Omega_k*/
{
{12.0677,     2.2593,    0.2378},
{-19.976,     5.4377,    0.0888},
{0.1834,   0.8099,    0.0482} };

double ohmic_bath_J(double omega, double coupling_str, double wc) {
    double pk, Gk, Ok;
    int i;
    double sum;

    sum = 0.0;
    for (i = 0; i < 3; i++) {
        pk = ohmic_params[i][0] * coupling_str * wc * wc * wc * wc;
        Gk = ohmic_params[i][1] * wc;
        Ok = ohmic_params[i][2] * wc;
        sum = sum + pk * omega / (((omega + Ok) * (omega + Ok) + Gk * Gk) * ((omega - Ok) * (omega - Ok) + Gk * Gk));
    }

    return sum;
}

/* J(iw) = i*Ji(w) */
double ohmic_bath_Ji(double omega, double coupling_str, double wc) {
    double pk, Gk, Ok, dtmp;
    int i;
    double sum;

    sum = 0.0;
    for (i = 0; i < 3; i++) {
        pk = ohmic_params[i][0] * coupling_str * wc * wc * wc * wc;
        Gk = ohmic_params[i][1] * wc;
        Ok = ohmic_params[i][2] * wc;
        dtmp = Gk * Gk + Ok * Ok - omega * omega;
        sum += pk * omega / (dtmp * dtmp + 4 * omega * omega * Ok * Ok);
    }

    return sum;
}

void set_ohmic_bath(vector<gsl_matrix_complex*>& S, vector<vector<gsl_complex>>& Gamma,
    vector<double>& Lambda, vector<vector<gsl_complex>>& Alpha, vector<vector<gsl_complex>>& Alpha_t,
    double beta, double coupling_str, double wc, int site, int K_m, int sys_size) {
    double pk, Gk, Ok;
    gsl_complex ztmp;
    double max;

    double l = coupling_str * wc / pi / h_bar;
    Lambda.push_back(l);
    vector<gsl_complex> a;
    vector<gsl_complex> g;

    for (int i = 0; i < 3; i++) {
        pk = ohmic_params[i][0] * coupling_str * wc * wc * wc * wc;
        Gk = ohmic_params[i][1] * wc;
        Ok = ohmic_params[i][2] * wc;
        ztmp = gsl_complex_coth(gsl_complex_rect(beta * Ok * h_bar / 2.0, -beta * Gk * h_bar / 2.0));
        ztmp = gsl_complex_add(ztmp, gsl_complex_rect(1.0, 0.0));
        ztmp = gsl_complex_mul_real(ztmp, pk / Ok / Gk / 8.0);
        a.push_back(ztmp);
        g.push_back(gsl_complex_rect(Gk, Ok));

        ztmp = gsl_complex_coth(gsl_complex_rect(beta * Ok * h_bar / 2.0, beta * Gk * h_bar / 2.0));
        ztmp = gsl_complex_add(ztmp, gsl_complex_rect(-1.0, 0.0));
        ztmp = gsl_complex_mul_real(ztmp, pk / Ok / Gk / 8.0);
        a.push_back(ztmp);
        g.push_back(gsl_complex_rect(Gk, -Ok));
    }

    for (int i = 0; i < K_m; i++) {
        double mu = 2 * pi * (i + 1) / h_bar / beta;
        double j = ohmic_bath_Ji(mu, coupling_str, wc);
        a.push_back(gsl_complex_rect(-2 * j / beta / h_bar, 0.0));
        g.push_back(gsl_complex_rect(mu, 0.0));
    }

    Alpha.push_back(a);
    Gamma.push_back(g);

    vector<gsl_complex> a_t;
    for (int i = 0; i < a.size(); i++) {
        a_t.push_back(gsl_complex_conjugate(a[i]));
    }
    Alpha_t.push_back(a_t);

    gsl_matrix_complex* m = gsl_matrix_complex_alloc(sys_size, sys_size);
    gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
    gsl_matrix_complex_set(m, site + 1, site + 1, gsl_complex_rect(1.0, 0.0));
    S.push_back(m);

}
