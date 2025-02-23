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


    /*
    // 4 modes
    a.push_back(gsl_complex_rect(-0.02577 / 2, -0.003757 / 2));
    a.push_back(gsl_complex_rect(-0.02577 / 2, 0.003757 / 2));
    g.push_back(gsl_complex_rect(0.08242264, 0.01481850));
    g.push_back(gsl_complex_rect(0.08242264, -0.01481850));

    a.push_back(gsl_complex_rect(-0.14826700 / 2, -0.02133838 / 2));
    a.push_back(gsl_complex_rect(-0.14826700 / 2, 0.02133838 / 2));
    g.push_back(gsl_complex_rect(0.05997344, 0.00681667));
    g.push_back(gsl_complex_rect(0.05997344, -0.00681667));

    a.push_back(gsl_complex_rect(-0.00031979 / 2, 0.00014004 / 2));
    a.push_back(gsl_complex_rect(-0.00031979 / 2, -0.00014004 / 2));
    g.push_back(gsl_complex_rect(0.02051238, 0.01309151));
    g.push_back(gsl_complex_rect(0.02051238, -0.01309151));

    a.push_back(gsl_complex_rect(0.17457964 / 2, 0.13968056 / 2));
    a.push_back(gsl_complex_rect(0.17457964 / 2, -0.13968056 / 2));
    g.push_back(gsl_complex_rect(0.06313848, 0.00147159));
    g.push_back(gsl_complex_rect(0.06313848, -0.00147159));
    */
    /*
    // 3 modes
    a.push_back(gsl_complex_rect(0.05308145, -0.10274389));
    a.push_back(gsl_complex_rect(0.05308145, 0.10274389));
    g.push_back(gsl_complex_rect(0.01269916, -0.01592092));
    g.push_back(gsl_complex_rect(0.01269916, 0.01592092));

    a.push_back(gsl_complex_rect(0.05767212, -0.11226478));
    a.push_back(gsl_complex_rect(0.05767212, 0.11226478));
    g.push_back(gsl_complex_rect(0.01327437, -0.01579495));
    g.push_back(gsl_complex_rect(0.01327437, 0.01579495));

    a.push_back(gsl_complex_rect(-0.11062709, 0.21471273));
    a.push_back(gsl_complex_rect(-0.11062709, -0.21471273));
    g.push_back(gsl_complex_rect(0.01298316, -0.01585916));
    g.push_back(gsl_complex_rect(0.01298316, 0.01585916));
    */
    /*
    // 2 modes
    a.push_back(gsl_complex_rect(0.00005415, 0.77408592));
    a.push_back(gsl_complex_rect(0.00005415, -0.77408592));
    g.push_back(gsl_complex_rect(0.01509037, 0.00000049));
    g.push_back(gsl_complex_rect(0.01509037, -0.00000049));

    a.push_back(gsl_complex_rect(0.00006917, -0.00011192));
    a.push_back(gsl_complex_rect(0.00006917, 0.00011192));
    g.push_back(gsl_complex_rect(0.02618192, -0.03353527));
    g.push_back(gsl_complex_rect(0.02618192, 0.03353527));
    */
    /*
    // 1 mode
    a.push_back(gsl_complex_rect(0.00012834, 1.58905417));
    a.push_back(gsl_complex_rect(0.00012834, -1.58905417));
    g.push_back(gsl_complex_rect(0.03104773, 0.00000227));
    g.push_back(gsl_complex_rect(0.03104773, -0.00000227));
    */


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
