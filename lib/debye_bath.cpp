#include "debye_bath.h"

void set_debye_bath(vector<gsl_matrix_complex*>& S, vector<vector<gsl_complex>>& Gamma,
					vector<double>& Lambda, vector<vector<gsl_complex>>& Alpha, vector<vector<gsl_complex>>& Alpha_t,
					double beta, double gamma, double lambda, int site, int K_m, int sys_size) {
	vector<gsl_complex> g;
	g.push_back({ gamma / h_bar, 0.0 });
	for (int j = 1; j < K_m; j++) {
		g.push_back({ 2 * pi * j / h_bar / beta, 0.0 });
	}
	Gamma.push_back(g);
	Lambda.push_back(lambda / h_bar);
	


	vector<gsl_complex> a;
	double z = beta * h_bar * GSL_REAL(Gamma[site][0]) / 2;
	a.push_back(gsl_complex_rect(GSL_REAL(Gamma[site][0]) * Lambda[site] * cos(z) / sin(z), 1 * GSL_REAL(Gamma[site][0]) * Lambda[site]));
	for (int j = 1; j < K_m; j++) {
		a.push_back(gsl_complex_rect(4 * GSL_REAL(Gamma[site][0]) * Lambda[site] / h_bar / beta * GSL_REAL(Gamma[site][j]) / (GSL_REAL(Gamma[site][j]) * GSL_REAL(Gamma[site][j]) - GSL_REAL(Gamma[site][0]) * GSL_REAL(Gamma[site][0])), 0.0));
	}
	Alpha.push_back(a);

	vector<gsl_complex> a_t;
	for (int j = 0; j < K_m; j++) {
		a_t.push_back(gsl_complex_conjugate(a[j]));
	}
	Alpha_t.push_back(a_t);

	gsl_matrix_complex* m = gsl_matrix_complex_alloc(sys_size, sys_size);
	gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
	gsl_matrix_complex_set(m, site + 1, site + 1, gsl_complex_rect(1.0, 0.0));
	S.push_back(m);

}