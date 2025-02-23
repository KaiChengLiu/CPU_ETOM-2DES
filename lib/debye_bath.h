#ifndef DEBYE_BATH_H
#define DEBYE_BATH_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "constant.h"

using namespace std;

void set_debye_bath(vector<gsl_matrix_complex*>& S, vector<vector<gsl_complex>>& Gamma,
					vector<double>& Lambda, vector<vector<gsl_complex>>& Alpha, vector<vector<gsl_complex>>& Alpha_t,
					double beta, double gamma, double lambda, int site, int K_m, int sys_size);

#endif
