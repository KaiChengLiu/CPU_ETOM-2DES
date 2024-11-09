#ifndef OHMIC_BATH_H
#define OHMIC_BATH_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "constant.h"

/* Here we use  E. Mangauda, C. Meiera and M. Desouter-Lecomte parameterization of the super Ohmic spectral density.
   J(w) = eta  * w * exp(-w/wc)
	   = \sum_{k=1..3} p_k*w/{[(w+Omega1_k)^2+Gamma1_k^2]*[(w-Omega1_k)^2+Gamma1_k^2]}
   where

   k    pk/(gamma*wc^4)  Gamma_k/wc  Omega_k/wc
   1    12.0677        2.2593      0.2378
   2   -19.9762        5.4377      0.0888
   3     0.1834        0.8099      0.0482
*/

using namespace std;

void set_ohmic_bath(vector<gsl_matrix_complex*>& S, vector<vector<gsl_complex>>& Gamma,
	vector<double>& Lambda, vector<vector<gsl_complex>>& Alpha, vector<vector<gsl_complex>>& Alpha_t,
	double beta, double coupling_str, double wc, int site, int K_m, int sys_size);

#endif
