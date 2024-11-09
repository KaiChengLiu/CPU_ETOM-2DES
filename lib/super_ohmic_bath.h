#ifndef SUPER_OHMIC_BATH_H
#define SUPER_OHMIC_BATH_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include "constant.h"

/* Here we use  E. Mangauda, C. Meiera and M. Desouter-Lecomte parameterization of the super Ohmic spectral density.
   J(w) = eta / 3! * w^3 / wc^2 * exp(-w/wc)
	   = \sum_{k=1..4} p_k*w/{[(w+Omega1_k)^2+Gamma1_k^2]*[(w-Omega1_k)^2+Gamma1_k^2]}
   where

   k    pk/(eta*wc^4)  Gamma_k/wc  Omega_k/wc
   1		10.6		  2.48        2.20
   2		7.00		  4.23        2.25
   3		-0.300		  0.0581      1.32
   4		-6.47		  9.87		  3.49
*/

using namespace std;

void set_super_ohmic_bath(vector<gsl_matrix_complex*>& S, vector<vector<gsl_complex>>& Gamma,
	vector<double>& Lambda, vector<vector<gsl_complex>>& Alpha, vector<vector<gsl_complex>>& Alpha_t,
	double beta, double coupling_str, double wc, int site, int K_m, int sys_size);

#endif
