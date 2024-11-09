#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include "param.h"
#include "TD_hamiltonian.h"
#include "utilize.h"
#include "polar.h"

using namespace std;

void total_ADO_dynamics(const param& key, const std::vector<gsl_matrix_complex*>& rho_copy, std::vector<gsl_matrix_complex*>& drho);

void total_ADO_dynamics_Ht(param& key, const std::vector<gsl_matrix_complex*>& rho_copy, std::vector<gsl_matrix_complex*>& drho, const gsl_matrix_complex* Ht);

void HEOM_Solver(param& key, std::vector<int>& sites);

void propagation_Ht(param& key, gsl_matrix_complex* H, const int nv1, const int nv2, const int nv3, vector<gsl_complex>& polarization);

void construct_Hal(param& key, gsl_matrix_complex* Hal, gsl_rng* r);

#endif

