#include <iostream>
#include <iomanip>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "param.h"
#include "utilize.h"
#include "dynamics.h"
#include "polar.h"

using namespace std;

double calculateMean(const vector<double>& data) {
	double sum = 0.0;
	for (double num : data) {
		sum += num;
	}
	return sum / data.size();
}

double calculateStandardDeviation(const vector<double>& data) {
	double mean = calculateMean(data);
	double variance = 0.0;

	for (double num : data) {
		variance += pow(num - mean, 2);
	}
	variance /= data.size();
	return sqrt(variance);
}

int main(int argc, char** argv)
{
	string filename(argv[1]);
	cout << "Now is CPU version" << '\n';
	cout << "Running " << filename << '\n';
	param k(filename);
	cout << "The bath type is " << k.bath_type << '\n';
	construct_ADO_set(k);

	int sys_size = k.sys_size;
	int total_size = k.sys_size * k.sys_size;
	gsl_matrix_complex* H = gsl_matrix_complex_alloc(sys_size, sys_size);

	const gsl_rng_type* T;
	gsl_rng* r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, k.seed);

	clock_t start = clock();
	vector<vector<gsl_complex>> p;

	for (int i = 0; i < k.n_sample; i++) {
		cout << "#Now running sample " << i + 1 << '\n';
		cout << "The program is running with importance value cutoff level " << k.cutoff_level << " and " << k.ado.size() << " ADOs" << '\n';
		cout << "The coherent time is " << k.pulses[1].tau0 - k.pulses[0].tau0 << " fs " << "and the population time is ";
		if (k.pulses[1].tau0 >= k.pulses[0].tau0) cout << k.pulses[2].tau0 - k.pulses[1].tau0 << " fs" << '\n';
		else cout << k.pulses[2].tau0 - k.pulses[0].tau0 << " fs" << '\n';

		gsl_matrix_complex_set_all(H, gsl_complex_rect(0.0, 0.0));
		construct_Hal(k, H, r);
		cout << "The disordered Hamiltonain is:" << '\n';
		print_Hal(H);

		polar_mat_set(k);
		polar_mat_ranrot(k, r);
		compute_pulse_interaction(k);
		cout << "The absolute dipole operator X_abs is:" << '\n';
		print_matrix_real(k.X_site);

		vector<gsl_complex> p1;
		vector<gsl_matrix_complex*> op1;
		propagation_Ht(k, H, 1, 1, 1, p1);

		vector<gsl_complex> p2;
		vector<gsl_matrix_complex*> op2;
		propagation_Ht(k, H, 1, 1, 0, p2);

		vector<gsl_complex> p3;
		vector<gsl_matrix_complex*> op3;
		propagation_Ht(k, H, 1, 0, 1, p3);

		vector<gsl_complex> pi;
		for (int j = 0; j < p1.size(); j++) {
			double real = GSL_REAL(p1[j]) - GSL_REAL(p2[j]) - GSL_REAL(p3[j]);
			double imag = GSL_IMAG(p1[j]) - GSL_IMAG(p2[j]) - GSL_IMAG(p3[j]);
			pi.push_back(gsl_complex_rect(real, imag));
		}
		p.push_back(pi);
	}
	clock_t end = clock();

	vector<gsl_complex> P(p[0].size(), gsl_complex_rect(0.0, 0.0));
	for (int i = 0; i < p.size(); i++) {
		for (int j = 0; j < p[i].size(); j++) {
			P[j].dat[0] += GSL_REAL(p[i][j]) / k.n_sample;
			P[j].dat[1] += GSL_IMAG(p[i][j]) / k.n_sample;
		}
	}

	string output = "out";
	string delimiter = "_";
	size_t pos = 0;
	std::string token;

	pos = filename.find(delimiter);
	filename.erase(0, pos + delimiter.length());

	while ((pos = filename.find(delimiter)) != std::string::npos) {
		token = filename.substr(0, pos);
		output += "_" + token;
		filename.erase(0, pos + delimiter.length());
	}

	pos = filename.find(".");
	if (pos != string::npos) {
		token = filename.substr(0, pos);
		output += "_" + token;
	}

	string f1;
	f1 = output + ".txt";
	ofstream file1("../2d output/" + f1);

	for (int i = 0; i < P.size(); i++) file1 << scientific << setprecision(6) << GSL_REAL(P[i]) << " " << GSL_IMAG(P[i]) << '\n';


	file1.close();
	k.param_free();
	return 0;
}


