#include "utilize.h"

void Build_ADO_Ivalue(const param& key, std::string& current, const int cur_L, std::vector<double>& I_val) {
	if (current.length() == key.K * key.K_m) {
		if (cur_L <= key.L) {
			double I = calculate_importance_val(key, current);
			I_val.push_back(I);
		}
		return;
	}

	for (char c = 0; c <= key.L; c++) {
		if (cur_L + c <= key.L) {
			current.push_back(c);
			Build_ADO_Ivalue(key, current, cur_L + c, I_val);
			current.pop_back();
		}
	}
}

void Build_ADO(param& key, std::string& current, const int cur_L, const std::vector<double>& I_val) {
	if (current.length() == key.K * key.K_m) {
		if (cur_L <= key.L) {
			double I = calculate_importance_val(key, current);
			if (I > I_val[key.cutoff_level]) {
				key.ado.push_back(current);
			}
		}
		return;
	}

	for (char c = 0; c <= key.L; c++) {
		if (cur_L + c <= key.L) {
			current.push_back(c);
			Build_ADO(key, current, cur_L + c, I_val);
			current.pop_back();
		}
	}
}

void Build_ADO_map(param& key) {
	for (int i = 0; i < key.ado.size(); i++) {
		key.ado_map[key.ado[i]] = i;
	}
}


double calculate_importance_val(const param& key, std::string arr) {
	double res = 1;
	int arr_level = 0;
	for (char ele : arr) arr_level += (int)ele;
	if (arr_level == 0) return res;

	for (int i = 0; i < arr.length(); i++) {
		double c1 = 0;
		double c2 = 0;
		int idx1 = (int)arr[i] / key.K;
		int idx2 = (int)arr[i] % key.K_m;
		c1 = GSL_REAL(key.alpha[0][idx2]) / gsl_complex_abs(key.gamma[0][idx2]);
		for (int j = 0; j <= idx2; j++) {
			c2 += gsl_complex_abs(key.gamma[0][j]);
		}
		if (c2 != 0) res *= c1 / c2;
		else res *= c1;
	}
	return abs(res);
}

void construct_ADO_set(param& key) {
	key.K_m = key.K_m + key.K_extra;
	string s1;
	string s2;
	std::vector<double> I;
	I.push_back(-1);
	Build_ADO_Ivalue(key, s1, 0, I);
	sort(I.begin(), I.end());
	auto last = unique(I.begin(), I.end());
	I.erase(last, I.end());
	if (key.cutoff_level >= I.size() - 1 || key.cutoff_level < 0) {
		cout << "cutoff level should in the range 0 to " << I.size() - 2 << "\n";
		exit(EXIT_FAILURE);
	}
	Build_ADO(key, s2, 0, I);
	Build_ADO_map(key);
	for (int i = 0; i < key.ado.size(); i++) {
		gsl_matrix_complex* m = gsl_matrix_complex_alloc(key.sys_size, key.sys_size);
		gsl_matrix_complex_set_all(m, gsl_complex_rect(0.0, 0.0));
		key.rho.push_back(m);
	}
	gsl_matrix_complex_set(key.rho[0], 0, 0, gsl_complex_rect(1.0, 0.0));

}

void print_matrix_real(const gsl_matrix_complex* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			gsl_complex ele = gsl_matrix_complex_get(M, i, j);
			printf("%.4f\t", GSL_REAL(ele));
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void print_Hal(const gsl_matrix_complex* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			gsl_complex ele = gsl_matrix_complex_get(M, i, j);
			printf("%.f\t", GSL_REAL(ele) * 5308);
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

void print_matrix(const gsl_matrix* M) {
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			double ele = gsl_matrix_get(M, i, j);
			std::cout << ele << " ";
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

double delta(int a, int b) {
	return a == b ? 1 : 0;
}
