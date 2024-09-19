
#include <vector>
#include "vegas.h"
double torusfunc(std::vector<double> const x, const double wgt) {
	/*
	double SUM = 0, P = 1;
	for (int i = 4; i <= 19; i++) {
		SUM += x[i];
	}
	for (int i = 20; i <= 29; i++) {
		P*= x[i];
	}
	return 4*x[0]*x[2]*x[2]*exp(2*x[0]*x[2]) * exp(SUM) * P / pow((1 + x[1] + x[3]),2);
	*/
	//return 1 / (pow(x[0] - 0.5, 2) + 1e-12) + 1 / ((pow(x[1] - 0.7, 2) + 1e-12) * (pow(x[2] - 0.1, 2) + 1e-12) * (pow(x[0] - 0.7, 2) + 1e-12));
	return 1 / (pow(x[0] + x[1] - 0.5, 2) + 1e-12);
}
int main() {
	int const dim = 2;
	double tgral, sd, chi2a;
	std::vector<double> regn(dim*2);
	for (int i = 0; i < dim*2; i++) {
		if( i< dim) regn[i] = 0.;else regn[i] = 1.;
	}

	vegas(regn, torusfunc, 0, 1000*1000, 10, 0, tgral, sd, chi2a);
	vegas(regn, torusfunc, 1, 1000*1000, 10, 0, tgral, sd, chi2a);
	vegas(regn, torusfunc, 2, 1000*1000, 5, 0, tgral, sd, chi2a);
	


	return 0;
}
