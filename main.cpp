#include <vector>
#include <iostream>
#include <iomanip>

#include "Vegas.cpp"
double torusfunc(std::vector<double> const &x) {
	//return 1 / (pow(x[0] - 0.5, 2) + 1e-12) + 1 / ((pow(x[1] - 0.7, 2) + 1e-12) * (pow(x[2] - 0.1, 2) + 1e-12) * (pow(x[0] - 0.1, 2) + 1e-12));
	
	double S = 0.,P=1.;
	for (int i = 4; i <=19; i++) {
		S = S + x[i];
	}
	for (int i = 20; i <= 29; i++) {
		P= P* x[i];
	}
	return 4 * x[0] * x[2] * x[2] * exp(2 * x[0] * x[2]) * exp(S) * P / ((1 + x[1] + x[3]) * (1 + x[1] + x[3]));
	//return x[1] * x[1] / (pow(x[0] - 0.5, 2) + 1e-12) + x[0] * x[0] / (pow(x[1] - 0.5, 2) + 1e-12);
}


int main() {
	
	std::vector<double> lb(30,0);
	std::vector<double> ub(30,1);
	VEGASResult<double, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>>
		A = Vegas(torusfunc, lb, ub, 20, 1000,10000, 1e-4, 1e-4, 0.5, 10,20);
	std::cout << "-----------------------------------------\n";
	std::cout << "Uoc tinh tich luy cuoi cung\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.integral_estimate << std::endl;

	std::cout << "Sai so\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.standard_deviation << std::endl;

	std::cout << "Chi binh/iter-1\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.chi_squared_average << std::endl;

}
