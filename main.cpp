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
	return (4 * x[0] * x[2] * x[2] * exp(2 * x[0] * x[2]) * exp(S) * P / ((1 + x[1] + x[3]) * (1 + x[1] + x[3])))/((pow(x[0]+0.5, 2)+1e-12)* (pow(x[1] - 0.5, 2) + 1e-12)* (pow(x[29] - 0.5, 2) + 1e-12)* (pow(x[15] - 0.5, 2) + 1e-12));
	//return x[1] * x[1] / (pow(x[0] - 0.5, 2) + 1e-12) + x[0] * x[0] / (pow(x[1] - 0.5, 2) + 1e-12);
}
const int DIM = 6; 

const double NORM_AC = 1.0 / std::pow(0.17720931990702889842, DIM);
const double NORM_B = 1.0 / std::pow(0.17724538509027909508, DIM);
double f(const std::vector<double>& x) {
	double dx2a = 0;
	for (int d = 0; d < DIM; ++d) {
		dx2a += std::pow(x[d] - 0.25, 2);
	}
	double dx2b = 0;
	for (int d = 0; d < DIM; ++d) {
		dx2b += std::pow(x[d] - 0.5, 2);
	}
	double dx2c = 0;
	for (int d = 0; d < DIM; ++d) {
		dx2c += std::pow(x[d] - 0.75, 2);
	}

	return (
		std::exp(-100.0 * dx2a) * NORM_AC
		+ std::exp(-100.0 * dx2b) * NORM_B
		+ std::exp(-100.0 * dx2c) * NORM_AC
		) / 3.0;
}
std::complex<double> ff2(const std::vector<double>& x) {
	//return 1 / (pow((sqrt(x[0] * x[0] + x[1] * x[1] ) - 0.5), 2)+1e-5);
	return 1 / (pow(x[0] + x[1] - 0.5, 2) + 1e-12);
}
const double M_PI = 3.141592653589793;
const double r_s = 2.0;
const double k_F = std::pow(9 * M_PI / 4, 1.0 / 3) / r_s;
const double T = 0.02 * k_F * k_F;
const std::complex<double> delta (0.0, 0.002 * k_F * k_F);
const std::complex<double> Omega (0.0,0.0);
const double q = 0.1 * k_F;

std::complex<double> fermi_dirac(std::complex<double> x) {
	if (real(x) > 700.0)
		return std::complex<double>(0, 0);
	else
		return  std::complex<double>(1, 0) / (std::exp(x / T) + std::complex<double>(1,0));
}

double epsilon_k(double kx, double ky, double kz) {
	double k_squared = kx * kx + ky * ky + kz * kz;
	return k_squared - k_F * k_F;
}

double epsilon_k_q(double kx, double ky, double kz) {
	double k_squared = (kx - q) * (kx - q) + (ky-q) * (ky-q) + (kz-q) * (kz-q);
	return k_squared - k_F * k_F;
}

std::complex<double> integrand(const std::vector<double>& x) {

	double kx = x[0];
	double ky = x[1];
	double kz = x[2];
	std::complex<double> eps_k(epsilon_k(kx, ky, kz), 0.);
	std::complex<double> eps_k_q(epsilon_k_q(kx, ky, kz), 0.);

	std::complex<double> num = fermi_dirac(eps_k_q) - fermi_dirac(eps_k);

	std::complex<double> denom = Omega - eps_k_q + eps_k + delta;
	std::complex<double> constant(-2.0 / std::pow(2 * M_PI, 3.0), 0.0);
	return constant * num / denom;
}
std::complex<double> func(const std::vector<double>& x) {

	std::complex<double> A(x[0], 0);
	
	return std::complex<double>(1.,1.) / (pow(A- std::complex<double>(0.5, 0.),2)+ std::complex < double>(1e-12, 0))
		+ std::complex<double>(0., 1.) / (pow(A - std::complex<double>(0.2, 0.), 2) + std::complex < double>(1e-12, 0));
	//return std::cos(std::complex<double>(1., 1.) + A * A) / sin(std::complex<double>(x[0], 2.));
	//return std::complex<double>(1., 0.) /sqrt(A - std::complex<double>(0.5, 0));
}
int main() {
	
	std::vector<double> lb(3,-3*k_F);
	std::vector<double> ub(3, 3 * k_F);
	VEGASResult<std::complex<double>, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>>
	A = Vegas(integrand, lb, ub, 30, 1000,1000000,1.0, 10,10);
	std::cout << "-----------------------------------------\n";
	std::cout << "Uoc tinh tich luy cuoi cung\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.integral_estimate << std::endl;

	std::cout << "Sai so\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.standard_deviation << std::endl;

	std::cout << "Chi binh/iter-1\n";
	std::cout << std::setw(35) << std::fixed << std::setprecision(5) << std::right << A.chi_squared_average << std::endl;

}
