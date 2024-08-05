#include <iostream>
#include <vector>

#include <random>
#include <numeric>
#include <complex>

template<typename T1, typename T2, typename T3, typename T4, typename T5>
class VEGASResult {
public:
    T1 integral_estimate;
    T2 standard_deviation;
    T3 chi_squared_average;
    T4 adaptive_grid;
    T5 grid_spacing;

    VEGASResult(T1 integral_estimate, T2 standard_deviation, T3 chi_squared_average, T4 adaptive_grid, T5 grid_spacing)
        : integral_estimate(integral_estimate), standard_deviation(standard_deviation), chi_squared_average(chi_squared_average), adaptive_grid(adaptive_grid), grid_spacing(grid_spacing) {}
};

static VEGASResult<std::complex<double>, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>>
Vegas(std::complex<double> fxn(const std::vector<double>&), const std::vector<double> lb, const std::vector<double> ub,
    int maxiter = 10, const int nbins = 1000, int ncalls = 10000,
    double alpha = 1.5, const int noncumulative = 10, const int adaptive = 10) {

    const int ndim = lb.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<std::vector<double>> x(nbins + 1, std::vector<double>(ndim));
    std::vector<std::vector<double>> delx(nbins, std::vector<double>(ndim));

    for (int dim = 0; dim < ndim; ++dim) {
        double step = (ub[dim] - lb[dim]) / nbins;
        for (int i = 0; i <= nbins; ++i) {
            x[i][dim] = lb[dim] + step * i;
            if (i < nbins) {
                delx[i][dim] = step / nbins;
            }
        }
    }

    std::complex<double> Itot = 0.0;
    double suminversigma = 0.0;
    double sd = 0.0;
    std::vector<std::complex<double>> integrals;
    std::vector<double> sigma_squares;
    int iter = 1;
    std::vector<std::vector<double>> ymat(ncalls, std::vector<double>(ndim));
    std::vector<std::vector<double>> xmat(ncalls, std::vector<double>(ndim));
    auto y2i = [&](double y) { return static_cast<int>(std::floor(nbins * y)); };
    auto delta = [&](double y) { return nbins * y - y2i(y); };
    auto y2x = [&](double y, int dim) { return x[y2i(y)][dim] + delx[y2i(y)][dim] * delta(y); };
    auto J = [&](double y, int dim) { return nbins * delx[y2i(y)][dim]; };
    std::vector<std::vector<int>> imat(ncalls, std::vector<int>(ndim));
    std::vector<double> Js(ncalls);
    std::vector<std::complex<double>> Jsf(ncalls);
    while (iter <= maxiter) {
        for (int i = 0; i < ncalls; i++) {
            Js[i] = 1.0;
            for (int dim = 0; dim < ndim; dim++) {
                ymat[i][dim] = dis(gen);
                xmat[i][dim] = y2x(ymat[i][dim], dim);
                Js[i] *= J(ymat[i][dim], dim);
                imat[i][dim] = y2i(ymat[i][dim]);
            }
            std::complex<double>A(Js[i], 0.0);
            std::complex<double>B = fxn(xmat[i]);
            Jsf[i] =  A* B;
        }
        std::complex<double> integral_mc = std::accumulate(Jsf.begin(), Jsf.end(), std::complex<double>(0.0, 0.0)) / std::complex<double> ( ncalls,0);
        double variance_mc = (std::accumulate(Jsf.begin(), Jsf.end(), 0., [](double sum, std::complex<double> x) { return sum + abs(x*x); }) / static_cast<double>(ncalls) - abs(integral_mc*integral_mc)) / (ncalls - 1.) + std::numeric_limits<double>::epsilon();
      
        integrals.push_back(integral_mc);
        sigma_squares.push_back(variance_mc);

        if (iter <= adaptive) {
            std::vector<std::vector<double>> d(nbins, std::vector<double>(ndim, 0.0));
            double ni = static_cast<double>(ncalls) / nbins;
            for (int i = 0; i < ncalls; ++i) {
                for (int dim = 0; dim < ndim; ++dim) {
                    d[imat[i][dim]][dim] += (std::norm(Jsf[i]) / ni);
                }
            }
            std::vector<std::vector<double>> dreg = d;
            for (int dim = 0; dim < ndim; ++dim) {
                double sumd = std::accumulate(dreg.begin(), dreg.end(), 0.0, [dim](double sum, const std::vector<double>& vec) { return sum + vec[dim]; });
                dreg[0][dim] = (7 * d[0][dim] + d[1][dim]) / 8;
                for (int j = 1; j < nbins - 1; j++) {
                    dreg[j][dim] = (d[j - 1][dim] + 6 * d[j][dim] + d[j + 1][dim]) / 8;
                }
                dreg[nbins - 1][dim] = (d[nbins - 2][dim] + 7 * d[nbins - 1][dim]) / 8;
                for (int j = 0; j < nbins; j++) {
                    dreg[j][dim] /= sumd;
                }
            }
            for (int dim = 0; dim < ndim; ++dim) {
                for (int j = 0; j < nbins; ++j) {
                    dreg[j][dim] = std::pow((1 - dreg[j][dim]) / std::log(1 / (dreg[j][dim])), alpha);
                }
            }
            std::vector<std::vector<double>> newx(nbins + 1, std::vector<double>(ndim, 0.0));

            for (int dim = 0; dim < ndim; dim++) {
                double Sd = 0.0;
                double delta_d = std::accumulate(dreg.begin(), dreg.end(), 0.0, [dim](double sum, const std::vector<double>& vec) { return sum + vec[dim]; }) / nbins;
                int i = 0, j = 0;
                newx[0][dim] = x[0][dim];
                newx[nbins][dim] = x[nbins][dim];
                i++;
                while (i < nbins) {
                    while (Sd < delta_d) {
                        Sd += dreg[j][dim];
                        j++;
                    }
                    Sd -= delta_d;
                    newx[i][dim] = x[j][dim] - (Sd / dreg[j - 1][dim]) * (delx[j - 1][dim]);
                    i++;
                }
            }
            x = newx;
            for (int i = 0; i < nbins; ++i) {
                for (int dim = 0; dim < ndim; ++dim) {
                    delx[i][dim] = x[i + 1][dim] - x[i][dim];
                }
            }
        }

        if (iter > noncumulative) {
             suminversigma = std::accumulate(sigma_squares.begin() + noncumulative, sigma_squares.end(), 0.0, [](double Sum, double XX) { return Sum + 1. / XX; })+ std::numeric_limits<double>::epsilon();
            Itot = std::inner_product(std::next(integrals.begin(), noncumulative), integrals.end(), std::next(sigma_squares.begin(), noncumulative), std::complex<double>(0.0, 0.0),
                std::plus<>(),
                [](std::complex<double> a, double b) { return a / b; }) / suminversigma;
            sd = 1. / sqrt(suminversigma);
            double chi_squared = 0.0;
            for (size_t i = noncumulative; i < integrals.size(); ++i) {
                chi_squared += std::norm(integrals[i] - Itot) / sigma_squares[i];
            }

            std::cout << iter << "th iteration" << "\nNon-cumulative estimate = " << integral_mc << std::endl;
            std::cout << "Cumulative estimate = " << Itot << std::endl;
            std::cout << "Sd = " << sd << std::endl;
            std::cout << "Chi/iter-1 = " << chi_squared / (iter - 1 - noncumulative + std::numeric_limits<double>::epsilon()) << std::endl;
        }

        iter++;
    }
    double chi_squared = 0.0;
    for (size_t i = noncumulative; i < integrals.size(); ++i) {
        chi_squared += std::norm(integrals[i] - Itot) / sigma_squares[i];
    }

    return VEGASResult<std::complex<double>, double, double, std::vector<std::vector<double>>, std::vector<std::vector<double>>>(Itot, sd, chi_squared / (iter - 2 - noncumulative), x, delx);
}
