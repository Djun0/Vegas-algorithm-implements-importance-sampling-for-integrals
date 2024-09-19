#include<vector>
#include"vegas.h"
#include <iostream>
#include <iomanip>
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dist(0., 1.);

void rebin(const double rc, const int nd, double r[NDMX], std::vector<double> xin, double xi[MXDIM][NDMX], const int j)
{
	int i, k = 0;
	double dr = 0.0, xn = 0.0, xo = 0.0;
	for (i = 0; i < nd - 1; i++) {
		while (rc > dr)
			dr += r[k++];

		if (k > 1) xo = xi[j][k - 2];
		xn = xi[j][k - 1];
		dr -= rc;
		xin[i] = xn - (xn - xo) * dr / r[k - 1];
	}
	for (i = 0; i < nd - 1; i++) xi[j][i] = xin[i];
	xi[j][nd - 1] = 1.0;
}
void print_Result(int nprn,int ndim,double tgral,double chi2a,double sd) {
	std::cout << std::fixed << std::setprecision(4);
	if (nprn >= 0) {
		std::cout << " iteration no. " << std::setw(3) << (it + 1);
		std::cout << " : integral = " << std::setw(14) << ti;
		std::cout << " +/- " << std::setw(9) << tsi << std::endl;
		std::cout << " all iterations: " << " integral =";
		std::cout << std::setw(14) << tgral << "+-" << std::setw(9) << sd;
		std::cout << " chi**2/IT n =" << std::setw(9) << chi2a << std::endl;

		if (nprn != 0) {
			for (int k = 0; k < ndim; k++) {
				print_AxisAndWeight(k);
			}
		}

	}

}

void print_AxisAndWeight(int dim) {
	
		std::cout << "\n Du lieu truc x " << std::setw(2) << dim << std::endl;
		std::cout << std::setw(29) << "Truc x" << "|Trong so\n";
		std::cout << std::setw(29) << 0 << "|\n";
		for (int n = 0; n < nd; n++) {
			std::cout << std::setw(30) << "|" << di[n][dim] << "\n";
			std::cout << std::setw(29) << xi[dim][n] << "|\n";
		}
	

}

void print_InputParameter(const int nprn,int ndim,int itmx, std::vector<double> regn) {
	if (nprn >= 0) {
		std::cout << " Input parameters for vegas";
		std::cout << " ndim= " << std::setw(4) << ndim;
		std::cout << " ncall= " << std::setw(8) << calls << std::endl;
		std::cout << std::setw(34) << " it=" << std::setw(5) << it;
		std::cout << " itmx=" << std::setw(5) << itmx << std::endl;
		std::cout << std::setw(34) << " nprn=" << std::setw(5) << nprn;
		std::cout << " ALPH=" << std::setw(9) << ALPH << std::endl;
		std::cout << std::setw(34) << " mds=" << std::setw(5) << mds;
		std::cout << " nd=" << std::setw(5) << nd << std::endl;
		for (j = 0; j < ndim; j++) {
			std::cout << std::setw(30) << " x1[" << std::setw(2) << j;
			std::cout << "]= " << std::setw(11) << regn[j] << " xu[";
			std::cout << std::setw(2) << j << "]= ";
			std::cout << std::setw(11) << regn[j + ndim] << std::endl;
		}
	}

}
void caculator_function(int ndim, std::vector<double> regn, double fxn(std::vector<double> I, const double  A))
{
	fb = f2b = 0.0;
	for (k = 0; k < npg; k++) {
		wgt = xjac;
		for (j = 0; j < ndim; j++) {


			xn = (kg[j] - dist(gen)) * dxg + 1.0;

			ia[j] = std::max(std::min(int(xn), NDMX), 1);
			if (ia[j] > 1) {
				xo = xi[j][ia[j] - 1] - xi[j][ia[j] - 2];
				rc = xi[j][ia[j] - 2] + (xn - ia[j]) * xo;
			}
			else {
				xo = xi[j][ia[j] - 1];
				rc = (xn - ia[j]) * xo;
			}
			x[j] = regn[j] + rc * dx[j];
			wgt *= xo * xnd;
		}
		f = wgt * fxn(x, wgt);
		f2 = f * f;
		fb += f;
		f2b += f2;
		for (j = 0; j < ndim; j++) {
			di[ia[j] - 1][j] += f;
			if (mds >= 0) d[ia[j] - 1][j] += f2;
		}
	}
	f2b = sqrt(f2b * npg);
	f2b = (f2b - fb) * (f2b + fb);
	if (f2b <= 0.0) f2b = TINY;
	ti += fb;
	tsi += f2b;
}
void smoothAndFitterAndImproveGrid(int ndim) {
	for (j = 0; j < ndim; j++) {
		//Tinh chỉnh lưới. 
		// Quá trình tinh chỉnh được giảm bớt để tránh những thay đổi nhanh chóng, gây mất ổn định và cũng được nén trong phạm vi bởi số mũ ALPH.
		xo = d[0][j];
		xn = d[1][j];
		d[0][j] = (xo + xn) / 2.0;
		dt[j] = d[0][j];
		for (i = 2; i < nd; i++) {
			rc = xo + xn;
			xo = xn;
			xn = d[i][j];
			d[i - 1][j] = (rc + xn) / 3.0;
			dt[j] += d[i - 1][j];
		}
		d[nd - 1][j] = (xo + xn) / 2.0;
		dt[j] += d[nd - 1][j];
	}
	for (j = 0; j < ndim; j++) {
		rc = 0.0;
		for (i = 0; i < nd; i++) {
			if (d[i][j] < TINY) d[i][j] = TINY;
			r[i] = pow((1.0 - d[i][j] / dt[j]) /
				(log(dt[j]) - log(d[i][j])), ALPH);
			rc += r[i];
		}
		rebin(rc / xnd, nd, r, xin, xi, j);
	}
}
void vegas(std::vector<double> regn, double fxn(std::vector<double> I, const double  A), const int init, const int ncall, const int itmx, const int nprn, double& tgral, double& sd, double& chi2a) {
	clock_t time_req{};
	int ndim = regn.size() / 2;// so chieu
	// Với init n=0, tiến thành khởi tạo lưới và thích ứng không lưu kết quả ước tính tích phân
	if (init <= 0) {
		ndo = 1;
		mds  =0;//mds = 1 hoặc 0; 1 để áp dụng lấy mẫu phân tầng
		
		for (j = 0; j < ndim; j++) xi[j][0] = 1.0;
	}
	if (init <= 1) si = swgt = schi = 0.0;
	if (init <= 2) {
		nd = NDMX;//lấy số khoảng tối đa
		ng = 1;
		if (mds != 0) {//Thiết lập để phân tầng.

			ng = int(pow(ncall / 2.0 + 0.25, 1.0 / ndim));// số khoảng tối ưu
			mds = 1;

			//kiểm tra xem số khoảng tối ưu ngx2 có lơn hơn số khoảng tối đa
			if ((2 * ng - NDMX) >= 0) {
				// tính toán lại ng,nd
				mds = -1;
				npg = ng / NDMX + 1;
				nd = ng / npg;//nd=NDMX+1
				ng = npg * nd;
			}
		}
		
		k = pow(ng, ndim);//for (k = 1, i = 0; i < ndim; i++) k *= ng;

		npg = std::max((int)ncall / k, 2);//M=2*N^n  -> M/N^n ==2

		//tính toán lại số phép tính
		calls = (double)npg * k;
		//gia số cho mỗi khoảng lưới ban đầu
		dxg = 1.0 / ng;
		//moi chinh lai cho chi2
		//for (dv2g = 1., i = 0; i < ndim; i++) dv2g *= dxg;
		//dv2g = sqrt(calls * dv2g) / npg / npg / (npg - 1.0);
		dv2g = (double)1 / (npg - 1.);//edit them vao, an 2 line code tren
		xnd = nd;
		dxg *= xnd;

		// tính toán Jaco
		xjac = 1.0 / calls;
		for (j = 0; j < ndim; j++) {
			dx[j] = regn[j + ndim] - regn[j];// cận trên - cận dưới
			xjac *= dx[j];
		}
		if (nd != ndo) {
			for (i = 0; i < std::max(nd, ndo); i++) r[i] = 1.0;
			for (j = 0; j < ndim; j++)
				rebin(ndo / xnd, nd, r, xin, xi, j);
			ndo = nd;
		}
		// in ra các tham số đầu vào
		print_InputParameter(nprn,ndim,itmx,regn);
		
	}
	
	for (it = 0; it < itmx; it++) {
		
		
		ti = tsi = 0.0;
		//khởi tạo trọng số =0
		for (j = 0; j < ndim; j++) {
			kg[j] = 1;
			for (i = 0; i < nd; i++) d[i][j] = di[i][j] = 0.0;
		}

		for (;;) {
			//tính ti,tsi
			caculator_function(ndim,regn,fxn);
		
			if (mds < 0) {
				//Sử dụng lấy mẫu phân tầng.
					for (j = 0; j < ndim; j++) d[ia[j] - 1][j] += f2b;
			}
			for (k = ndim - 1; k >= 0; k--) {
				kg[k] %= ng;
				if (++kg[k] != 1) break;
			}
			if (k < 0) break;
		}

		tsi *= dv2g; 
			wgt = 1. / tsi;
		si += wgt * ti;
		schi += wgt * ti * ti;
		swgt += wgt;
		tgral = (double)si / swgt;
		chi2a = (double)(schi - si * tgral) / (it - 0.9999);
		if (chi2a < 0.0) chi2a = 0.0;
		sd = sqrt(1.0 / swgt);
		tsi = sqrt(tsi);
		//in ra kết quả mỗi lần lặp và tổng các lần lặp
		print_Result(nprn, ndim, tgral, chi2a, sd);

		//nén và làm mịn trọng số sau đó rebin
		smoothAndFitterAndImproveGrid(ndim);

		
	}
	time_req = clock() - time_req;
	printf("Processor time taken for multiplication: %f "
		"seconds\n",
		(float)time_req / CLOCKS_PER_SEC);
}





