#include <vector>
#include <random>
#pragma once

static const int NDMX = 1000, MXDIM = 30, RANSEED = 5330;
static const double ALPH =0.5, TINY = 1.0e-30;
static int i, it, j, k, mds, nd, ndo, ng, npg;
static double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti;
static double tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt;
static std::vector<int> ia(MXDIM), kg(MXDIM);
static std::vector<double> dt(MXDIM), dx(MXDIM), x(MXDIM), xin(NDMX);
static double d[NDMX][MXDIM], di[NDMX][MXDIM], xi[MXDIM][NDMX], r[NDMX];



struct VegasMap
{//Vector chứa thông tin các khoảng trục 
	double  xmap[MXDIM][NDMX];
	//Vector chứa thông tin trọng số trên một trục
	double   wmap[MXDIM][NDMX];
	VegasMap() :
		xmap(), // Khởi tạo x_map với Dim vector con
		wmap() // Khởi tạo del_x với Dim vector con
	{}
};
static VegasMap vegasmap=VegasMap();
//chứa danh sách các thông tin qua mỗi lần lặp
static std::vector<VegasMap> List_Map;
void print_Result(int nprn, int ndim, double tgral, double chi2a, double sd);
void caculator_function(int ndim, std::vector<double> regn, double fxn(std::vector<double> I, const double  A));
void print_InputParameter(const int nprn, int ndim, int itmx, std::vector<double> regn);
void smoothAndFitterAndImproveGrid(int ndim);
//hàm in ra thông tin trục
void print_AxisAndWeight(int dim);
void rebin(const double rc, const int nd, double r[NDMX], std::vector<double> xin, double xi[MXDIM][NDMX], const int j);
void vegas(std::vector<double> regn, double fxn(std::vector<double> I,const double  A), const int init, const int ncall, const int itmx, const int nprn, double & tgral, double & sd, double& chi2a);
// Giari thích tên biến
//  vector regn xác định khối chữ nhật , Dim phần tử  đầu tiên là các tọa độ dưới ,  Dim phần tử  còn lại chứa các tọa độ trên
// itmx là số lần lặp
//ncall là số lượng mẫu tối đa trong một lần đánh giá trên toàn bộ miền
//INDMX so khoang toi da tren moi truc
// MXDIM la so chieu toi da
//ALPHA he so lam nen, ALPHA mac dinh la 1.5 . ALPHA càng lớn thì tốc độ hội tụ càng nhanh nhưng không nhất quán, giảm ALPHA khi thông tin về hàm kém trongg quá trình cải thiện lưới
//TINY để cộng vào kết quả khi xảy ra sai số chặt cụt, ví dụ: số có độ chính xác 4 chữ số sau dấu phẩy, kết quả mong muốn: 0.00001   -> máy tính: 0.0000

//vector di chứa kết quả trọng số các khoảng


