//////////////////////////////////////////////////////////////////
// Functions for calculating the Green's function in planar DLW
// Using the transverse operator method (Baturin 2013 PRSTAB)
// Conversion to cylindrical coordinates for use with DiWaCAT fields
// T Overton 2019
//////////////////////////////////////////////////////////////////


#ifndef Planar_GF
#define _USE_MATH_DEFINES
#include <ctime>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>
#include <fstream>
#include <numeric>
#include <omp.h>


using std::vector;

// ---------------------------------------------------- //
// -----  Functions for Eigenmode Calculations  ------- //
// ---------------------------------------------------- //

double f1(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)  //Rescaled LSE Symmetric Function
{//Rescaled LSE Symmetric Dispersion Equation
	double kx = M_PI*n / w;
	double kch = sqrt((1.0 - B2)*(kmd / (c - b)*kmd / (c - b) + kx*kx) / (Ep*Mu*B2 - 1.0) + kx*kx);
	double kchdkmd = sqrt((1.0 - B2)*(1 / (c - b) * 1 / (c - b) + kx / kmd*kx / kmd) / (Ep*Mu*B2 - 1.0) + kx / kmd*kx / kmd);
	return exp(-b*kch)*(Ep*kchdkmd*sinh(kch*b)*cos(kmd) - 1 / (c - b)*sin(kmd)*cosh(kch*b));
}
double f2(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSE Asymmetric Dispersion Equation
	double kx = M_PI*n / w;
	double kch = sqrt((1 - B2)*(kmd / (c - b)*kmd / (c - b) + kx*kx) / (Ep*Mu*B2 - 1) + kx*kx);
	double kchdkmd = sqrt((1.0 - B2)*(1 / (c - b) * 1 / (c - b) + kx / kmd*kx / kmd) / (Ep*Mu*B2 - 1.0) + kx / kmd*kx / kmd);
	return exp(-b*kch)*(Ep*kchdkmd*cosh(kch*b)*cos(kmd) - 1 / (c - b)*sin(kmd)*sinh(kch*b));
}

double f3(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSM Symmetric Dispersion Equation
	double kx = M_PI*n / w;
	double kch = sqrt((1 - B2)*(kmd / (c - b)*kmd / (c - b) + kx*kx) / (Ep*Mu*B2 - 1) + kx*kx);
	double kchdkmd = sqrt((1.0 - B2)*(1 / (c - b) * 1 / (c - b) + kx / kmd*kx / kmd) / (Ep*Mu*B2 - 1.0) + kx / kmd*kx / kmd);

	return exp(-b*kch)*(Mu*kchdkmd*sinh(kch*b)*sin(kmd) + 1 / (c - b)*cos(kmd)*cosh(kch*b));
}
double f4(double kmd, double n, double Ep, double Mu, double B2, double b, double c, double w)
{//Rescaled LSM Asymmetric Dispersion Equation
	double kx = M_PI*n / w;
	double kch = sqrt((1 - B2)*(kmd / (c - b)*kmd / (c - b) + kx*kx) / (Ep*Mu*B2 - 1) + kx*kx);
	double kchdkmd = sqrt((1.0 - B2)*(1 / (c - b) * 1 / (c - b) + kx / kmd*kx / kmd) / (Ep*Mu*B2 - 1.0) + kx / kmd*kx / kmd);
	return exp(-b*kch)*(1 / Mu / kchdkmd / (c - b)*sinh(kch*b)*cos(kmd) + sin(kmd)*cosh(kch*b));
}

double YE_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//LSE Symmetric Eigenfunction (y-component)
	double YEs = 0.0;
	if (y >= -c && y < -b)
	{
		YEs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YEs;
}
double YE_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//LSE Asymmetric Eigenfunction (y-component)
	double YEa = 0.0;
	if (y >= -c && y < -b)
	{
		YEa = -1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return  YEa;
}

double dyYEc_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{  //Derivative of LSE Symmetric Eigenfunction (y-component), placeholder for complex conjugate
	double YEs = 0.0;
	if (y >= -c && y < -b)
	{
		YEs = -kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEs = kch*1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEs = kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YEs;
}
double dyYEc_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Derivative of LSE Asymmetric Eigenfunction (y-component), placeholder for complex conjugate
	double YEa = 0.0;
	if (y >= -c && y < -b)
	{
		YEa = kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEa = kch*1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEa = kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YEa;
}

double IYE_S(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{ //Antiderivative of LSE Symmetric Eigenfunction (y-component)
	double YEs = 0.0;
	if (y >= -c && y < -b)
	{
		YEs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEs = 1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEs = -1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c - y));

	}
	return YEs;
}
double IYE_A(double y, double lambE, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{ //Antiderivative of LSE Asymmetric Eigenfunction (y-component)
	double YEa = 0.0;
	if (y >= -c && y < -b)
	{
		YEa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YEa = 1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YEa = -1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Ep*cos(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YEa;
}

double XE(double x, int n, double w)
{//x-dependence of LSE eigenfunctions
	return sin((M_PI*n / w)*x);
}
double dxXE(double x, int n, double w)
{//derivative of x-dependence of LSE eigenfunctions
	return M_PI*n / w*cos((M_PI*n / w)*x);
}

double YH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//LSM Symmetric Eigenfunction (y-component)
	double YHs = 0.0;
	if (y >= -c && y < -b)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YHs;
}
double YH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//LSE Asymmetric Eigenfunction (y-component)
	double YHa = 0.0;
	if (y >= -c && y < -b)
	{
		YHa = -1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YHa;
}

double YHc_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Same as YH_S currently, placeholder for complex conjugate
	double YHs = 0.0;
	if (y >= -c && y < -b)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHs = 1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YHs;
}
double YHc_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Same as YH_A currently, placeholder for complex conjugate
	double YHa = 0.0;
	if (y >= -c && y < -b)
	{
		YHa = -1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*sin(kmd*(c - y));
	}
	return YHa;
}

double dyYHc_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Derivative of LSM Symmetric Eigenfunction (y-component), placeholder for complex conjugate
	double YHs = 0.0;
	if (y >= -c && y < -b)
	{
		YHs = kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHs = kch*1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHs = kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHs;
}
double dyYHc_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{ //Derivative of LSM Asymmetric Eigenfunction (y-component), placeholder for complex conjugate
	double YHa = 0.0;

	if (y >= -c && y < -b)
	{
		YHa = -kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHa = kch*1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHa = kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHa;
}
double dyYH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Derivative of LSM Symmetric Eigenfunction (y-component)
	double YHs = 0.0;
	if (y >= -c && y < -b)
	{
		YHs = kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHs = kch*1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHs = kmd*1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHs;
}
double dyYH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Derivative of LSM Asymmetric Eigenfunction (y-component)
	double YHa = 0.0;

	if (y >= -c && y < -b)
	{
		YHa = -kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHa = kch*1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHa = kmd*1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHa;
}

double IYH_S(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Antiderivative of LSM Symmetric Eigenfunction (y-component)
	double YHs = 0.0;
	if (y >= -c && y < -b)
	{
		YHs = -1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHs = 1.0 / 2.0*(1.0 - exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHs = -1.0 / 2.0*(1.0 + exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHs;
}
double IYH_A(double y, double lambH, int n, double Ep, double Mu, double b, double c, double w, double B2, double kch, double kmd)
{//Antiderivative of LSM Asymmetric Eigenfunction (y-component)
	double YHa = 0.0;
	if (y >= -c && y < -b)
	{
		YHa = 1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c + y));
	}
	if (y >= -b && y <= b)
	{
		YHa = 1.0 / 2.0*(1.0 + exp(-2.0*kch*y));
	}
	if (y > b && y <= c)
	{
		YHa = -1.0 / 2.0*(1.0 - exp(-2.0*kch*b)) / (Mu*sin(kmd*(c - b)))*cos(kmd*(c - y));
	}
	return YHa;
}

double XH(double x, int n, double w)
{//x-dependence of LSM eigenfunctions
	return cos(M_PI*n / w*x);
}

// --------------- Wavenumber ------------//

double ky(double y, double lamb, double b, double c, double w, double Ep, double Mu, double B2, int n)
{//calculation of frequency inside and outside vaccum
	double Ky = 0;
	double kx2 = M_PI*n / w*M_PI*n / w;
	if ((y < -b && y > -c) || (y > b && y < c))
		Ky = sqrt((Ep*Mu*B2 - 1)*lamb - kx2);
	if (y >= -b && y <= b)
		Ky = sqrt((1 - B2)*lamb + kx2);
	return Ky;
}

// -------------------------------------------------------- //
// --------------  Find Eigenfunction roots --------------- //
// -------------------------------------------------------- //

void FindRoots(vector<vector<double> >& zeroes, double(*f)(double, double, double, double, double, double, double, double), int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w, double acc)
{
	//resize eigenvalue matrix
	zeroes.resize(sI);
	for (int i = 0; i < sI; ++i)
	{
		zeroes[i].resize(sN);
	}
	double n = 1;                                        //Iteration over parameter n
	double x1 = 1.0;                                     //Lower x guess
	double x2 = 2.0;                                     //Larger x guess
	double x3 = 0.0;                                     //Current x guess
	double xlast = 0.0;                                  //Previous x guess
	double PN = 0.0;                                     //Test for (P)ositive or (N)egative value of function
	double PNlast = 0.0;                                 //Hold previous value of function
	double step = 0.0;                                   //Step size for guessing next root
	int newRoot = 0;                                     //Flag to reset root finder
	for (int nA = 1; nA< sN + 1; nA++)
	{
		n = double(nA);
		x1 = step;
		x2 = x1;
		x3 = x1;
		PN = f(x1, n, Ep, Mu, B2, b, c, w);
		PNlast = f(x1, n, Ep, Mu, B2, b, c, w);
		step = .1;
		xlast = x1;
		int j = 0;
		newRoot = 1;
		//Begin rootfinding search
		while (j < sI)
		{
			if ((PN < 0 && PNlast > 0) || (PN > 0 && PNlast < 0))
			{//If a sign change was detected we have found an apprpriate bracket for root finding
				if (PN < 0)
				{
					x1 = xlast;
					x2 = x3;
				}
				else
				{
					x1 = x3;
					x2 = xlast;
				}
				while (abs(PN) > acc || PN != PN)
				{//Use bisection method to find root within requested accuracy                    
					if (PN > 0)
					{
						x1 = x3;
					}
					else
					{
						x2 = x3;
					}
					x3 = (x1 + x2) / 2.0;
					PN = f(x3, n, Ep, Mu, B2, b, c, w);
				}
				//Store root and reset root finder  near next root
				zeroes[j][nA - 1] = x3;
				j++;
				x3 += M_PI;
				newRoot = 0;
			}
			xlast = x3;
			if (newRoot == 0)
			{//Decide which direction to search for new root
				PN = f(x3, n, Ep, Mu, B2, b, c, w);
				if ((PN < 0 && f(x3 + .1, n, Ep, Mu, B2, b, c, w) > PN) || (PN > 0 && f(x3 + .1, n, Ep, Mu, B2, b, c, w) < PN))
					step = .1;
				else
					step = -.1;
				newRoot = 1;
			}
			x3 += step;
			PNlast = PN;
			PN = f(x3, n, Ep, Mu, B2, b, c, w);
		}
	}
}

// ----------------------------------------------------------------- //
// ---------- Combine for Complete Eigenvalue Calculations --------- //
// ----------------------------------------------------------------- //

std::tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> EigenvaluesCalculator(int sI, int sN, vector<vector<double>>& zeroesES, vector<vector<double>>& zeroesEA, vector<vector<double>>& zeroesHS, vector<vector<double>>& zeroesHA, double Ep, double Mu, double B2, double b, double c, double w, double acc)
{
	//Create and resize vectors to hold eigenvalues of allowed frequency modes
	zeroesES.resize(sI);
	zeroesHS.resize(sI);
	zeroesEA.resize(sI);
	zeroesHA.resize(sI);
	for (int i = 0; i < sI; ++i)
	{
		zeroesES[i].resize(sN);
		zeroesEA[i].resize(sN);
		zeroesHS[i].resize(sN);
		zeroesHA[i].resize(sN);
	}
	//Use dispersion relation to find scaled eigenvalues for each mode
	FindRoots(zeroesES, &f1, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSE Symmetric Modes
	FindRoots(zeroesEA, &f2, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSE Asymmetric Modes
	FindRoots(zeroesHS, &f3, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSM Symmetric Modes
	FindRoots(zeroesHA, &f4, sN, sI, Ep, Mu, B2, b, c, w, acc);  //LSM Asymmetric Modes

																 //Rescale roots to eigenvalues of each mode
	for (int N = 0; N<sN; N++)
	{
		for (int I = 0; I< sI; I++)
		{
			zeroesES[I][N] = ((zeroesES[I][N] / (c - b)*zeroesES[I][N] / (c - b) + (N + 1)*M_PI / w*(N + 1)*M_PI / w) / (Ep*Mu*B2 - 1.0));
			zeroesEA[I][N] = ((zeroesEA[I][N] / (c - b)*zeroesEA[I][N] / (c - b) + (N + 1)*M_PI / w*(N + 1)*M_PI / w) / (Ep*Mu*B2 - 1.0));
			zeroesHS[I][N] = ((zeroesHS[I][N] / (c - b)*zeroesHS[I][N] / (c - b) + (N + 1)*M_PI / w*(N + 1)*M_PI / w) / (Ep*Mu*B2 - 1.0));
			zeroesHA[I][N] = ((zeroesHA[I][N] / (c - b)*zeroesHA[I][N] / (c - b) + (N + 1)*M_PI / w*(N + 1)*M_PI / w) / (Ep*Mu*B2 - 1.0));
		}
	}
	return std::make_tuple(zeroesES, zeroesEA, zeroesHS, zeroesHA);
}

// ----------------------------------------------------------- //
// -------------- Green's Function Calculator ---------------- //
// ----------------------------------------------------------- //

std::tuple<double,double,double> CalcWakeElement(double& Fz, double& Fx, double& Fy, double xi, double x0, double yi, double y0, double zi, double z0, vector<vector<double> > zeroesES, vector<vector<double> > zeroesEA, vector<vector<double> > zeroesHS, vector<vector<double> > zeroesHA, int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w)
{// Calculate force at each desired spot in space
 //Defining additional vectors to help split up computation
 // (a,b,c) component depends on value of (x,y,z)
	double aFz;
	double aFx;
	double aFy;

	double bFzES = 0;
	double bFzEA = 0;
	double bFzHS = 0;
	double bFzHA = 0;
	double bFxES = 0;
	double bFxEA = 0;
	double bFxHS = 0;
	double bFxHA = 0;
	double bFyES = 0;
	double bFyEA = 0;
	double bFyHS = 0;
	double bFyHA = 0;
	double cFzES = 0;
	double cFzEA = 0;
	double cFzHS = 0;
	double cFzHA = 0;
	double cFxES = 0;
	double cFxEA = 0;
	double cFxHS = 0;
	double cFxHA = 0;
	double cFyES = 0;
	double cFyEA = 0;
	double cFyHS = 0;
	double cFyHA = 0;

	//Create Mesh in x
	double x = xi;

	//Create Mesh in y
	double y = yi;

	//Create Mesh in zeta
	double zeta = zi - z0;

	//Loop over x modes
	for (int z = 1; z < sN + 1; z++)
	{
		double kx2 = M_PI*z / w*M_PI*z / w;
		//Calculate force components only dependent on x
		aFz = XE(x0, z, w)*XE(x, z, w);
		aFx = XE(x0, z, w)*dxXE(x, z, w);
		aFy = aFz;

		//Loop over y modes for each x mode
		for (int a = 0; a<sI; a++)
		{
			//Calculate frequency for given eigenvalue
			double kchES = sqrt((1.0 - B2)*zeroesES[a][z - 1] + kx2);
			double kchEA = sqrt((1.0 - B2)*zeroesEA[a][z - 1] + kx2);
			double kchHS = sqrt((1.0 - B2)*zeroesHS[a][z - 1] + kx2);
			double kchHA = sqrt((1.0 - B2)*zeroesHA[a][z - 1] + kx2);
			double kmdES = sqrt((Ep*Mu*B2 - 1.0)*zeroesES[a][z - 1] - kx2);
			double kmdEA = sqrt((Ep*Mu*B2 - 1.0)*zeroesEA[a][z - 1] - kx2);
			double kmdHS = sqrt((Ep*Mu*B2 - 1.0)*zeroesHS[a][z - 1] - kx2);
			double kmdHA = sqrt((Ep*Mu*B2 - 1.0)*zeroesHA[a][z - 1] - kx2);

			//Calculate eigenfunction normalization factors
			double x = kchES*b;
			double As = (2.0*x / (-w / 2.0*(Ep*(1.0 - Ep*Mu*B2)*1.0 / 4.0*(1.0 + 2.0*exp(-2.0*x) + exp(-4.0*x))*(1.0 / (Ep*cos(kmdES*(c - b))))*(1.0 / (Ep*cos(kmdES*(c - b))))*(c - b + sin(2.0*kmdES*(c - b)) / (2.0*kmdES))*2.0*x + (1.0 - B2)*b*(1.0 / 2.0*(1.0 - exp(-4.0*x)) + 2.0*x*exp(-2.0*x)))));
			double Aas = (1.0 / (-w / 2.0*(Ep*(1.0 - Ep*Mu*B2)*1.0 / 4.0*(1.0 - 2.0*exp(-2.0*kchEA*b) + exp(-4.0*b*kchEA)) / (Ep*cos(kmdEA*(c - b))) / (Ep*cos(kmdEA*(c - b)))*(c - b + sin(2.0*kmdEA*(c - b)) / (2.0*kmdEA)) + (1.0 - B2)*(1.0 / 2.0*(1.0 - exp(-4.0*kchEA*b)) / (2.0*kchEA) - b*exp(-2.0*b*kchEA)))));
			double Bs = (1.0 / (-w / 2.0*(Mu*(1.0 - Ep*Mu*B2)*1.0 / 4.0*(1.0 + 2.0*exp(-2.0*kchHS*b) + exp(-4.0*b*kchHS)) / (Mu*sin(kmdHS*(c - b))) / (Mu*sin(kmdHS*(c - b)))*(c - b - sin(2.0*kmdHS*(c - b)) / (2.0*kmdHS)) + (1.0 - B2)*(1.0 / 2.0*(1.0 - exp(-4.0*kchHS*b)) / (2.0*kchHS) + b*exp(-2.0*b*kchHS)))));
			double Bas = (1.0 / (-w / 2.0*(Mu*(1.0 - Ep*Mu*B2)*1.0 / 4.0*(1.0 - 2.0*exp(-2.0*kchHA*b) + exp(-4.0*b*kchHA)) / (Mu*sin(kmdHA*(c - b))) / (Mu*sin(kmdHA*(c - b)))*(c - b - sin(2.0*kmdHA*(c - b)) / (2.0*kmdHA)) + (1.0 - B2)*(1.0 / 2.0*(1.0 - exp(-4.0*kchHA*b)) / (2.0*kchHA) - b*exp(-2.0*b*kchHA)))));

			//FIX THIS 
			//Calculate force component only dependent on zeta and frequency
			bFzES = cos(sqrt(zeroesES[a][z - 1])*abs(zeta));
			bFxES = sin(sqrt(zeroesES[a][z - 1])*abs(zeta));
			bFyES = bFxES;
			bFzEA = cos(sqrt(zeroesEA[a][z - 1])*abs(zeta));
			bFxEA = sin(sqrt(zeroesEA[a][z - 1])*abs(zeta));
			bFyEA = bFxEA;
			bFzHS = cos(sqrt(zeroesHS[a][z - 1])*abs(zeta));
			bFxHS = sin(sqrt(zeroesHS[a][z - 1])*abs(zeta));
			bFyHS = bFxHS;
			bFzHA = cos(sqrt(zeroesHA[a][z - 1])*abs(zeta));
			bFxHA = sin(sqrt(zeroesHA[a][z - 1])*abs(zeta));
			bFyHA = bFxHA;

			//Calculate force component dependent on y
			//INSIDE VACCUM
			double ZeroES = zeroesES[a][z - 1];
			double ZeroEA = zeroesEA[a][z - 1];
			double ZeroHS = zeroesHS[a][z - 1];
			double ZeroHA = zeroesHA[a][z - 1];
			if (abs(y)< b && abs(y0)<b)
			{
				//ELECTRIC COMPONENTS
				//SYMETRIC
				cFzES = exp(kchES*(y0 + y - 2.0*b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*IYE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2);
				cFxES = exp(kchES*(y0 + y - 2.0*b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*IYE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2) / sqrt(ZeroES);
				cFyES = exp(kchES*(y0 + y - 2.0*b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*YE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2) / sqrt(ZeroES);

				//ASYMETRIC
				cFzEA = exp(kchEA*(y0 + y - 2.0*b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*IYE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2);
				cFxEA = exp(kchEA*(y0 + y - 2.0*b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*IYE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2) / sqrt(ZeroEA);
				cFyEA = exp(kchEA*(y0 + y - 2.0*b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*YE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2) / sqrt(ZeroEA);

				//MAGNETIC COMPONENTS
				//SYMETRIC
				cFzHS = exp(kchHS*(y0 + y - 2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*YH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2);
				cFxHS = exp(kchHS*(y0 + y - 2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*YH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2) / sqrt(ZeroHS);
				cFyHS = exp(kchHS*(y0 + y - 2.0*b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*dyYH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2) / sqrt(ZeroHS);
				//ASYMETRIC
				cFzHA = exp(kchHA*(y0 + y - 2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*YH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2);
				cFxHA = exp(kchHA*(y0 + y - 2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*YH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2) / sqrt(ZeroHA);
				cFyHA = exp(kchHA*(y0 + y - 2.0*b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*dyYH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2) / sqrt(ZeroHA);
			}
			else
				if (abs(y)> b && abs(y0)<b)//Inside Dielectric
				{
					//ELECTRIC COMPONENTS
					//SYMETRIC
					cFzES = exp(kchES*(y0 - b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*IYE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2);
					cFxES = exp(kchES*(y0 - b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*IYE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2) / sqrt(ZeroES);
					cFyES = exp(kchES*(y0 - b))*As*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*ky(y, ZeroES, b, c, w, Ep, Mu, B2, z)*(dyYEc_S(y0, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)*YE_S(y, ZeroES, z, Ep, Mu, b, c, w, B2, kchES, kmdES)) / (ZeroES + kx2) / sqrt(ZeroES);

					//ASYMETRIC
					cFzEA = exp(kchEA*(y0 - b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*IYE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2);
					cFxEA = exp(kchEA*(y0 - b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*IYE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2) / sqrt(ZeroEA);
					cFyEA = exp(kchEA*(y0 - b))*Aas*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*ky(y, ZeroEA, b, c, w, Ep, Mu, B2, z)*(dyYEc_A(y0, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)*YE_A(y, ZeroEA, z, Ep, Mu, b, c, w, B2, kchEA, kmdEA)) / (ZeroEA + kx2) / sqrt(ZeroEA);

					//MAGNETIC COMPONENTS
					//SYMETRIC
					cFzHS = exp(kchHS*(y0 - b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*YH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2);
					cFxHS = exp(kchHS*(y0 - b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*YH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2) / sqrt(ZeroHS);
					cFyHS = exp(kchHS*(y0 - b))*Bs*Mu*B2*kx2*YHc_S(y0, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS)*dyYH_S(y, ZeroHS, z, Ep, Mu, b, c, w, B2, kchHS, kmdHS) / (ZeroHS + kx2) / sqrt(ZeroHS);
					//ASYMETRIC
					cFzHA = exp(kchHA*(y0 - b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*YH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2);
					cFxHA = exp(kchHA*(y0 - b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*YH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2) / sqrt(ZeroHA);
					cFyHA = exp(kchHA*(y0 - b))*Bas*Mu*B2*kx2*YHc_A(y0, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA)*dyYH_A(y, ZeroHA, z, Ep, Mu, b, c, w, B2, kchHA, kmdHA) / (ZeroHA + kx2) / sqrt(ZeroHA);
				}
			//Multiply to get total force at each point
			//The minus sign is because it's an electron bunch
			Fx += 4.0*M_PI*aFx*(bFxES*cFxES + bFxEA*cFxEA + bFxHS*cFxHS + bFxHA*cFxHA);
			Fy += 4.0*M_PI*aFy*(bFyES*cFyES + bFyEA*cFyEA + bFyHS*cFyHS + bFyHA*cFyHA);
			Fz += -4.0*M_PI*aFz*(bFzES*cFzES + bFzEA*cFzEA + bFzHS*cFzHS + bFzHA*cFzHA);
		}
	}
	return std::make_tuple(Fx,Fy,Fz);
}

// --------------------------------------------------------------------- //
// ------- Calculate Fields from Green's Function Convolution ---------- //
// --------------------------------------------------------------------- //

std::tuple<vector<double>, vector<double>, vector<double>> TotalForceMeshHDF5(vector <double>& Fz, vector <double>& Fx, vector <double>& Fy, int NParticle, double xi, double x0, double yi, double y0, double zi, double zetamin, double zetamax, double z0, vector<double>& x, vector<double>& y, vector<double>& zeta, vector<vector<double> >& zeroesES, vector<vector<double> >& zeroesEA, vector<vector<double> >& zeroesHS, vector<vector<double> >& zeroesHA, int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w, vector<double> q)
{
	//Want a conversion factor from cgs to SI at the end of this
	//Define at the front rather than inside a loop to save time
	//Times the force by a constant to convert from cgs to SI units
	//Remember we've used coulombs for force so 3e9 is C to esu and then second factor is statvolt/cm to V/m
	//We've been using all the densities in cm^3 so theres a 10^6 in there too (and so on in each dimension) -- simplify with 10^(2*dimensions)
	//This is from the integral calculation (i.e. the first 3 terms are the cgs to SI conversion and the last is the convolution conversion from cm to m)
	double k = 2.99792458e9 * 1e-4 * (1 / 2.99792458) * pow(1e2, 3) * 1e3;

	//Using position vectors, find the repeating values
	//Use this to calulcate the Green's functions we will need (using x1-x2 etc.)
	//For each mesh point, iterate over the distances between points to get the force values

	//Find repeated x, y, and z
	//Sort the vector, find the repeated values and erase the rest
	vector<double> xValues = x;
	std::sort(xValues.begin(), xValues.end());
	auto lastx = std::unique(xValues.begin(), xValues.end());
	xValues.erase(lastx, xValues.end());

	vector<double> yValues = y;
	std::sort(yValues.begin(), yValues.end());
	auto lasty = std::unique(yValues.begin(), yValues.end());
	yValues.erase(lasty, yValues.end());

	vector<double> zValues = zeta;
	std::sort(zValues.begin(), zValues.end());
	auto lastz = std::unique(zValues.begin(), zValues.end());
	zValues.erase(lastz, zValues.end());

	//Set up a vector to store the calculations
	//We'll reserve the size to save time re-constructing the vector
	//For each z row there are xValues.size() * yValues.size() points. Each have xValues.size() * yValues.size() points to calculate and there are zValues.size() rows
	//So size to reserve = zValues.size() * (xValues.size() * yValues.size())^2

	int xSize = (int)xValues.size();
	int ySize = (int)yValues.size();
	int TransverseGridSize = xSize*xSize*ySize*ySize;
	int NumberOfCalculations = (int)((zValues.size()) * TransverseGridSize);

	//Need both the field position and corresponding greens point transversely -> 2D vector
	vector<vector<double>> CalculationsList(NumberOfCalculations, vector<double>(3,0));
	vector<vector<double>> xDistanceList(NumberOfCalculations, vector<double>(2, 0));
	vector<vector<double>> yDistanceList(NumberOfCalculations, vector<double>(2, 0));
	vector<double> zDistanceList(NumberOfCalculations, 0);;
	/*
	CalculationsList.reserve(NumberOfCalculations);
	xDistanceList.reserve(NumberOfCalculations);
	yDistanceList.reserve(NumberOfCalculations);
	zDistanceList.reserve(NumberOfCalculations);
	*/
	//Now do the calculations - working row by row
	std::cout << "Greens Function Calculations:" << std::endl;
	auto GreensStartTime = std::chrono::high_resolution_clock::now();
	double zInterval = zValues[1] - zValues[0];
	double yInterval = yValues[1] - yValues[0];
	double xInterval = xValues[1] - xValues[0];

	for (size_t i = 1; i<zValues.size(); i++) {
		//Distance between the rows
		double zDistance = zValues[i] - zValues[0];
		//For every point in the row
		for (size_t j = 0; j<yValues.size(); j++) {
			double yCalcPoint = yValues[j];
			for (size_t k = 0; k<xValues.size(); k++) {
				double xCalcPoint = xValues[k];
				for (size_t l = 0; l<yValues.size(); l++)
				{
					double yGreen = yValues[l];
					for (size_t m = 0; m<xValues.size(); m++)
					{
						double xGreen = xValues[m];
						int IndexPosition = (i * TransverseGridSize) + (j * (ySize * xSize * xSize)) + (k * (ySize * xSize)) + l * xSize + m;
						//Calculate the Green's function for every point in the row behind
						std::tuple<double,double,double> potential = {0,0,0};
						potential = CalcWakeElement(std::get<0>(potential), std::get<1>(potential), std::get<2>(potential), xCalcPoint, xGreen, yCalcPoint, yGreen, zValues[i], zValues[0], zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w);
						vector<double> Calculations = { std::get<0>(potential), std::get<1>(potential), std::get<2>(potential)};
						/*
						CalculationsList.push_back(Calculations);
						xDistanceList.push_back({ xCalcPoint,xGreen });
						yDistanceList.push_back({ yCalcPoint,yGreen });
						zDistanceList.push_back(std::round(zDistance / zInterval));
						if (isnan(deltaFy) || isnan(deltaFx) || isnan(deltaFz)) {
							std::cout << "NaN Calculated. Try choosing ConvergenceCalculate=false" << std::endl;
							exit(-1);
						}
						*/
						CalculationsList[IndexPosition] = Calculations;
						xDistanceList[IndexPosition] = { xCalcPoint,xGreen };
						yDistanceList[IndexPosition] = { yCalcPoint,yGreen };
						zDistanceList[IndexPosition] = std::round(zDistance / zInterval);
					}
				}
			}
		}
		double fractiondone = (double)i / (zValues.size());
		std::cout << 100 * fractiondone << "%, ";
		auto TimeNow = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> TimeLeft = (TimeNow - GreensStartTime)*(1 / fractiondone - 1);
		std::cout << "Time Remaining: " << TimeLeft.count() << "s" << std::endl;
	}
	auto GreensEndTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> GreensTime = GreensEndTime - GreensStartTime;
	std::cout << "Time for Green's calculations:" << GreensTime.count() << "s" << std::endl;

	//Now go point by point in the mesh - find the distance between mesh points that matches and times by charge
	std::cout << "Force Calculation:" << std::endl;
	auto CalculationTimeStart = std::chrono::high_resolution_clock::now();
	//We should count the number of errors
	int ErrorCount(0);
	//std::cout<<TransverseGridSize <<std::endl;

	for (int f = 0; f<NParticle; f++)
	{
		//Set up the positions we're calculating the force at
		double xP = x[f];
		double yP = y[f];
		double zP = zeta[f];
		//Set up the forces for the point
		double FxSum(0);
		double FySum(0);
		double FzSum(0);

		//Now sum over the calculations from CalcForceMesh, weighted with the density value
		//Number of Particles used for the convolution (takes forever if nP = nParticles)
		for (int iLoop = 0; iLoop<NParticle; iLoop++)
		{
			//Give the position of the small element of the bunch and initialise a matrix for the small force element
			double zB = zeta[iLoop];
			double xB = x[iLoop];
			double yB = y[iLoop];

			double zDifference = std::round((zP - zB) / zInterval);

			if (zDifference>0 && q[iLoop]>0)
			{
				double deltaFz = 0;
				double deltaFy = 0;
				double deltaFx = 0;

				//Use the CalcForceMesh with no x,y,zdiv to get the single point
				//This is the Green's function calculation                   
				//Check if (xP-xB, yP-yB, zP-zB) has been calculated                   
				//Find the positions of the mesh points - take advantage of the ordered nature of the Green's calculations
				//Going inwards we find yP, then xP, then yB, then xB
				int zPosition = zDifference - 1;
				int yPPosition = std::round((yP - yValues[0]) / yInterval);
				int xPPosition = std::round((xP - xValues[0]) / xInterval);
				int yBPosition = std::round((yB - yValues[0]) / yInterval);
				int xBPosition = std::round((xB - xValues[0]) / xInterval);

				int IndexPosition = (zPosition*TransverseGridSize) + (yPPosition*(ySize*xSize*xSize)) + (xPPosition*(ySize*xSize)) + yBPosition*xSize + xBPosition;

				if ((xDistanceList[IndexPosition][0] == xP) && (xDistanceList[IndexPosition][1] == xB) && (yDistanceList[IndexPosition][0] == yP) && (yDistanceList[IndexPosition][1] == yB)) {
					deltaFx = CalculationsList[IndexPosition][0];
					deltaFy = CalculationsList[IndexPosition][1];
					deltaFz = CalculationsList[IndexPosition][2];
				}
				else {
					ErrorCount++;
				}
				FzSum += (q[iLoop] * deltaFz);
				FxSum += (q[iLoop] * deltaFx);
				FySum += (q[iLoop] * deltaFy);
			}
		}
		//Times each by the conversion factor to SI
		Fx[f] = k*FxSum;
		Fy[f] = k*FySum;
		Fz[f] = k*FzSum;
		if (f % 50 == 0) {
			double fractiondone = (double)f / y.size();
			std::cout << 100 * fractiondone << "%" << std::endl;
		}

	}
	auto CalculationEndTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> CalculationTime = CalculationEndTime - CalculationTimeStart;
	std::cout << "Force Calculation time: " << CalculationTime.count() << "s" << std::endl;
	std::cout << "Calculation Mesh Point Errors: " << ErrorCount << std::endl;
	return std::make_tuple(Fx,Fy,Fz);
}

std::tuple<vector<double>, vector<double>, vector<double>> TotalForceMeshHDF5VerticalPlate(vector <double>& Fz, vector <double>& Fx, vector <double>& Fy, int NParticle, double xi, double x0, double yi, double y0, double zi, double zetamin, double zetamax, double z0, vector<double>& x, vector<double>& y, vector<double>& zeta, vector<vector<double> >& zeroesES, vector<vector<double> >& zeroesEA, vector<vector<double> >& zeroesHS, vector<vector<double> >& zeroesHA, int sN, int sI, double Ep, double Mu, double B2, double b, double c, double w, vector<double> q)
{
	std::tuple<vector<double>, vector<double>, vector<double>> fields = TotalForceMeshHDF5(Fz, Fy, Fx, NParticle, xi, x0, yi, y0, zi, zetamin, zetamax, z0, y, x, zeta, zeroesES, zeroesEA, zeroesHS, zeroesHA, sN, sI, Ep, Mu, B2, b, c, w, q);
	return fields;
}


#define Planar_GF

#endif /* Planar_GF */