///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for calculating the Green's function in circular DLW
// Using the method from Ng Phys Rev D and based on DiWakeCy (https://github.com/NIUaard/DiWakeCyl/tree/master)
// Conversion to cylindrical coordinates for use with DiWaCAT fields
// T Overton 06/23
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef Circ_GF
#define _USE_MATH_DEFINES
#include <ctime>
#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>
#include <fstream>
#include <numeric>

using std::vector;

const double epsilon0 = 8.854187817620e-12;
const double qelec = 1.601e-19;
const double cms = 299792458;
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

// ---------------------------------------------- //
//---- Converting Mesh Points to Cylindrical ---- //
// ---------------------------------------------- //

void CartMeshToCyl(vector<double>& rMesh, vector <double>& thetaMesh, vector<double> xvalues, vector<double> yvalues) {
	//Make sure the radial mesh is the right size
	rMesh.resize(xvalues.size());
	thetaMesh.resize(xvalues.size());
	for (size_t i = 0; i < xvalues.size(); i++) {
		double x = xvalues[i];
		double y = yvalues[i];
		if (x == 0 && y == 0) {
			rMesh[i] = 0;
			thetaMesh[i] = 0;
		}
		else if(y == 0) {
			rMesh[i] = abs(x);
			if (x > 0) {
				thetaMesh[i] = M_PI * 0.5;
			}
			else {
				thetaMesh[i] = M_PI * 1.5;
			}
		}
		else if (x == 0) {
			rMesh[i] = abs(y);
			if (y > 0) {
				thetaMesh[i] = 0;
			}
			else {
				thetaMesh[i] = M_PI;
			}
		}
		else {
			rMesh[i] = abs(sqrt(x*x + y*y));
			double theta = atan(abs(x / y));
			if (x > 0 && y>0) {
				thetaMesh[i] = theta;
			}
			else if (x < 0 && y>0) {
				thetaMesh[i] = 2 * M_PI - theta;
			}
			else if (x > 0 && y<0) {
				thetaMesh[i] = M_PI - theta;
			}
			else if (x < 0 && y<0) {
				thetaMesh[i] = M_PI + theta;
			}
		}
	}
}

void CylForceToCart(vector<double> &Fx, vector<double> &Fy, vector<double> Fr, vector<double> Ftheta, vector<double> theta) {
	if (Fr.size() != theta.size()) {
		std::cout << "Force and angular position array must be of equal length.";
		exit(-1);
	}
	Fx.resize(Fr.size());
	Fy.resize(Fr.size());
	for (size_t i = 0; i < theta.size(); i++) {
		Fy[i] = (Fr[i] * cos(theta[i])) - (Ftheta[i] * sin(theta[i]));
		Fx[i] = (Fr[i] * sin(theta[i])) - (Ftheta[i] * cos(theta[i]));
	}
}


//------------------------------------------------//
// ------- Functions for Dispersion Eqn --------- //
// ---------------------------------------------- //

//---------------Bessel Functions -------------//
// Taken from https://www.atnf.csiro.au/computing/software/gipsy/sub/bessel.c
static double bessj0(double x)
// Evaluate Bessel function of first kind and order 0 at input x
{
	double ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7
			+ y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))));
		ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718
			+ y*(59272.64853 + y*(267.8532712 + y*1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z*z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4
			+ y*(-0.2073370639e-5 + y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y*(0.1430488765e-3
			+ y*(-0.6911147651e-5 + y*(0.7621095161e-6
				- y*0.934935152e-7)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
	}
	return ans;
}
static double bessj1(double x)
// Evaluate Bessel function of first kind and order 1 at input x 
{
	double ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1
			+ y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
		ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74
			+ y*(99447.43394 + y*(376.9991397 + y*1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z*z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4
			+ y*(0.2457520174e-5 + y*(-0.240337019e-6))));
		ans2 = 0.04687499995 + y*(-0.2002690873e-3
			+ y*(0.8449199096e-5 + y*(-0.88228987e-6
				+ y*0.105787412e-6)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}
double cyl_bessel_j(int n, double x)
// Evaluate Bessel function of first kind and order n at input x. The function can also be called for n = 0 and n = 1
{
	int    j, jsum, m;
	double ax, bj, bjm, bjp, sum, tox, ans;
	ax = fabs(x);
	if (n == 0)
		return(bessj0(ax));
	if (n == 1)
		return(bessj1(ax));
	if (ax == 0.0)
		return 0.0;
	else if (ax > (double)n) {
		tox = 2.0 / ax;
		bjm = bessj0(ax);
		bj = bessj1(ax);
		for (j = 1; j<n; j++) {
			bjp = j*tox*bj - bjm;
			bjm = bj;
			bj = bjp;
		}
		ans = bj;
	}
	else {
		tox = 2.0 / ax;
		m = 2 * ((n + (int)sqrt(ACC*n)) / 2);
		jsum = 0;
		bjp = ans = sum = 0.0;
		bj = 1.0;
		for (j = m; j>0; j--) {
			bjm = j*tox*bj - bjp;
			bjp = bj;
			bj = bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum = !jsum;
			if (j == n) ans = bjp;
		}
		sum = 2.0*sum - bj;
		ans /= sum;
	}
	return  x < 0.0 && n % 2 == 1 ? -ans : ans;
}

static double bessy0(double x)
// Evaluate Bessel function of second kind and order 0 at input x //
{
	double z;
	double xx, y, ans, ans1, ans2;

	if (x < 8.0) {
		y = x*x;
		ans1 = -2957821389.0 + y*(7062834065.0 + y*(-512359803.6
			+ y*(10879881.29 + y*(-86327.92757 + y*228.4622733))));
		ans2 = 40076544269.0 + y*(745249964.8 + y*(7189466.438
			+ y*(47447.26470 + y*(226.1030244 + y*1.0))));
		ans = (ans1 / ans2) + 0.636619772*bessj0(x)*log(x);
	}
	else {
		z = 8.0 / x;
		y = z*z;
		xx = x - 0.785398164;
		ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4
			+ y*(-0.2073370639e-5 + y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y*(0.1430488765e-3
			+ y*(-0.6911147651e-5 + y*(0.7621095161e-6
				+ y*(-0.934945152e-7))));
		ans = sqrt(0.636619772 / x)*(sin(xx)*ans1 + z*cos(xx)*ans2);
	}
	return ans;
}

static double bessy1(double x)
// Evaluate Bessel function of second kind and order 1 at input x .//
{
	double z;
	double xx, y, ans, ans1, ans2;
	if (x < 8.0) {
		y = x*x;
		ans1 = x*(-0.4900604943e13 + y*(0.1275274390e13
			+ y*(-0.5153438139e11 + y*(0.7349264551e9
				+ y*(-0.4237922726e7 + y*0.8511937935e4)))));
		ans2 = 0.2499580570e14 + y*(0.4244419664e12
			+ y*(0.3733650367e10 + y*(0.2245904002e8
				+ y*(0.1020426050e6 + y*(0.3549632885e3 + y)))));
		ans = (ans1 / ans2) + 0.636619772*(bessj1(x)*log(x) - 1.0 / x);
	}
	else {
		z = 8.0 / x;
		y = z*z;
		xx = x - 2.356194491;
		ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4
			+ y*(0.2457520174e-5 + y*(-0.240337019e-6))));
		ans2 = 0.04687499995 + y*(-0.2002690873e-3
			+ y*(0.8449199096e-5 + y*(-0.88228987e-6
				+ y*0.105787412e-6)));
		ans = sqrt(0.636619772 / x)*(sin(xx)*ans1 + z*cos(xx)*ans2);
	}
	return ans;
}
double cyl_neumann(int n, double x)
{
	// Evaluate Bessel function of second kind and order  n for input x. (n >= 0) //
	int j;
	double by, bym, byp, tox;

	if (n == 0)
		return(bessy0(x));
	if (n == 1)
		return(bessy1(x));

	tox = 2.0 / x;
	by = bessy1(x);
	bym = bessy0(x);
	for (j = 1; j<n; j++) {
		byp = j*tox*by - bym;
		bym = by;
		by = byp;
	}
	return by;
}

double bessel_j_prime(int n, double x) {
	if (n == 0) {
		return -1 * cyl_bessel_j(n + 1, x);
	}
	double bessel_below = cyl_bessel_j(n - 1, x);
	double bessel_above = cyl_bessel_j(n + 1, x);
	return 0.5*(bessel_below - bessel_above);
}
double neumann_j_prime(int n, double x) {
	if (n == 0) {
		return -1 * cyl_neumann(n + 1, x);
	}
	double neumann_below = cyl_neumann(n - 1, x);
	double neumann_above = cyl_neumann(n + 1, x);
	return 0.5*(neumann_below - neumann_above);
}

double R_rg(double x, double a, double b, int n) {
	return bessel_j_prime(n,x*a)*cyl_neumann(n,x*b) - neumann_j_prime(n,x*a)*cyl_bessel_j(n,x*b);
}
double Rp_rg(double x, double a, double b, int n) {
	return bessel_j_prime(n, x*a)*neumann_j_prime(n, x*b) - neumann_j_prime(n, x*a)*bessel_j_prime(n, x*b);
}
double P_rg(double x, double a, double b, int n) {
	return cyl_bessel_j(n, x*a)*cyl_neumann(n, x*b) - cyl_bessel_j(n, x*b)*cyl_neumann(n, x*a);
}
double Pp_rg(double x, double a, double b, int n) {
	return cyl_bessel_j(n,x*a)*neumann_j_prime(n, x*b) - bessel_j_prime(n, x*b)*cyl_neumann(n, x*a);
}
int Heaviside(double x) {
	if (x > 0) return 1;
	else return 0;
}

// ---------------------------------------------- //
// -------- Dispersion Equation ----------------- //
// ---------------------------------------------- //

double DispersionEqn(double s, double b, double c, int n, double mu, double epsilon) {
	//Dispersion eqn from Ng's paper - x here correspeonds to s in that paper
	double x = s*c;
	double xi = b / c;
	if (n == 0) {
		double P = P_rg(s, c, b, n);
		double Pp = Pp_rg(s, c, b, n);
		return (x*Pp + (x*x*xi*P)*1/(2*epsilon));
	}
	else {
		double P = P_rg(s, c, b, n);
		double Pp = Pp_rg(s, c, b, n);
		double R = R_rg(s, c, b, n);
		double Rp = Rp_rg(s, c, b, n);
		return ((x*x*xi*xi * 1 / (n + 1) - n*(mu*epsilon + 1))*P*R + x*xi*(epsilon*Pp*R + mu*Rp*P));
	}
}

// ------------------------------------------ //
// ------------- Mode Calculation ----------- //
// ------------------------------------------ //

std::tuple<vector<vector<double>>, vector<vector<double>>> FindModes(double b, double c, double mu, double epsilon, int nTheta, int nRadial, double kStep, bool PrintProgress) {
	//k is a stepping parameter used to scan the dispersion function for a zero
	//Also the accuracy of mode calculations
	auto StartTime = std::chrono::high_resolution_clock::now();
	double kTol = 1e-6;
	vector<vector<double>> WaveVec, WaveAmp;
	WaveVec.resize(nTheta);
	WaveAmp.resize(nTheta);
	for (int i = 0; i < nTheta; ++i)
	{
		WaveVec[i].resize(nRadial);
		WaveAmp[i].resize(nRadial);
	}
	//To reset the k value after each set of modes
	double kInitial = kStep;
	for (int i = 0; i < nTheta; i++) {
		int CurrentRadialMode = 0;
		double k = kInitial;
		vector<double> Roots(nRadial);
		//Solve to find the next root of the dispersion eqn in turn
		while (CurrentRadialMode < nRadial) {
			double kTolerance = kStep;
			double Dkmin = DispersionEqn(k, b, c, i, mu, epsilon);
			double Dkmax = DispersionEqn(k + kStep, b, c, i, mu, epsilon);
			if (Dkmax > 0 && Dkmin < 0) {
				while (kTolerance > kTol) {
					kTolerance *= 0.5;
					Dkmax = DispersionEqn(k + kTolerance, b, c, i, mu, epsilon);
					if (Dkmax < 0) {
						Dkmin = Dkmax;
						Dkmax = DispersionEqn(k + 2 * kTolerance, b, c, i, mu, epsilon);
						k = k + kTolerance;
					}
				}
				if (Dkmax < 0) k = k - kTolerance;
				Roots[CurrentRadialMode] = k;
				//Save the sign of the gradient of dispersion function for the amplitude later
				WaveAmp[i][CurrentRadialMode] = +1;
				CurrentRadialMode++;
			}
			if (Dkmax < 0 && Dkmin > 0) {
				while (kTolerance > kTol) {
					kTolerance *= 0.5;
					Dkmax = DispersionEqn(k + kTolerance, b, c, i, mu, epsilon);
					if (Dkmax > 0) {
						Dkmin = Dkmax;
						Dkmax = DispersionEqn(k + 2 * kTolerance, b, c, i, mu, epsilon);
						k = k + kTolerance;
					}
				}
				if (Dkmax > 0) k = k - kTolerance;
				Roots[CurrentRadialMode] = k;
				//Save the sign of the gradient of dispersion function for the amplitude later
				WaveAmp[i][CurrentRadialMode] = -1;
				CurrentRadialMode++;
			}
			k += kStep;
		}

		//Amplitude of each mode
		double CGStoSI = -1 / (4 * M_PI*epsilon0);
		double FieldToWake = 1 / qelec;
		if (i == 0) {
			double ConstantFactor = 4 * qelec * 1 / (c*b*epsilon);
			for (int j = 0; j < nRadial; j++) {
				double delta = 2* kTol;
				double D_s(0); //Set to a value just so we enter the while loop at least once
				bool DsFound = false;
				while (DsFound == false) {
					D_s = abs((DispersionEqn(Roots[j] + delta, b, c, i, mu, epsilon) - DispersionEqn(Roots[j], b, c, i, mu, epsilon)) * 1 / (c*delta));
					if (isinf(1/D_s) || isnan(1 / D_s)) {
						delta *= 10;
					}
					else DsFound = true;
				}
				WaveAmp[i][j] = WaveAmp[i][j] * CGStoSI * ConstantFactor * FieldToWake * Roots[j]* c* P_rg(Roots[j], b, c, i) * 1 / D_s;
				WaveVec[i][j] = Roots[j] * 1 / (sqrt(epsilon*mu - 1));
			}
		}
		if (i > 0) {
			for (int j = 0; j < nRadial; j++) {
				double delta = 2* kTol;
				double D_s = 10; //Set to a value just so we enter the while loop at least once
				bool DsFound = false;
				while (DsFound == false) {
					D_s = abs((DispersionEqn(Roots[j] + delta, b, c, i, mu, epsilon) - DispersionEqn(Roots[j], b, c, i, mu, epsilon)) * 1 / (c*delta));
					if (isinf(1 / D_s) || isnan(1 / D_s)) {
						delta *= 10;
					}
					else DsFound = true;
				}
				WaveAmp[i][j] = WaveAmp[i][j] * CGStoSI * P_rg(Roots[j], b, c, i) * R_rg(Roots[j],b,c,i) * 1 / D_s;
				WaveVec[i][j] = Roots[j] * 1 / (sqrt(epsilon*mu - 1));
			}
		}
		if (PrintProgress == true) {
			auto CurrentTime = std::chrono::high_resolution_clock::now();
			double fractiondone = (double)(i + 1)/nTheta;
			std::chrono::duration<double> TimeLeft = (CurrentTime - StartTime)*(1 / fractiondone);
			std::cout << "Time Remaining: " << TimeLeft.count() << "s" << std::endl;
		}
	}
	return std::make_tuple(WaveVec, WaveAmp);
}

vector<vector<double>> FindModesSingleNT(int nTheta, int nRadial, double b, double c, double mu, double epsilon, double kStep) {
	//Same function as FindModes but only for a single value of nTheta
	vector<double> ModeAmp(nRadial);
	vector<double> Roots(nRadial);
	int CurrentRadialMode = 0;
	double k = kStep;
	double kTol = 1e-6;
	//Solve to find the next root of the dispersion eqn in turn
	while (CurrentRadialMode < nRadial) {
		double kTolerance = kTol;
		double Dkmin = DispersionEqn(k, b, c, nTheta, mu, epsilon);
		double Dkmax = DispersionEqn(k + kStep, b, c, nTheta, mu, epsilon);
		if (Dkmax > 0 && Dkmin < 0) {
			while (kTolerance > 1e-6) {
				kTolerance *= 0.5;
				Dkmax = DispersionEqn(k + kTolerance, b, c, nTheta, mu, epsilon);
				if (Dkmax < 0) {
					Dkmin = Dkmax;
					Dkmax = DispersionEqn(k + 2 * kTolerance, b, c, nTheta, mu, epsilon);
					k = k + kTolerance;
				}
			}
			Roots[CurrentRadialMode] = k;
			ModeAmp[CurrentRadialMode] = +1;
			CurrentRadialMode++;
		}
		if (Dkmax < 0 && Dkmin > 0) {
			while (kTolerance > 1e-6) {
				kTolerance *= 0.5;
				Dkmax = DispersionEqn(k + kTolerance, b, c, nTheta, mu, epsilon);
				if (Dkmax > 0) {
					Dkmin = Dkmax;
					Dkmax = DispersionEqn(k + 2 * kTolerance, b, c, nTheta, mu, epsilon);
					k = k + kTolerance;
				}
			}
			Roots[CurrentRadialMode] = k;
			ModeAmp[CurrentRadialMode] = -1;
			CurrentRadialMode++;
		}
		k += kStep;
	}
	//Amplitude of each mode
	double CGStoSI = -1 / (4 * M_PI*epsilon0);
	double FieldToWake = 1 / qelec;
	if (nTheta == 0) {
		double ConstantFactor = 4 * qelec * 1 / (c*b*epsilon);
		for (int j = 0; j < nRadial; j++) {
			double delta = kTol * 2;
			double D_s(0); //Set to a value just so we enter the while loop at least once
			bool DsFound = false;
			while (DsFound == false) {
				D_s = abs((DispersionEqn(Roots[j] + delta, b, c, nTheta, mu, epsilon) - DispersionEqn(Roots[j], b, c, nTheta, mu, epsilon)) * 1 / (c*delta)); //c here and in mode amp is cgs to si for dispersion eqn
				if (isinf(1 / D_s) || isnan(1 / D_s)) {
					delta *= 0.1;
				}
				else DsFound = true;
			}
			ModeAmp[j] = ModeAmp[j] * CGStoSI * ConstantFactor * FieldToWake * Roots[j] * c * P_rg(Roots[j], b, c, nTheta) * 1 / D_s;
			Roots[j] = Roots[j] * 1 / (sqrt(epsilon*mu - 1));
		}
	}
	if (nTheta > 0) {
		for (int j = 0; j < nRadial; j++) {
			double delta = kTol * 2;
			double D_s(0); //Set to a value just so we enter the while loop at least once
			bool DsFound = false;
			while (DsFound == false) {
				D_s = abs((DispersionEqn(Roots[j] + delta, b, c, nTheta, mu, epsilon) - DispersionEqn(Roots[j], b, c, nTheta, mu, epsilon)) * 1 / (c*delta));
				if (isinf(1 / D_s) || isnan(1 / D_s)) {
					delta *= 0.1;
				}
				else DsFound = true;
			}
			ModeAmp[j] = ModeAmp[j] * CGStoSI  * P_rg(Roots[j], b, c, nTheta) * R_rg(Roots[j], b, c, nTheta) * 1 / D_s;
			Roots[j] = Roots[j] * 1 / (sqrt(epsilon*mu - 1));
		}
	}
	return{ Roots,ModeAmp };
}

std::tuple<vector<vector<double>>, vector<vector<double>>> ModeConvergence(double rMax, double b, double c, double mu, double epsilon, double kStep, double ModeAccuracy) {
	int nR(10), nTheta(2);
	double InitialThetaConvergenceSum = rMax / b + pow(rMax / b, 2);
	while (ModeAccuracy * InitialThetaConvergenceSum < pow(rMax / b, nTheta + 1)) {
		nTheta++;
		InitialThetaConvergenceSum += pow(rMax / b, nTheta);
	}
	std::tuple<vector<vector<double>>, vector<vector<double>>> modes;
	bool ModeConverged = false;
	std::cout << "Calculating Convergence" << std::endl;
	while (ModeConverged == false) {
		bool RModeConverged = false;
		bool TModeConverged = false;
		vector<vector<double>> ModeAmplitude, WaveVectors;
		//FindModes(WaveVectors, ModeAmplitude, b, c, mu, epsilon, nTheta, nR, kStep,false);
		while (RModeConverged == false) {
			vector<double> AmplitudesZero = FindModesSingleNT(0, nR, b, c, mu, epsilon, kStep)[1];
			vector<double> AmplitudesTheta = FindModesSingleNT(nTheta, nR, b, c, mu, epsilon, kStep)[1];
			double nRSumZero = 0;
			double nRSumTheta = 0;
			for (int i = 0; i < nR; i++) {
				nRSumZero += abs(AmplitudesZero[i]);
				nRSumTheta += abs(AmplitudesTheta[i]);
			}
			double RConvergence = (abs(AmplitudesZero[nR - 1])+ abs(AmplitudesZero[nR - 2])) / nRSumZero;
			//RConvergence += (abs(AmplitudesTheta[nR - 1]) + abs(AmplitudesTheta[nR - 2])) / nRSumTheta;
			if (RConvergence < ModeAccuracy) RModeConverged = true;
			if (RModeConverged == false) {
				std::cout << "Radial Convergence = " << 100 * RConvergence << "%" << std::endl;
				nR++;
			}
		}
		modes = FindModes(b, c, mu, epsilon, nTheta, 1, kStep, false);
		double nTSum(0);
		for (int i = 1; i < nTheta; i++) {
			nTSum += abs(std::get<1>(modes)[i][0] * pow(rMax / b, i));
			//std::cout << nTSum << std::endl;
		}
		while (TModeConverged == false) {
			vector<vector<double>> NextModeSolution = FindModesSingleNT(nTheta, 1, b, c, mu, epsilon, kStep);
			nTSum += abs(NextModeSolution[1][0] * pow(rMax / b, nTheta));
			double TConvergence = (NextModeSolution[1][0] * pow(rMax / b, nTheta) * 1 / nTSum);
			if (TConvergence < ModeAccuracy && nTSum > 0) TModeConverged = true;
			if (TModeConverged == false) {
				std::cout << "Azimuthal Convergence = " << 100 * TConvergence << "%" << std::endl;
				nTheta++;
			}
		}
		std::cout << "Calculating Full Mode Composition" << std::endl;
		std::cout << "Number of Radial Modes: " << nR << std::endl;
		std::cout << "Number of Azimuthal Modes: " << nTheta << std::endl;
		modes = FindModes(b, c, mu, epsilon, nTheta, nR, kStep, true);
		if (RModeConverged == true && TModeConverged == true) ModeConverged = true;
	}
	return modes;
}


std::tuple<double, double, double> CalcWakeCyl(double &Fz, double &Fr, double &Ftheta, vector<vector<double>> RootAmp, vector<vector<double>> WaveVec, double r0, double r, double theta0, double theta, double b, double c, double z0, double z, double mu, double epsilon) {
	double zz = z - z0;
	if (zz < 0) {
		return std::make_tuple(Fz,Fr,Ftheta);
	}
	if (r > b || r0 > b) {
		return std::make_tuple(Fz,Fr,Ftheta);
	}
	for (size_t j = 0; j < RootAmp.size(); j++) {
		//Loop through theta modes
		for (size_t i = 0; i < RootAmp[0].size(); i++) {
			//Loop through R modes
			double Wz(0), Wr(0), Wtheta(0);
			if (j == 0) {
				Wz = -1 * RootAmp[j][i] * cos(WaveVec[j][i] * zz);
			}
			else {
				double constLong = 8 * 1 / (c*c) * pow(r0*r * 1 / (b*b), j) * WaveVec[j][i] * (c*sqrt(mu*epsilon - 1)); //c*sqrt() is to convert from frequency to wavenumber
				double constR = 1 / (c*c) * pow(r0 * 1 / b, j) * pow(r * 1 / b, j - 1)* (c/b) * 8 * j * sqrt(mu*epsilon - 1); //Minus sign is because it's a force field -> need to account for electron charge
				double constTheta = -1 / (c*c) * pow(r0 * 1 / b, j) * pow(r * 1 / b, j - 1) * (c / b) * 8 * j * sqrt(mu*epsilon - 1); //Positive sign is because it's a force field -> need to account for electron charge - potential is negative due to cos differential
				Wz = constLong * RootAmp[j][i] * cos(WaveVec[j][i] * zz) * cos(j*(theta0 - theta));
				Wr = constR*RootAmp[j][i] * sin(WaveVec[j][i] * zz) * cos(j*(theta0 - theta));//* 1 / (WaveVec[j][i] * sqrt(epsilon*mu - 1)*c) 
				Wtheta = constTheta*RootAmp[j][i]  * sin(WaveVec[j][i] * zz) * sin(j*(theta0 - theta));//* 1 / (WaveVec[j][i] * sqrt(epsilon*mu - 1)*c)
			}
			Fz += Wz;
			Fr += Wr;
			Ftheta += Wtheta;
		}
	}
	return std::make_tuple(Fz, Fr, Ftheta);
}

std::tuple<vector<double>, vector<double>, vector<double>> TotalForceMeshCircular(vector<double> &Fz, vector<double> &Fr, vector<double> &Ftheta, vector<double> r, vector<double> theta, vector<double> zeta, vector<double> q, vector<vector<double>> ModeAmplitude, vector<vector<double>> WaveNums, double b, double c, double mu, double epsilon) {

	//Using a regular mesh - we want to find the repeating values
	//For transverse positions we need to find the repeating pairs - find the points with the same z value
	//Find repeating z
	//Positions are in cm because of planar functions - so need *0.01 for mesh positions

	int NParticle = r.size();

	vector<double> zValues = zeta;
	std::sort(zValues.begin(), zValues.end());
	auto lastz = std::unique(zValues.begin(), zValues.end());
	zValues.erase(lastz, zValues.end());

	vector<int> TransverseMeshIndices(0);
	for (size_t i = 0; i < zeta.size(); i++) {
		if (zeta[i] == zValues[0]) {
			TransverseMeshIndices.push_back(i);
		}
	}
	vector<double> rValues(0);
	vector<double> thetaValues(0);
	for (size_t i = 0; i < TransverseMeshIndices.size(); i++) {
		rValues.push_back(r[TransverseMeshIndices[i]]*0.01);
		thetaValues.push_back(theta[TransverseMeshIndices[i]]);
	}

	//Set up a vector to store the calculations
	//Reserving the space is more efficient
	int zSize = (int)zValues.size();
	int TransverseSize = (int)rValues.size();
	int NumberOfCalculations = zSize*TransverseSize*TransverseSize;
	vector<vector<double>> CalculationsList;
	vector<vector<double>> rDistanceList;
	vector<vector<double>> thetaDistanceList;
	vector<double> zDistanceList;

	CalculationsList.reserve(NumberOfCalculations);
	rDistanceList.reserve(NumberOfCalculations);
	thetaDistanceList.reserve(NumberOfCalculations);
	zDistanceList.reserve(NumberOfCalculations);

	//Now do the calculations going row by row
	std::cout << "Green's Functions Calculations:" << std::endl;
	auto GreensStartTime = std::chrono::high_resolution_clock::now();
	double zInterval = (zValues[1] - zValues[0]);

	for (size_t i = 1; i < zValues.size(); i++) {
		double zDistance = (zValues[i] - zValues[0])*0.01;
		//For every point in the row	
		for (int j = 0; j < TransverseSize; j++) {
			double rWitness = rValues[j];
			double thetaWitness = thetaValues[j];
			for (int k = 0; k < TransverseSize; k++) {
				double rSource = rValues[k];
				double thetaSource = thetaValues[k];
				//Calculate the Green's function with these source and witness positions
				double deltaFz(0), deltaFr(0), deltaFt(0);
				std::tuple<double, double, double> values = CalcWakeCyl(deltaFz, deltaFr, deltaFt, ModeAmplitude, WaveNums, rSource, rWitness, thetaSource, thetaWitness, b, c, 0, zDistance, mu, epsilon);
				vector<double> Calculations = { std::get<1>(values), std::get<2>(values), std::get<0>(values) };
				CalculationsList.push_back(Calculations);
				rDistanceList.push_back({ rWitness, rSource });
				thetaDistanceList.push_back({ thetaWitness, thetaSource });
				zDistanceList.push_back(zDistance);
				if (isnan(deltaFz) || isnan(deltaFr) || isnan(deltaFt)) {
					std::cout << "NaN Calculated. Try using ConvergenceCalculate=false" << std::endl;
					exit(-1);
				}
			}
		}
		double fractiondone = (double)i / (zValues.size());
		std::cout << 100 * fractiondone << "%. ";
		auto TimeNow = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> TimeLeft = (TimeNow - GreensStartTime)*(1 / fractiondone - 1);
		std::cout << "Time Remaining: " << TimeLeft.count() << "s" << std::endl;
	}
	auto GreensEndTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> GreensTime = GreensEndTime - GreensStartTime;
	std::cout << "Time for Green's Calculations: " << GreensTime.count() << "s" << std::endl;

	//Now go point by point in the mesh and find the corresponding already calculated wake potential to calculate the force
	std::cout << "Force Calculation" << std::endl;
	int ErrorCount(0);
	auto CalculationStart = std::chrono::high_resolution_clock::now();
	for (int f = 0; f < NParticle; f++) {
		//Positions calculating the force at
		double rP = r[f]*0.01;
		double tP = theta[f];
		double zP = zeta[f];
		//Set up the forces
		double FzSum(0), FrSum(0), FtSum(0);

		//Now sum over the calculation points from before
		for (int iLoop = 0; iLoop < NParticle; iLoop++) {
			double rB = r[iLoop]*0.01;
			double tB = theta[iLoop];
			double zB = zeta[iLoop];
			double zDifference = (zP - zB);
			if (zDifference > 0 && abs(q[iLoop]) > 0) {
				double deltaFz(0), deltaFr(0), deltaFt(0);
				int zIndex = std::round(zDifference/zInterval) - 1;
				int WitnessTransverseIndex(-1), SourceTransverseIndex(-1);
				int CalcPointIndex(-1);
				for (int i = 0; i < TransverseSize; i++) {
					if (rP == rValues[i] && tP == thetaValues[i]) {
						WitnessTransverseIndex = i;
					}
				}
				for (int i = 0; i < TransverseSize; i++) {
					if (rB == rValues[i] && tB == thetaValues[i]) {
						SourceTransverseIndex = i;
					}
				}
				if (SourceTransverseIndex != -1 && WitnessTransverseIndex != -1) {
					CalcPointIndex = zIndex*TransverseSize*TransverseSize + WitnessTransverseIndex*TransverseSize + SourceTransverseIndex;
				}
				
				if (CalcPointIndex != -1) {
					if (rDistanceList[CalcPointIndex][0] == rP && rDistanceList[CalcPointIndex][1] == rB && thetaDistanceList[CalcPointIndex][0] == tP && thetaDistanceList[CalcPointIndex][1] == tB) {
						deltaFr = CalculationsList[CalcPointIndex][0];
						deltaFt = CalculationsList[CalcPointIndex][1];
						deltaFz = CalculationsList[CalcPointIndex][2];
					}
					else ErrorCount++;
					if (std::round(zDistanceList[CalcPointIndex]*100/zInterval) != std::round((zDifference)/zInterval)) ErrorCount++;
				}
				else ErrorCount++;
				FzSum += q[iLoop] * deltaFz;
				FrSum += q[iLoop] * deltaFr;
				FtSum += q[iLoop] * deltaFt;
			}
		}
		Fz[f] = FzSum;
		Fr[f] = FrSum;
		Ftheta[f] = FtSum;
		if (f % 100 == 0) {
			double fractiondone = (double)f / r.size();
			std::cout << 100 * fractiondone << "%" << std::endl;
		}
	}
	auto CalculationEndTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> CalculationTime = CalculationEndTime - CalculationStart;
	std::cout << "Force Calculation Time: " << CalculationTime.count() << "s" << std::endl;
	std::cout << "Number of Unfound Points: " << ErrorCount << std::endl;
	
	return std::make_tuple(Fz, Fr, Ftheta);
}



#define Circ_GF

#endif /* Circ_GF */