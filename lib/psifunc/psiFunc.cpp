// ***************************************************************************
// psiFunc.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "psiFunc.h"

using namespace std;

double gammapdf(double x, double k, double theta) {
	return exp((k-1)*log(x)-x/theta-lgamma(k)-k*log(theta));
}
		
double digamma(double x) {
	double euler_mascheroni = 0.57721566490153286060;
	double r;
	double ret;
	double x2;

	if(x <= 0) {
		cerr << "The input argument must be positive!" << endl;
		exit(-1);
	}

	x2 = x;
	ret = 0.0;
	
	if(x2 <= 0.00001) {
		ret = -euler_mascheroni-1.0/x2;
		return ret;
	}
	
	while(x2 < 8.5) {
		ret -= 1.0/x2;
		x2 += 1.0;
	}
	
	r = 1.0/x2;
	ret += log(x2)-0.5*r;
	r *= r;
	ret -= r*(1.0/12-r/120+r*r/252);

	return ret;
}

double trigamma(double x) {
	double a = 0.0001;
	double b = 5.0;
	double b2 =  0.1666666667;
	double b4 = -0.03333333333;
	double b6 =  0.02380952381;
	double b8 = -0.03333333333;
	double ret;
	double y;
	double z;
	
	z = x;
	
	if(x <= a) {
		ret = 1.0/(x*x);
		return ret;
	}
	
	ret = 0.0;

	while(z < b) {
		ret += 1.0/(z*z);
		z += 1.0;
	}
	
	y = 1.0/(z*z);

	ret += 0.5*y+(1+y*(b2+y*(b4+y*(b6+y*b8))))/z;

	return ret;
}

double psi(int degree, double x) {
	if(!(degree == 0 || degree == 1)) {
		cerr << "Error: The input degree must be one of the values in [0,1]!" << endl;
		exit(-1);
	}
	if(x <= 0) {
		cerr << "Error: The input argument must be positive!" << endl;
		exit(-1);
	}
	if(degree == 0)
		return digamma(x);
	else
		return trigamma(x);
}
