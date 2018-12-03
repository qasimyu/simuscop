// ***************************************************************************
// psiFunc.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _PSIFUNC_H
#define _PSIFUNC_H

/*** gamma pdf ***/
double gammapdf(double x, double k, double theta);

/*** digamma function ***/
double digamma(double x);

/*** trigamma function ***/
double trigamma(double x);

/*** polygamma function ***/
double psi(int degree, double x);


#endif

