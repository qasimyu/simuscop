// ***************************************************************************
// MyDefine.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MYDEFINE_H
#define _MYDEFINE_H

#include <string>

#include "Config.h"
#include "Genome.h"
#include "Profile.h"
#include "ThreadPool.h"
#include "SeqWriter.h"
#include "Matrix.h"

using namespace std;

class MyDefine {		
	public:
		MyDefine();
		~MyDefine();
};

extern string current_version;

extern Config config;
extern Genome genome;
extern Profile profile;
extern ThreadPool* threadPool;
extern SeqWriter* swp;

//*** declaration of functions ***//
Matrix<double> norm_trans(Matrix<double> &T, double thres);

double normpdf(double x, double mu, double sigma);
double stupdf(double x, double mu, double nu, double sigma);

void QuickSort(double *a, int n);

double median(double *p, int n);
double median(Matrix<double> &mat);

double stdev(double *p, int n, int flag);

double min(double *p, int n, int &indx);
double max(double *p, int n, int &indx);

int randsrc(Matrix<int>& alphabet, Matrix<double>& aprob, bool iscdf);
Matrix<int> randsrc(int m, int n, Matrix<int> alphabet, Matrix<double> aprob, bool iscdf);
int randIndx(Matrix<double>& aprob, bool iscdf);
Matrix<int> randIndx(int m, int n, Matrix<double> aprob, bool iscdf);
double randomDouble(double start, double end);
long randomInteger(long start, long end);

string trim(const string &str, const char *charlist = " \t\r\n");
string abbrOfChr(string chr);

//void mergeChrFastqFiles(string popu, string chr);
//void mergePopuFastqFiles(string popu, vector<string>& chromosomes);
//void mergeFastqFiles(string fn);

//*** end of declaration ***//


#endif

