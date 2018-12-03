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

extern double ZERO_FINAL;

extern Config config;
extern Genome genome;
extern Profile profile;
extern ThreadPool* threadPool;
extern SeqWriter* swp;

//*** declaration of functions ***//
Matrix<double> norm_trans(Matrix<double> &T, double thres);

double normpdf(double x, double mu, double sigma);
double stupdf(double x, double mu, double nu, double sigma);

//****** quick sort ******//
template <class numtype>
int Partition(numtype *a, int low, int high) {
	int i = low, j = high;
	numtype pivot = a[low];
	while(i < j) {
		while(i<j && a[j]>=pivot) j--;
		a[i] = a[j];
		while(i<j && a[i]<=pivot) i++;
		a[j] = a[i];
	}
	a[i] = pivot;
    return i;
}

template <class numtype>
void QSort(numtype *a, int low, int high) {
	int pivotloc;
	if(low < high) {
		pivotloc = Partition(a, low, high);
		QSort(a, low, pivotloc-1);
		QSort(a, pivotloc+1, high);
	}
}

template <class numtype>
void QuickSort(numtype *a, int n){
	QSort(a, 0, n-1);
}

//****** calculate median value ******//
template <class numtype>
numtype median(numtype *p, int n) {
	QuickSort(p, n);
	if(n%2 != 0) {
		return p[n/2];
	}
	else {
		return (p[n/2]+p[n/2-1])/2;
	}
}

//****** calculate median value of a matrix ******//
template <class numtype>
numtype median(Matrix<numtype> &mat) {
	int n = mat.getROWS()*mat.getCOLS();
	Matrix<numtype> T = mat.reshape(1, n);
	numtype *p = T.getEntrance()[0];
	return median(p, n);
}

//****** calculate median value of a vector ******//
template <class numtype>
numtype median(vector<numtype> &x) {
	int n = x.size();
	numtype *p = new numtype[n];
	for(int i = 0; i < n; i++) {
		p[i] = x[i];
	}
	numtype ret = median(p, n);
	delete[] p;
	return ret;
}

//****** calculate standard deviation ******//
template <class numtype>
double stdev(numtype *p, int n, int flag) {
	double average = 0;
	int i;
	for(i = 0; i < n; i++) {
		average += p[i];
	}
	average /= n;
	double sigma = 0;
	for(i = 0; i < n; i++) {
		sigma += pow(p[i]-average,2);
	}
	if(flag == 0) {
		sigma = sqrt(sigma/(n-1));
	}
	else {
		sigma = sqrt(sigma/n);
	}
	return sigma;
}

//****** calculate variance ******//
template <class numtype>
double var(vector<numtype> &x) {
	int n = x.size();
	numtype *p = new numtype[n];
	for(int i = 0; i < n; i++) {
		p[i] = x[i];
	}
	double s = stdev(p, n, 1);
	delete[] p;
	return s*s;
}

//****** minimum element of an array ******//
template <class numtype>
numtype min(numtype *p, int n,int &indx) {
	numtype minValue = p[0];
	indx = 0;
	for(int i = 1; i < n; i++) {
		if(p[i] < minValue) {
			minValue = p[i];
			indx = i;
		}
	}
	return minValue;
}

//****** maximum element of an array ******//
template <class numtype>
numtype max(numtype *p, int n, int &indx) {
	numtype maxValue = p[0];
	indx = 0;
	for(int i = 1; i < n; i++) {
		if(p[i] > maxValue) {
			maxValue = p[i];
			indx = i;
		}
	}
	return maxValue;
}

//****** percentile ******//
template <class numtype>
numtype prctile(vector<numtype> &x, double perc) {
	int n = x.size();
	if(perc < 0) {
		perc = 0;
	}
	if(perc > 1) {
		perc = 1;
	}
	numtype *p = new numtype[n];
	for(int i = 0; i < n; i++) {
		p[i] = x[i];
	}
	QuickSort(p, n);
	int j = n*perc;
	j = max(0, min(j, n-1));
	numtype ret = p[j];
	delete[] p;
	return ret;
}

int randsrc(Matrix<int>& alphabet, Matrix<double>& aprob, bool iscdf);
Matrix<int> randsrc(int m, int n, Matrix<int> alphabet, Matrix<double> aprob, bool iscdf);
int randIndx(Matrix<double>& aprob, bool iscdf);
int randIndx(double *cdf, int ac);
Matrix<int> randIndx(int m, int n, Matrix<double> aprob, bool iscdf);
double randomDouble(double start, double end);
long randomInteger(long start, long end);

string trim(const string &str, const char *charlist = " \t\r\n");
string abbrOfChr(string chr);

int getIndexOfBase(char base);
bool getNextLine(ifstream& ifs, string& line, int& lineNum);
char* getComplementSeq(char* sequence);
int calculateGCPercent(char* sequence);
double calculateGCContent(char* sequence);

//void mergeChrFastqFiles(string popu, string chr);
//void mergePopuFastqFiles(string popu, vector<string>& chromosomes);
//void mergeFastqFiles(string fn);

//*** end of declaration ***//


#endif

