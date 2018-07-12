// ***************************************************************************
// MyDefine.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <vector>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

string current_version = "1.0";

Config config;
Genome genome;
Profile profile;
ThreadPool* threadPool;
SeqWriter* swp;

//****** normalize matrix ******//
Matrix<double> norm_trans(Matrix<double> &T, double thres) {
	Matrix<double> ret(T);
	double **p1 = ret.getEntrance();
	double ZERO_FINAL = 2.2204e-16;
	for(size_t i = 0; i < ret.getROWS(); i++) {
		Matrix<double> temp = ret.Row(i);
		double tmp = temp.sum();
		if(p1[i][i] < thres*tmp) {
			temp.getEntrance()[0][i] = 0;
			Matrix<double> a = temp*((1-thres)/temp.sum());
			a.getEntrance()[0][i] = thres;
			ret.setRow(i, a);
		}
		else {
			Matrix<double> a = temp/(tmp+ZERO_FINAL);
			ret.setRow(i, a);
		}
	}
	return ret;
}

//****** pdf of normal distribution ******//
double normpdf(double x, double mu, double sigma) {
	double PI = 3.1415926;
	return exp(-pow(x-mu,2)/(2*pow(sigma,2)))/(sqrt(2*PI)*sigma);
}

//****** pdf of student't distribution ******//
double stupdf(double x, double mu, double nu, double sigma) {
	double PI = 3.1415926;
	double ret = lgamma((nu+1)*0.5)-(nu+1)*0.5*log(1+pow(x-mu,2)/(nu*sigma*sigma))
				-lgamma(nu*0.5)-0.5*log(nu)-log(sigma)-0.5*log(PI);
	ret = exp(ret);
	return ret;
}

//****** calculate median value of a matrix ******//
double median(Matrix<double> &mat) {
	int n = mat.getROWS()*mat.getCOLS();
	Matrix<double> T = mat.reshape(1,n);
	double *p = T.getEntrance()[0];
	return median(p,n);
}

//****** calculate median value ******//
double median(double *p, int n) {
	QuickSort(p,n);
	if(n%2 != 0) {
		return p[n/2];
	}
	else {
		return (p[n/2]+p[n/2-1])/2;
	}
}

//****** quick sort ******//
int Partition(double *a, int low, int high) {
	int i = low, j = high;
	double pivot = a[low];
	while(i < j) {
		while(i<j && a[j]>=pivot) j--;
		a[i]=a[j];
		while(i<j && a[i]<=pivot) i++;
		a[j]=a[i];
	}
	a[i]=pivot;
    return i;
}

void QSort(double *a, int low, int high) {
	int pivotloc;
	if(low < high) {
		pivotloc = Partition(a, low, high);
		QSort(a, low, pivotloc-1);
		QSort(a, pivotloc+1, high);
	}
}

void QuickSort(double *a, int n){
	QSort(a,0,n-1);
}

//****** calculate standard deviation ******//
double stdev(double *p, int n, int flag) {
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
	if(flag == 0)
		sigma = sqrt(sigma/(n-1));
	else
		sigma = sqrt(sigma/n);
	return sigma;
}

//****** minimum element of an array ******//
double min(double *p, int n,int &indx) {
	double minValue = p[0];
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
double max(double *p, int n, int &indx) {
	double maxValue = p[0];
	indx = 0;
	for(int i = 1; i < n; i++) {
		if(p[i] > maxValue) {
			maxValue = p[i];
			indx = i;
		}
	}
	return maxValue;
}

//****** produce a matrix containing random intergers according to the probability distribution ******//
Matrix<int> randsrc(int m, int n, Matrix<int> alphabet, Matrix<double> aprob, bool iscdf) {
	Matrix<int> ret;
	if(m <= 0 || n <= 0) {
		return ret;
	}
	
	assert(alphabet.getROWS() == 1 && alphabet.getCOLS() > 0);
	assert(aprob.getROWS() == 1 && aprob.getCOLS() == alphabet.getCOLS());
	
	Matrix<double> prob;
	if(iscdf) {
		prob = aprob;
	}
	else {
		prob = aprob.cumsum();
	}

	srand(time(0));
	
	ret.resize(m, n, false);
	int** p = ret.getEntrance();
	int i, j, k;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			p[i][j] = randsrc(alphabet, prob, true);
		}
	}
	
	return ret;
}

//****** produce a random integer according to the probability distribution ******//
int randsrc(Matrix<int>& alphabet, Matrix<double>& aprob, bool iscdf) {
	assert(alphabet.getROWS() == 1 && alphabet.getCOLS() > 0);
	assert(aprob.getROWS() == 1 && aprob.getCOLS() == alphabet.getCOLS());
	
	int ac = aprob.getCOLS();
	double** p1;
	if(iscdf) {
		p1 = aprob.getEntrance();
	}
	else {
		Matrix<double> prob = aprob.cumsum();
		p1 = prob.getEntrance();
	}
	
	int** p2 = alphabet.getEntrance();
	double r = randomDouble(0, 1);
	for(size_t k = 0; k < ac; k++) {
		if(r <= p1[0][k]) {
			return p2[0][k];
		}
	}
	return p2[0][ac-1];
}

//****** produce a matrix containing random index according to the probability distribution ******//
Matrix<int> randIndx(int m, int n, Matrix<double> aprob, bool iscdf) {
	Matrix<int> ret;
	if(m <= 0 || n <= 0) {
		return ret;
	}
	
	assert(aprob.getROWS() == 1 && aprob.getCOLS() >= 1);
	
	Matrix<double> prob;
	if(iscdf) {
		prob = aprob;
	}
	else {
		prob = aprob.cumsum();
	}

	srand(time(0));
	
	ret.resize(m, n, false);
	int** p = ret.getEntrance();
	int i, j, k;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			p[i][j] = randIndx(prob, true);
		}
	}
	
	return ret;
}

int randIndx(Matrix<double>& aprob, bool iscdf) {
	int ac = aprob.getCOLS();
	double** p1;
	Matrix<double> prob;
	if(iscdf) {
		p1 = aprob.getEntrance();
	}
	else {
		prob = aprob.cumsum();
		p1 = prob.getEntrance();
	}
	
	double r = randomDouble(0, 1);
	for(size_t k = 0; k < ac; k++) {
		if(r <= p1[0][k]) {
			return k;
		}
	}
	return ac-1;
}

//****** produce random double ******//
double randomDouble(double start, double end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}

//****** produce random int ******//
long randomInteger(long start, long end) {
	return start+(end-start)*(rand()/(RAND_MAX+1.0));
}

//****** trim string ******//
string trim(const string &str, const char *charlist) {
	string ret(str);
	size_t indx = ret.find_first_not_of(charlist);
	if(indx != string::npos) {
		ret.erase(0, indx);
		indx = ret.find_last_not_of(charlist);
		ret.erase(indx+1);
	}
	else {
		ret.erase();
	}
	return ret;
}

//****** abbreviation of chromosome name ******//
string abbrOfChr(string chr) {
	string abbr_chr = chr;
	size_t i = abbr_chr.find("chrom");
	if(i == string::npos) {
		i = abbr_chr.find("chr");
		if(i != string::npos) {
			abbr_chr = abbr_chr.substr(i+3,abbr_chr.size()-3);
		}
	}
	else {
		abbr_chr = abbr_chr.substr(i+5,abbr_chr.size()-5);
	}
	return abbr_chr;
}


