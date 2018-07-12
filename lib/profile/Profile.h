// ***************************************************************************
// Profile.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _PROFILE_H
#define _PROFILE_H

#include <vector>
#include <map>
#include <string>
#include <pthread.h>

#include "Matrix.h"

using namespace std;

class Profile {
	private:
		int meanReadLength;
		int estimatedCoverage;
		
		int binCount;
		int randomCount;
		Matrix<int> baseAlphabet;
		Matrix<int> qualityAlphabet;
		Matrix<int> insertSizeAlphabet;
		Matrix<int>* gcAlphabets;

		int minBaseQuality;
		int maxBaseQuality;
		
		double stdInsertSize;
		
		string nucleotides;
		
		Matrix<double>** conditionalSubstitutionProbs1;
		Matrix<double>** conditionalSubstitutionProbs2;
		Matrix<double> initSubstitutionProbs1;
		Matrix<double> initSubstitutionProbs2;
		Matrix<double> baseQualityDist;
		Matrix<double> errorBaseQualityDist;
		Matrix<double> insertSizeDist;
		Matrix<double>* gcDist;
		
		Matrix<int>** conditionalSubstitutionIndxs1;
		Matrix<int>** conditionalSubstitutionIndxs2;
		Matrix<int> initSubstitutionIndxs1;
		Matrix<int> initSubstitutionIndxs2;
		Matrix<int> baseQualityIndxs;
		Matrix<int> errorBaseQualityIndxs;
		Matrix<int> insertSizeIndxs;
		Matrix<int>* gcIndxs;
		
		//GC-content bias
		map<int, vector<double> > readCountsByGC;
		double gcFactors[101];
		double gcStds[101];
		Matrix<double>* gcValues;
		
		void init();
		int getIndexOfNucleotide(char nucleotide);
		
		int countGC(string chr, long position);
		void evaluateGC();
		
		void saveResults(string bamFile, string outFile);
		void initRandNumberPool();
		void initSamplingPool();
		
		void load(string modelFile);
		void initGCFactors();
		void initParas(bool isLoaded);
		void estimateCoverage();
		
		bool getNextLine(ifstream& ifs, string& line, int& lineNum);
		
	public:
		Profile();
		Profile(string &nucleotides);
		~Profile();
		
		int calculateGCPercent(char* sequence);
		
		int processRead(char* read);
		
		int yieldInsertSize();
		double getStdInsertSize();
		int getMaxInsertSize();
		double getGCFactor(int gc);
		int getReplicates(int gc);
		int getInitBaseIndx1(int base);
		int getInitBaseIndx2(int base);
		int getCondBaseIndx1(int preBase, int base, int binIndx);
		int getCondBaseIndx2(int preBase, int base, int binIndx);
		int getBaseQuality(int base);
		int getErrorBaseQuality(int base);
		
		void train(string modelFile);
		void train(string bamFile, string samtools, string outFile);
		void predict(char* refSeq, char* results, int num, int isRead1);
		
};


#endif

