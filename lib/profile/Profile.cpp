// ***************************************************************************
// Profile.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
//#include <cctype>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <random>
#include <unistd.h>

#include "split.h"
#include "MyDefine.h"
#include "Profile.h"

Profile::Profile(void) {
	nucleotides = "ACTGN";
	init();
}

Profile::Profile(string &nucleotides) {
	this->nucleotides = nucleotides;
	init();
}

Profile::~Profile() {
	int i, j;
		
	baseAlphabet.clear();
	qualityAlphabet.clear();
	insertSizeAlphabet.clear();
	if(gcAlphabets != NULL) {
		for(i = 0; i < 101; i++) {
			gcAlphabets[i].clear();
		}
		delete[] gcAlphabets;
	}
	
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < binCount; j++) {
			if(conditionalSubstitutionProbs1 != NULL) {
				conditionalSubstitutionProbs1[i][j].clear();
			}
			if(conditionalSubstitutionIndxs1 != NULL) {
				conditionalSubstitutionIndxs1[i][j].clear();
			}
			if(conditionalSubstitutionProbs2 != NULL) {
				conditionalSubstitutionProbs2[i][j].clear();
			}
			if(conditionalSubstitutionIndxs2 != NULL) {
				conditionalSubstitutionIndxs2[i][j].clear();
			}
		}
		if(conditionalSubstitutionProbs1 != NULL) {
			delete[] conditionalSubstitutionProbs1[i];
		}
		if(conditionalSubstitutionIndxs1 != NULL) {
			delete[] conditionalSubstitutionIndxs1[i];
		}
		if(conditionalSubstitutionProbs2 != NULL) {
			delete[] conditionalSubstitutionProbs2[i];
		}
		if(conditionalSubstitutionIndxs2 != NULL) {
			delete[] conditionalSubstitutionIndxs2[i];
		}
	}
	if(conditionalSubstitutionProbs1 != NULL) {
		delete[] conditionalSubstitutionProbs1;
	}
	if(conditionalSubstitutionIndxs1 != NULL) {
		delete[] conditionalSubstitutionIndxs1;
	}
	if(conditionalSubstitutionProbs2 != NULL) {
		delete[] conditionalSubstitutionProbs2;
	}
	if(conditionalSubstitutionIndxs2 != NULL) {
		delete[] conditionalSubstitutionIndxs2;
	}

	initSubstitutionIndxs1.clear();
	initSubstitutionIndxs2.clear();
	baseQualityIndxs.clear();
	errorBaseQualityIndxs.clear();
	insertSizeIndxs.clear();

	if(gcDist != NULL) {
		for(i = 0; i < 101; i++) {
			gcDist[i].clear();
		}
		delete[] gcDist;
	}
	if(gcIndxs != NULL) {
		for(i = 0; i < 101; i++) {
			gcIndxs[i].clear();
		}
		delete[] gcIndxs;
	}
	if(gcValues != NULL) {
		for(i = 0; i < 101; i++) {
			gcValues[i].clear();
		}
		delete[] gcValues;
	}
}

void Profile::init() {
	minBaseQuality = 33;
	maxBaseQuality = 126;
	randomCount = 5000;
	binCount = 50;

	meanReadLength = 0;

	int N = nucleotides.length();
	
	conditionalSubstitutionProbs1 = new Matrix<double>*[N];
	conditionalSubstitutionProbs2 = new Matrix<double>*[N];
	for(size_t i = 0; i < N; i++) {
		conditionalSubstitutionProbs1[i] = new Matrix<double>[binCount];
		conditionalSubstitutionProbs2[i] = new Matrix<double>[binCount];
		for(size_t j = 0; j < binCount; j++) {
			Matrix<double> csProbs(N, N, true);
			conditionalSubstitutionProbs1[i][j] = csProbs;
			conditionalSubstitutionProbs2[i][j] = csProbs;
		}
	}
	initSubstitutionProbs1.resize(N, N, true);
	initSubstitutionProbs2.resize(N, N, true);
	
	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	baseQualityDist.resize(N, baseQualtiyCount, true);
	errorBaseQualityDist.resize(N, baseQualtiyCount, true);
	
	insertSizeDist.resize(1, 1000, true);
	stdInsertSize = 0;

	conditionalSubstitutionIndxs1 = NULL;
	conditionalSubstitutionIndxs2 = NULL;
	
	gcAlphabets = NULL;
	gcDist = NULL;
	gcIndxs = NULL;
	gcValues = NULL;
}

int Profile::getIndexOfNucleotide(char nucleotide) {
	for(size_t i = 0; i < nucleotides.length(); i++) {
		if(nucleotide == nucleotides[i]) {
			return i;
		}
	}
	return nucleotides.length()-1;
}

int Profile::calculateGCPercent(char* sequence) {
	int gcCount = 0;
	int nCount = 0;
	if(sequence == NULL || sequence[0] == '\0') {
		return 0;
	}
	int n = strlen(sequence);
	for(size_t i = 0; i < n; i++) {
		if(sequence[i] == 'G' || sequence[i] == 'C') {
			gcCount++;
		}
		else if(sequence[i] == 'N') {
			nCount++;
		}
	}
	if(nCount > n/2) {
		return 0;
	}
	return 100*gcCount/(n-nCount);
}

int Profile::processRead(char* read) {
	if(read == NULL || read[0] == '\0') {
		return 0;
	}

	static int wgs = (genome.getInTargets().empty())? 1:0;
	
	static unsigned long readCount = 0;
	static unsigned long maxCount = 200000000;
	static vector<string>& chromosomes = genome.getChroms();
	vector<string>::iterator v_it;
	
	int i, j, n;
	
	/***check the number of fields of read***/
	char* elems[20];
	i = split(read, '\t', elems);
	if(i < 11) {
		cerr << "Error: malformed read " <<
				", there should be 11 mandatory fields" << endl;
		cerr << read << endl;
		exit(1);
	}
	
	/***mandatory fields***/
	int flag = atoi(elems[1]);
	string chr(elems[2]);
	long position = atol(elems[3]);
	int mapQuality = atoi(elems[4]);
	char* cigar = elems[5];
	int tlen = atoi(elems[8]);
	char* readSeq = elems[9];
	char* baseQuality = elems[10];
	
	if(position == 0) {
		return 0;
	}

	if(mapQuality < 15) {
		return 0;
	}

	v_it = find(chromosomes.begin(), chromosomes.end(), chr);
	if(v_it == chromosomes.end()) {
		return 0;
	}
	
	/***update read counts distribution associated with GC-content***/
	
	i = countGC(chr, position);
	if(i == 0) {
		return 0;	
	}
	
	
	n = strlen(cigar);
	for(i = 0; i < n-1; i++) {
		if(!(cigar[i] >= '0' && cigar[i] <= '9')) {
			break;
		}
	}
	if(i < n-1 || cigar[i] != 'M') {
		return 0;
	}
	
	if(strcmp(readSeq, "*") == 0) {
		return 0;
	}
	
	char* refSeq = genome.getSubSequence(chr, position-1, strlen(readSeq));
	char* altSeq = genome.getSubAltSequence(chr, position-1, strlen(readSeq));
	
	int isRead1 = 1;
	if(tlen < 0) {
		reverse(refSeq, refSeq+strlen(refSeq));
		refSeq = Segment::getComplementSeq(refSeq);
		reverse(altSeq, altSeq+strlen(altSeq));
		altSeq = Segment::getComplementSeq(altSeq);
		reverse(readSeq, readSeq+strlen(readSeq));
		readSeq = Segment::getComplementSeq(readSeq);
		reverse(baseQuality, baseQuality+strlen(baseQuality));
		isRead1 = 0;
	}	
	
	/***update initSubstitutionProbs and conditionalSubstitutionProbs***/
	int indx, indx1, indx2, binIndx;
	indx1 = getIndexOfNucleotide(refSeq[0]);
	indx2 = getIndexOfNucleotide(readSeq[0]);
	if(indx1 != -1 && indx2 != -1) {
		if(altSeq[0] == readSeq[0]) {
			indx1 = indx2;
		}
		if(isRead1) {
			i = initSubstitutionProbs1.get(indx1, indx2);
			initSubstitutionProbs1.set(indx1, indx2, i+1);
		}
		else {
			i = initSubstitutionProbs2.get(indx1, indx2);
			initSubstitutionProbs2.set(indx1, indx2, i+1);
		}
	}
	n = strlen(refSeq);
	for(i = 1; i < n; i++) {
		indx = getIndexOfNucleotide(refSeq[i]);
		
		if(readSeq[i-1] == altSeq[i-1]) {
			indx1 = getIndexOfNucleotide(altSeq[i-1]);
		}
		else {
			indx1 = getIndexOfNucleotide(refSeq[i-1]);
		}
		
		//indx1 = getIndexOfNucleotide(readSeq[i-1]);
		indx2 = getIndexOfNucleotide(readSeq[i]);
		binIndx = i*binCount/n;
		if(indx != -1 && indx1 != -1 && indx2 != -1) {
			if(altSeq[i] == readSeq[i]) {
				indx = indx2;
			}
			if(isRead1) {
				j = conditionalSubstitutionProbs1[indx][binIndx].get(indx1, indx2);
				conditionalSubstitutionProbs1[indx][binIndx].set(indx1, indx2, j+1);
			}
			else {
				j = conditionalSubstitutionProbs2[indx][binIndx].get(indx1, indx2);
				conditionalSubstitutionProbs2[indx][binIndx].set(indx1, indx2, j+1);
			}
		}
	}
	
	/***update insert sizeset distribution***/
	//if(config.isPairedEnd() && tlen > 0) {
	if(tlen > 0) {
		if(tlen > insertSizeDist.getCOLS()-1) {
			insertSizeDist.resize(1, tlen+1, true);
		}
		insertSizeDist.set(0, tlen, insertSizeDist.get(0, tlen)+1);
	}
	
	/***update base quality distribution***/
	if(strlen(baseQuality) == strlen(readSeq)) {
		n = strlen(readSeq);
		for(i = 0; i < n; i++) {
			indx1 = getIndexOfNucleotide(readSeq[i]);
			if(indx1 == -1) {
				continue;
			}
			int k = baseQuality[i];
			//int k = baseQuality[i]-33;
			
			if(k >= minBaseQuality && k <= maxBaseQuality) {
				k -= minBaseQuality;
				if(refSeq[i] == readSeq[i] || altSeq[i] == readSeq[i]) {
					baseQualityDist.set(indx1, k, baseQualityDist.get(indx1, k)+1);
				}
				else {
					errorBaseQualityDist.set(indx1, k, errorBaseQualityDist.get(indx1, k)+1);
				}
			}
		}
	}

	/***update mean read length***/
	readCount++;
	meanReadLength = ((readCount-1)*meanReadLength+strlen(readSeq))/readCount;	
	
	/*
	if(config.isVerbose() && readCount%1000000 == 0) {
		cerr << readCount << " reads processed!" << endl;
	}
	*/
	if(readCount%1000000 == 0) {
		cerr << readCount << " reads processed!" << endl;
	}
	
	delete[] refSeq;
	delete[] altSeq;
	
	if(wgs == 1 && readCount >= maxCount) {
		return 2;
	}
	/*
	else {
		return 0;
	}
	*/
	if(wgs == 0 && readCount >= 2*maxCount) {
		return 2;
	}
	return 1;
	
}

int Profile::countGC(string chr, long position) {
	/***evaluate GC-content effect on read counts***/	
	static string preChr = "";
	static long refLen;
	static unsigned int winSize = Segment::getFragmentSize();
	static long rightPos = 0;
	static long leftPos = -1;
	static int GC = -1;
	static double rc = 0;
	static int wgs = (genome.getInTargets().empty())? 1:0;

	static int targetIndx = -1;
	
	char* refSeq;

	position -= 1;

	if(chr.compare("X") == 0 || chr.compare("Y") == 0 || chr.compare("M") == 0 || chr.compare("MT") == 0) {
		return -1;
	}
	
	if(preChr.compare(chr) == 0) {
		if(refLen == 0) {
			return 0;
		}
		if(position < leftPos) {
			return 0;
		}
		if(position >= leftPos && position <= rightPos) {
			rc++;
			return 1;
		}
		if(GC != -1 && rc > 0) {
			if(wgs == 1) {
				readCountsByGC[GC].push_back(rc);
			}
			else {
				int targetSize = rightPos-leftPos+1;
				rc = winSize*rc/targetSize;
				readCountsByGC[GC].push_back(rc);
			}
		}
		if(wgs == 1) {
			rightPos += winSize;
			while(rightPos < position) {
				rightPos += winSize;
			}
			rightPos = min(rightPos, refLen-1);
			leftPos = rightPos-winSize+1;
			refSeq = genome.getSubSequence(chr, leftPos, winSize);
			rc = 1;
		}
		else {
			vector<Target>& targets = genome.getInTargets(chr);
			targetIndx++;
			for(; targetIndx < targets.size(); targetIndx++) {
				if(targets[targetIndx].epos-1 >= position) {
					break;
				}
			}
			if(targetIndx < targets.size()) {
				if(targets[targetIndx].epos <= refLen) {
					rightPos = targets[targetIndx].epos-1;
					leftPos = targets[targetIndx].spos;
					refSeq = genome.getSubSequence(chr, leftPos, rightPos-leftPos+1);
					if(leftPos <= position) {
						rc = 1;		
					}
					else {
						rc = 0;
					}
				}
				else {
					refSeq = NULL;
					rc = 0;
				}
			}
			else {
				refSeq = NULL;
				rc = 0;
				leftPos = refLen;
				rightPos = refLen;
			}
		}
		if(refSeq != NULL) {
			GC = calculateGCPercent(refSeq);
			delete[] refSeq;
		}
		else {
			GC = -1;
		}
		return rc;
	}
	else {
		if(preChr.compare("") != 0 && GC != -1) {
			if(wgs == 1) {
				readCountsByGC[GC].push_back(rc);
			}
			else {
				int targetSize = rightPos-leftPos+1;
				rc = winSize*rc/targetSize;
				readCountsByGC[GC].push_back(rc);
			}
		}
		
		preChr = chr;
		refLen = genome.getChromLen(chr);
		if(refLen == 0) {
			refSeq = NULL;
			rc = 0;
			leftPos = -1;
			rightPos = -1;
		}
		else {
			if(winSize > refLen) {
				cerr << "[GC evaluation] Window size greater than chromosome length of " << chr << ", adjusting to chromosome length: " << refLen << endl;
				winSize = refLen;
			}

			rightPos = -1;
			if(wgs == 1) {
				rightPos += winSize;
				while(rightPos < position) {
					rightPos += winSize;
				}
				rightPos = min(rightPos, refLen-1);
				leftPos = rightPos-winSize+1;
				refSeq = genome.getSubSequence(chr, leftPos, winSize);
				rc = 1;
			}
			else {
				vector<Target>& targets = genome.getInTargets(chr);
				for(targetIndx = 0; targetIndx < targets.size(); targetIndx++) {
					if(targets[targetIndx].epos-1 >= position) {
						break;
					}
				}
				if(targetIndx < targets.size()) {
					if(targets[targetIndx].epos <= refLen) {
						rightPos = targets[targetIndx].epos-1;
						leftPos = targets[targetIndx].spos;
						refSeq = genome.getSubSequence(chr, leftPos, rightPos-leftPos+1);
						if(leftPos <= position) {
							rc = 1;		
						}
						else {
							rc = 0;
						}
					}
					else {
						refSeq = NULL;
						rc = 0;
					}
				}
				else {
					refSeq = NULL;
					rc = 0;
					leftPos = refLen;
					rightPos = refLen;
				}
			}
		}
		if(refSeq != NULL) {
			GC = calculateGCPercent(refSeq);
			delete[] refSeq;
		}
		else {
			GC = -1;
		}
		return rc;
	}
	
}

void Profile::estimateCoverage() {
	int i, j;
	vector<double> tmp;
	map<int, vector<double> >::iterator it;
	for(it = readCountsByGC.begin(); it != readCountsByGC.end(); it++) {
		vector<double>& readCounts = (*it).second;
		tmp.insert(tmp.begin(), readCounts.begin(), readCounts.end());
	}
	sort(tmp.begin(), tmp.end());

	int medRC;
	long n = tmp.size();
	if(n%2 != 0) {
		medRC = tmp[n/2];
	}
	else {
		medRC = (tmp[n/2]+tmp[n/2-1])/2;
	}
	tmp.clear();

	//cerr << "meanReadLength:" << meanReadLength << endl;

	unsigned int winSize = Segment::getFragmentSize();
	estimatedCoverage = meanReadLength*medRC/winSize;
	//cerr << "estimatedCoverage:" << estimatedCoverage << endl;
	
}

void Profile::initGCFactors() {
	int i;
	for(i = 0; i < 101; i++) {
		gcFactors[i] = 5;
		gcStds[i] = 0;
	}
	estimatedCoverage = 0;
}

void Profile::initParas(bool isLoaded) {
	int i, j, k;
	/***normalize probability matrix***/
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < binCount; j++) {
			conditionalSubstitutionProbs1[i][j].normalize(0);
			conditionalSubstitutionProbs2[i][j].normalize(0);
		}
	}
	
	initSubstitutionProbs1.normalize(0);
	initSubstitutionProbs2.normalize(0);
	
	baseQualityDist.normalize(0);

	errorBaseQualityDist.normalize(0);
	
	if(!isLoaded) {
		int maxCount = 0;
		j = 0;
		for(i = 0; i < insertSizeDist.getCOLS(); i++) {
			if(insertSizeDist.get(0, i) > maxCount) {
				maxCount = insertSizeDist.get(0, i);
				j = i;
			}
		}
		
		for(i = j*5; i < insertSizeDist.getCOLS(); i++) {
			insertSizeDist.set(0, i, 0);
		}

		insertSizeDist.normalize(0);
		double meanTlen = 0;
		for(i = 0; i < insertSizeDist.getCOLS(); i++) {
			meanTlen += insertSizeDist.get(0, i)*i;
		}
		stdInsertSize = 0;
		for(i = 0; i < insertSizeDist.getCOLS(); i++) {
			stdInsertSize += insertSizeDist.get(0, i)*pow(i-meanTlen, 2);
		}
		stdInsertSize = sqrt(stdInsertSize);
		//cerr << "stdInsertSize: " << stdInsertSize << endl;
		//insertSizeDist.clear();
	}
	else {
		int N = nucleotides.length();
		baseAlphabet.resize(1, N, false);
		for(i = 0; i < N; i++) {
			baseAlphabet.set(0, i, i);
		}
		int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
		qualityAlphabet.resize(1, baseQualtiyCount, false);
		for(i = 0, j = minBaseQuality; i < baseQualtiyCount; i++, j++) {
			qualityAlphabet.set(0, i, j);
		}
	
		if(config.isPairedEnd() && stdInsertSize > 0) {
			int meanInsertSize = config.getIntPara("insertSize")+1;
			int intervalLen = 6 * stdInsertSize;
			int minInsertSize = max(meanInsertSize - intervalLen/2, config.getIntPara("readLength"));
			int maxInsertSize = 2*meanInsertSize - minInsertSize;
			//cerr << stdInsertSize << endl;
			//cerr << minInsertSize << '\t' << maxInsertSize << '\t' << meanInsertSize << endl;
			insertSizeAlphabet.resize(1, maxInsertSize-minInsertSize+1, false);
			for(i = 0; i < insertSizeAlphabet.getCOLS(); i++) {
				insertSizeAlphabet.set(0, i, minInsertSize++);
			}
		
			insertSizeDist.resize(1, insertSizeAlphabet.getCOLS(), false);
			insertSizeIndxs.resize(1, randomCount*insertSizeAlphabet.getCOLS(), false);
			for(i = 0; i < insertSizeDist.getCOLS(); i++) {
				double pdf = normpdf(insertSizeAlphabet.get(0, i), meanInsertSize, stdInsertSize);
				insertSizeDist.set(0, i, pdf);
			}
			insertSizeDist.normalize(0);
			//insertSizeCdf = insertSizeDist.cumsum();
			//insertSizeDist.clear();
		}
		
		conditionalSubstitutionIndxs1 = new Matrix<int>*[N];
		conditionalSubstitutionIndxs2 = new Matrix<int>*[N];
		for(i = 0; i < N; i++) {
			conditionalSubstitutionIndxs1[i] = new Matrix<int>[binCount];
			conditionalSubstitutionIndxs2[i] = new Matrix<int>[binCount];
			for(size_t j = 0; j < binCount; j++) {
				Matrix<int> csIndxs(N, randomCount*N, false);
				conditionalSubstitutionIndxs1[i][j] = csIndxs;
				conditionalSubstitutionIndxs2[i][j] = csIndxs;
			}
		}
		initSubstitutionIndxs1.resize(N, randomCount*N, false);
		initSubstitutionIndxs2.resize(N, randomCount*N, false);
	
		baseQualityIndxs.resize(N, randomCount*baseQualtiyCount, false);
		errorBaseQualityIndxs.resize(N, randomCount*baseQualtiyCount, false);
		
		gcValues = new Matrix<double>[101];
	}
}

void Profile::evaluateGC() {
	int i, j, k;
	/*calculate median value*/
	memset(gcFactors, 0, 101*sizeof(double));
	map<int, vector<double> >::iterator it;
	double minFactor = -1;
	int minGC, maxGC;
	for(it = readCountsByGC.begin(); it != readCountsByGC.end(); it++) {
		int gc = (*it).first;
		vector<double>& readCounts = (*it).second;
		sort(readCounts.begin(), readCounts.end());

		if(readCounts.size() > 10) {
			int n = readCounts.size()/100;
		
			int m = readCounts.size()-2*n;
			if(m%2 == 0) {
				gcFactors[gc] = (readCounts[n+m/2-1]+readCounts[n+m/2])/2;
			}
			else {
				gcFactors[gc] = readCounts[n+m/2];
			}
		}
		else {
			gcFactors[gc] = 0;
		}

		if(it == readCountsByGC.begin()) {
			minGC = gc;
		}
		maxGC = gc;

		if(minFactor < 0 || minFactor > gcFactors[gc]) {
			minFactor = gcFactors[gc];
		}
	}
	//double minFactor = *min_element(gcFactors, gcFactors+101);
	//cerr << minFactor << endl;

	/*interpolation*/
	Matrix<double> B(maxGC-minGC+1, 2, false);
	Matrix<double> y(maxGC-minGC+1, 1, false);
	for(i = minGC; i <= maxGC; i++) {
		if(gcFactors[i] == 0) {
			if(i > 0) {
				double y1 = gcFactors[i-1];
				double y2;
				for(j = i+1; j < 101; j++) {
					if(gcFactors[j] > 0) {
						y2 = gcFactors[j];
						break;
					}
				}
				if(j == 101) {
					y2 = minFactor;
					j = 100;
				}
				gcFactors[i] = (y2-y1)/(j-i+1)+y1;
			}
			else {
				gcFactors[i] = minFactor;
			}
		}
		B.set(i-minGC, 0, 1);
		B.set(i-minGC, 1, i);
		y.set(i-minGC, 0, gcFactors[i]);
		//cerr << "gc=" << i << ", mean=" << gcFactors[i] << endl;
	}
	
	/*fitting locally weighted linear regression*/
	double tau = 5;
	for(i = minGC; i <= maxGC; i++) {
	//for(i = 0; i < 101; i++) {
		Matrix<double> W(B.getROWS(), B.getROWS(), true);
		for(j = minGC; j <= maxGC; j++) {
		//for(j = 0; j < 101; j++) {
			double v = exp(-pow((double)(i-j), 2)/(2*tau));
			W.set(j-minGC, j-minGC, v);
			//W.set(j, j, v);
		}
		Matrix<double> beta = (B.transpose()*W*B).inverse()*B.transpose()*W*y;
		Matrix<double> y_predict = B.Row(i-minGC)*beta;
		//Matrix<double> y_predict = B.Row(i)*beta;
		gcFactors[i] = max(0.0, y_predict.get(0, 0));
	}
	
	/*calculate standrad deviation*/
	for(i = 0; i < 101; i++) {
		vector<double> readCounts = readCountsByGC[i];
		if(readCounts.empty()) {
			gcStds[i] = 0;
		}
		else {
			int n = 10*readCounts.size()/100;
			gcStds[i] = 0;
			double avg = 0;
			for(j = n; j < readCounts.size()-n; j++) {
				avg += readCounts[j];
			}
			avg /= (readCounts.size()-2*n);
			for(j = n; j < readCounts.size()-n; j++) {
				gcStds[i] += pow(readCounts[j]-avg, 2);
			}
			gcStds[i] = sqrt(gcStds[i]/(readCounts.size()-2*n));
			gcStds[i] *= pow(gcFactors[i]/avg, 2);
		}
		//cerr << "gc=" << i << ", mean=" << gcFactors[i] << ", stdev=" << gcStds[i] << endl;
	}
	
	readCountsByGC.clear();
}

bool Profile::getNextLine(ifstream& ifs, string& line, int& lineNum) {
	line = "";
	while(getline(ifs, line)) {
		lineNum++;
		if(!line.empty() && line.at(0) != '#') {
			break;
		}
	}
	if(line.empty()) {
		return false;
	}
	return true;
}

void Profile::load(string modelFile) {
	ifstream ifs;
	ifs.open(modelFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open file " << modelFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	
	estimatedCoverage = -1;
	nucleotides = "";
	binCount = -1;

	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	
	string errMsg = "Error: malformed model file "+modelFile+" @line ";
	
	// parse "coverage", "nucleotides" and "binCount"
	while(getNextLine(ifs, line, lineNum)) {
		vector<string> fields = split(line, ':');
		if(fields.size() != 2) {
			cerr << errMsg << lineNum << "\n" << line;
			exit(1);
		}
		if(trim(fields[0]).compare("estimated sequencing coverage") == 0) {
			estimatedCoverage = atoi(trim(fields[1]).c_str());
			if(estimatedCoverage <= 0) {
				cerr << errMsg << lineNum << "\n" << line;
				exit(1);
			}
		}
		else if(trim(fields[0]).compare("nucleotides") == 0) {
			nucleotides = trim(fields[1]);
			if(nucleotides.empty()) {
				cerr << errMsg << lineNum << "\n" << line;
				exit(1);
			}
		}
		else if(trim(fields[0]).compare("binCount") == 0) {
			binCount = atoi(trim(fields[1]).c_str());
			if(binCount <= 0) {
				cerr << errMsg << lineNum << "\n" << line;
				exit(1);
			}
		}
		else {
			cerr << errMsg << lineNum << "\n" << line;
			exit(1);
		}
		if(estimatedCoverage > 0 && !nucleotides.empty() && binCount > 0) {
			break;
		}
	}
	if(estimatedCoverage < 0 || nucleotides.empty() || binCount < 0) {
		cerr << "Error: malformed model file " << modelFile << endl;
		exit(1);
	}

	int N = nucleotides.length();
	
	int i, j, k, l;
	int paraLoadedCount = 0;
	vector<string> fields;
	while(getNextLine(ifs, line, lineNum)) {
		if(line.compare("[Initial Substitution Probs]") == 0) {
			for(j = 0; j < N*2; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != N) {
						cerr << errMsg << lineNum << "\n" << line;
						exit(1);
					}
					for(k = 0; k < N; k++) {
						double prob = atof(trim(fields[k]).c_str());
						if(j < N) {
							initSubstitutionProbs1.set(j, k, prob);
						}
						else {
							initSubstitutionProbs2.set(j-N, k, prob);
						}
					}
				}
				else {
					cerr << "Error: malformed model file " << modelFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Conditional Substitution Probs]") == 0) {
			for(i = 0; i < N; i++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, ':');
					if(fields.size() != 2 || trim(fields[0]).compare("Nucleotide") != 0) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					char nucle = trim(fields[1]).at(0);
					int nucleIndx = nucleotides.find(nucle);
					if(nucleIndx == string::npos) {
						cerr << "Error: unrecognized nucleotide @line " << lineNum <<
								" in model file " << modelFile << endl;
						cerr << line << endl;
						exit(1);
					}
					for(l = 0; l < binCount; l++) {
						if(getNextLine(ifs, line, lineNum)) {
							fields = split(line, ':');
							if(fields.size() != 2 || trim(fields[0]).compare("binIndx") != 0) {
								cerr << errMsg << lineNum << "\n" << line << endl;
								exit(1);
							}
							int binIndx = atoi(trim(fields[1]).c_str());
							if(binIndx < 0 || binIndx >= binCount) {
								cerr << errMsg << lineNum << "\n" << line << endl;
								exit(1);
							}
							for(j = 0; j < N*2; j++) {
								if(getNextLine(ifs, line, lineNum)) {
									fields = split(line, '\t');
									if(fields.size() != N) {
										cerr << errMsg << lineNum << "\n" << line << endl;
										exit(1);
									}
									for(k = 0; k < N; k++) {
										double prob = atof(trim(fields[k]).c_str());
										if(j < N) {
											conditionalSubstitutionProbs1[nucleIndx][binIndx].set(j, k, prob);
										}
										else {
											conditionalSubstitutionProbs2[nucleIndx][binIndx].set(j-N, k, prob);
										}
									}
								}
								else {
									cerr << "Error: malformed model file " << modelFile << endl;
									exit(1);
								}
							}
						}
						else {
							cerr << "Error: malformed model file " << modelFile << endl;
							exit(1);
						}
					}
				}
				else {
					cerr << "Error: malformed model file " << modelFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Base Quality Distribution]") == 0) {
			for(j = 0; j < N; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != baseQualtiyCount) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					for(k = 0; k < baseQualtiyCount; k++) {
						double prob = atof(trim(fields[k]).c_str());
						baseQualityDist.set(j, k, prob);
					}
				}
				else {
					cerr << "Error: malformed model file " << modelFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Error Base Quality Distribution]") == 0) {
			for(j = 0; j < N; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != baseQualtiyCount) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					for(k = 0; k < baseQualtiyCount; k++) {
						double prob = atof(trim(fields[k]).c_str());
						errorBaseQualityDist.set(j, k, prob);
					}
				}
				else {
					cerr << "Error: malformed model file " << modelFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Insert Size Standard Deviation]") == 0) {
			if(getNextLine(ifs, line, lineNum)) {
				if(line.empty()) {
					cerr << errMsg << lineNum << "\n" << line << endl;
					exit(1);
				}
				stdInsertSize = atof(trim(line).c_str());
			}
			else {
				cerr << "Error: malformed model file " << modelFile << endl;
				exit(1);
			}
			paraLoadedCount++;
		}
		else if(line.compare("[GC-content Factors]") == 0) {
			for(j = 0; j < 101; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != 3) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					int gc = atoi(fields[0].c_str());
					double mean = atof(fields[1].c_str());
					double stdev = atof(fields[2].c_str());
					gcFactors[gc] = mean;
					gcStds[gc] = stdev;
				}
				else {
					cerr << "Error: malformed model file " << modelFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else {
			continue;
		}
	}
	ifs.close();
	
	if(paraLoadedCount < 6) {
		cerr << "Error: corrupted model file " << modelFile <<
				", failed to load some parameters!" << endl;
		exit(1);
	}
}

void Profile::saveResults(string bamFile, string outFile) {
	ostream* ost;
	ofstream ofs;
	if(!outFile.empty()) {	
		ofs.open(outFile.c_str());
		if(!ofs.is_open()) {
			cerr << "Error: cannot open file to save model training results:\n" << outFile << endl;
			exit(-1);
		}
		ost = &ofs;
	}
	else {
		ost = &(std::cout);
	}
	
	time_t timel;  
	time(&timel);     
	(*ost) << "#model created at " << asctime(gmtime(&timel));
	(*ost) << "#reads: " << bamFile << endl << endl;
	(*ost) << "estimated sequencing coverage: " << estimatedCoverage << endl;
	(*ost) << "nucleotides: " << nucleotides << endl;
	(*ost) << "binCount: " << binCount << endl << endl;
	
	int i, j, k, l;
	double **p;

	(*ost) << "[Initial Substitution Probs]" << endl;
	p = initSubstitutionProbs1.getEntrance();
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < nucleotides.length(); j++) {
			if(j < nucleotides.length()-1) {
				(*ost) << p[i][j] << '\t';
			}
			else {
				(*ost) << p[i][j] << endl;
			}
		}
	}
	p = initSubstitutionProbs2.getEntrance();
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < nucleotides.length(); j++) {
			if(j < nucleotides.length()-1) {
				(*ost) << p[i][j] << '\t';
			}
			else {
				(*ost) << p[i][j] << endl;
			}
		}
	}
	
	(*ost) << "\n[Conditional Substitution Probs]" << endl;
	for(i = 0; i < nucleotides.length(); i++) {
		(*ost) << "Nucleotide: " << nucleotides[i] << endl;
		for(j = 0; j < binCount; j++) {
			p = conditionalSubstitutionProbs1[i][j].getEntrance();
			(*ost) << "binIndx: " << j << endl;
			for(k = 0; k < nucleotides.length(); k++) {
				for(l = 0; l < nucleotides.length(); l++) {
					if(l < nucleotides.length()-1) {
						(*ost) << p[k][l] << '\t';
					}
					else {
						(*ost) << p[k][l] << endl;
					}
				}
			}
			p = conditionalSubstitutionProbs2[i][j].getEntrance();
			for(k = 0; k < nucleotides.length(); k++) {
				for(l = 0; l < nucleotides.length(); l++) {
					if(l < nucleotides.length()-1) {
						(*ost) << p[k][l] << '\t';
					}
					else {
						(*ost) << p[k][l] << endl;
					}
				}
			}
		}	
	}
	
	(*ost) << "\n[Base Quality Distribution]" << endl;
	p = baseQualityDist.getEntrance();
	k = baseQualityDist.getCOLS();
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < k; j++) {
			if(j < k-1) {
				(*ost) << p[i][j] << '\t';
			}
			else {
				(*ost) << p[i][j] << endl;
			}
		}
	}
	
	(*ost) << "\n[Error Base Quality Distribution]" << endl;
	p = errorBaseQualityDist.getEntrance();
	k = errorBaseQualityDist.getCOLS();
	for(i = 0; i < nucleotides.length(); i++) {
		for(j = 0; j < k; j++) {
			if(j < k-1) {
				(*ost) << p[i][j] << '\t';
			}
			else {
				(*ost) << p[i][j] << endl;
			}
		}
	}
	
	(*ost) << "\n[Insert Size Standard Deviation]" << endl;
	(*ost) << stdInsertSize << endl;
	
	(*ost) << "\n[GC-content Factors]" << endl;
	for(i = 0; i < 101; i++) {
		(*ost) << i << '\t' << gcFactors[i] << '\t' << gcStds[i] << endl;
	}
	
	if(!outFile.empty()) {
		ofs.close();
	}
}

void Profile::initSamplingPool() {	
	unsigned int i, j, k, l, iters = 10000;

	double zero_final = 1e-8;
	
	//baseQuality sampling
	Matrix<double> baseQualityCdf = baseQualityDist.cumsum();
	for(j = 0; j < nucleotides.length(); j++) {
		Matrix<int> tmp = randIndx(1, baseQualityIndxs.getCOLS(), baseQualityCdf.Row(j), true);
		baseQualityIndxs.setRow(j, tmp);
	}
	baseQualityDist.clear();
	baseQualityCdf.clear();
	
	//errorBaseQuality sampling
	Matrix<double> errorBaseQualityCdf = errorBaseQualityDist.cumsum();
	for(j = 0; j < nucleotides.length(); j++) {
		Matrix<int> tmp = randIndx(1, errorBaseQualityIndxs.getCOLS(), errorBaseQualityCdf.Row(j), true);
		errorBaseQualityIndxs.setRow(j, tmp);
	}
	errorBaseQualityDist.clear();
	errorBaseQualityCdf.clear();
	
	//insertSize sampling
	if(config.isPairedEnd() && stdInsertSize > 0) {
		Matrix<double> insertSizeCdf = insertSizeDist.cumsum();
		Matrix<int> tmp = randIndx(1, insertSizeIndxs.getCOLS(), insertSizeCdf, true);
		insertSizeIndxs.setRow(0, tmp);
		insertSizeCdf.clear();
		insertSizeDist.clear();
	}
	
	//gc sampling
	int coverage = config.getIntPara("coverage");
	//gcAlphabets = new Matrix<int>[101];
	for(i = 0; i < 101; i++) {
		/*
		if(estimatedCoverage > 0) {
			gcFactors[i] *= (double) coverage/estimatedCoverage;
			gcStds[i] *= pow((double) coverage/estimatedCoverage, 2);
		}
		*/

		/*
		double intervalLen = 1.96 * 2 * gcStds[i];
		int minValue = floor(max(gcFactors[i] - intervalLen/2, 0.0));
		int maxValue = minValue+intervalLen;
		//cerr << "gc=" << i << ", std=" << gcStds[i] << ", minValue=" << minValue << ", maxValue=" << maxValue << ", meanValue=" << gcFactors[i] << endl;
		gcAlphabets[i].resize(1, maxValue-minValue+1, false);
		for(j = 0; j < gcAlphabets[i].getCOLS(); j++) {
			gcAlphabets[i].set(0, j, minValue++);
		}
		gcDist[i].resize(1, gcAlphabets[i].getCOLS(), false);
		if(gcAlphabets[i].getCOLS() == 1) {
			gcIndxs[i].resize(1, 1, false);
			gcDist[i].set(0, 0 ,1);
		}
		else {
			gcIndxs[i].resize(1, randomCount*gcAlphabets[i].getCOLS(), false);
			for(j = 0; j < gcAlphabets[i].getCOLS(); j++) {
				double pdf = normpdf(gcAlphabets[i].get(0, j), gcFactors[i], gcStds[i]);
				gcDist[i].set(0, j, pdf);
			}
			gcDist[i].normalize(0);
		}
		*/
	}
	/*
	for(k = 0; k < 101; k++) {
		Matrix<double> gcCdf = gcDist[k].cumsum();
		Matrix<int> tmp = randIndx(1, gcIndxs[k].getCOLS(), gcCdf, true);
		gcIndxs[k].setRow(0, tmp);
		gcCdf.clear();
	}
	*/
	for(k = 0; k < 101; k++) {
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		normal_distribution<double> normal(gcFactors[k], gcStds[k]);
	 	j = max(1.0, sqrt(gcStds[k]));
		gcValues[k].resize(1, j*randomCount, false);
		for(i = 0; i < j*randomCount; i++) {
			double v = normal(generator);
			while(v < 0) {
				v = normal(generator);
			}
			gcValues[k].set(0, i, v);
		}
	}

	//cs and is sampling
	int N = nucleotides.length();
	Matrix<double> isc1 = initSubstitutionProbs1.cumsum();
	Matrix<double> isc2 = initSubstitutionProbs2.cumsum();
	for(i = 0; i < N; i++) {
		for(j = 0; j < binCount; j++) {
			Matrix<double> csc = conditionalSubstitutionProbs1[i][j].cumsum();
			for(k = 0; k < N; k++) {
				Matrix<int> tmp = randIndx(1, conditionalSubstitutionIndxs1[i][j].getCOLS(), csc.Row(k), true);
				conditionalSubstitutionIndxs1[i][j].setRow(k, tmp);
			}
			if(config.isPairedEnd()) {
				if(stdInsertSize > 0) {
					csc = conditionalSubstitutionProbs2[i][j].cumsum();
					for(k = 0; k < N; k++) {
						Matrix<int> tmp = randIndx(1, conditionalSubstitutionIndxs2[i][j].getCOLS(), csc.Row(k), true);
						conditionalSubstitutionIndxs2[i][j].setRow(k, tmp);
					}
				}
				else {
					conditionalSubstitutionIndxs2[i][j].clear();
				}
				
			}
			else {
				conditionalSubstitutionIndxs2[i][j].clear();
			}
			conditionalSubstitutionProbs1[i][j].clear();
			conditionalSubstitutionProbs2[i][j].clear();
			csc.clear();
		}
		delete[] conditionalSubstitutionProbs1[i];
		delete[] conditionalSubstitutionProbs2[i];
		
		Matrix<int> tmp = randIndx(1, initSubstitutionIndxs1.getCOLS(), isc1.Row(i), true);
		initSubstitutionIndxs1.setRow(i, tmp);
		if(config.isPairedEnd()) {
			if(stdInsertSize > 0) {
				tmp = randIndx(1, initSubstitutionIndxs2.getCOLS(), isc2.Row(i), true);
				initSubstitutionIndxs2.setRow(i, tmp);
			}
			else {
				initSubstitutionIndxs2.clear();
			}
		}
	}
	isc1.clear();
	isc2.clear();
	delete[] conditionalSubstitutionProbs1;
	delete[] conditionalSubstitutionProbs2;
	conditionalSubstitutionProbs1 = NULL;
	conditionalSubstitutionProbs2 = NULL;
	initSubstitutionProbs1.clear();
	initSubstitutionProbs2.clear();
}

void Profile::train(string modelFile) {
	load(modelFile);
	cerr << "profile was loaded from file " << modelFile << endl;
	initParas(true);
	initSamplingPool();
}

void Profile::train(string bamFile, string samtools, string outFile) {
	if(samtools.empty()) {
		samtools = "samtools";
	}
	string cmd = samtools + " view -F 0xD04 -q 20 " + bamFile;
	FILE* fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "cannot open BAM file " << bamFile << endl;
		exit(-1);
    	}
	
	char buf[2500000];
	long count = 0;
	int ret;
	while(fgets(buf , 2500000 , fp)) {
		buf[strlen(buf)-1] = '\0';
		ret = processRead(buf);
		if(ret == 2) {
			count++;			
			break;
		}
		else {
			count += ret;
		}
	}
	fclose(fp);

	int wgs = (genome.getInTargets().empty())? 1:0;
	long minReadsRequired = 1000000;
	if((wgs == 1 && count < minReadsRequired) || (wgs == 0 && count < 2*minReadsRequired)) {
		cerr << "\nWarning: no enough reads to evaluate GC-content effects!" << endl;		
		initGCFactors();
	}
	else {	
		estimateCoverage();
		evaluateGC();
	}
	initParas(false);
	saveResults(bamFile, outFile);
}

int Profile::yieldInsertSize() {
	if(insertSizeAlphabet.getEntrance() == NULL) {
		return config.getIntPara("insertSize");
	}
	int i = threadPool->randomInteger(0, insertSizeIndxs.getCOLS());
	i = insertSizeIndxs.get(0, i);
	return insertSizeAlphabet.get(0, i);
}

double Profile::getStdInsertSize() {
	return stdInsertSize;
}

int Profile::getMaxInsertSize() {
	if(insertSizeAlphabet.getEntrance() == NULL) {
		return config.getIntPara("insertSize");
	}
	int k = insertSizeAlphabet.getCOLS();
	return insertSizeAlphabet.get(0, k-1);
}

double Profile::getGCFactor(int gc) {
	if(gc < 0 || gc > 100) {
		return 0;
	}
	/*
	double weight = 0;
	for(int i = 0; i < gcDist[gc].getCOLS(); i++) {
		weight += gcDist[gc].get(0, i)*gcAlphabets[gc].get(0, i);
	}
	return weight;
	*/
	int i = threadPool->randomInteger(0, gcValues[gc].getCOLS());
	return gcValues[gc].get(0, i);
}

int Profile::getReplicates(int gc) {
	if(gc < 0 || gc > 100) {
		return 0;
	}
	int i = threadPool->randomInteger(0, gcIndxs[gc].getCOLS());
	i = gcIndxs[gc].get(0, i);
	return gcAlphabets[gc].get(0, i);
}

int Profile::getInitBaseIndx1(int base) {
	int j = base;
	int i = threadPool->randomInteger(0, initSubstitutionIndxs1.getCOLS());
	return initSubstitutionIndxs1.get(j, i);
}

int Profile::getInitBaseIndx2(int base) {
	int i, j = base;
	if(initSubstitutionIndxs2.getEntrance() == NULL) {
		i = threadPool->randomInteger(0, initSubstitutionIndxs1.getCOLS());
		return initSubstitutionIndxs1.get(j, i);
	}
	else {
		i = threadPool->randomInteger(0, initSubstitutionIndxs2.getCOLS());
		return initSubstitutionIndxs2.get(j, i);
	}
}

int Profile::getCondBaseIndx1(int preBase, int base, int binIndx) {
	int i = preBase, j = base;
	int k = threadPool->randomInteger(0, conditionalSubstitutionIndxs1[j][binIndx].getCOLS());
	return conditionalSubstitutionIndxs1[j][binIndx].get(i, k);
}

int Profile::getCondBaseIndx2(int preBase, int base, int binIndx) {
	int i = preBase, j = base, k;
	if(conditionalSubstitutionIndxs2[j][binIndx].getEntrance() == NULL) {
		k = threadPool->randomInteger(0, conditionalSubstitutionIndxs1[j][binIndx].getCOLS());
		return conditionalSubstitutionIndxs1[j][binIndx].get(i, k);
	}
	else {
		k = threadPool->randomInteger(0, conditionalSubstitutionIndxs2[j][binIndx].getCOLS());
		return conditionalSubstitutionIndxs2[j][binIndx].get(i, k);
	}
}

int Profile::getBaseQuality(int base) {
	int j = base;
	int i = threadPool->randomInteger(0, baseQualityIndxs.getCOLS());
	i = baseQualityIndxs.get(j, i);
	return qualityAlphabet.get(0, i);
}

int Profile::getErrorBaseQuality(int base) {
	int j = base;
	int i = threadPool->randomInteger(0, errorBaseQualityIndxs.getCOLS());
	i = errorBaseQualityIndxs.get(j, i);
	return qualityAlphabet.get(0, i);
}

void Profile::predict(char* refSeq, char* results, int num, int isRead1) {
	if(refSeq == NULL || refSeq[0] == '\0') {
		results[0] = '\0';
		return;
	}
	
	int n = strlen(refSeq);
	int nucleIndxs[n];
	char readSeq[n];
	char baseQuality[n];
	
	//int iters = 100;
	int iters = 0;
	int i, j, k, quality;
	int preIndx, indx, indx1, binIndx;
	
	for(i = 0; i < n; i++) {
		nucleIndxs[i] = getIndexOfNucleotide(refSeq[i]);
	}
	
	int m = 0;
	for(i = 0; i < num; i++) {
		for(j = 0; j < n; j++) {
			indx = getIndexOfNucleotide(refSeq[j]);
			//indx = nucleIndxs[j];
			//k = indx;
			if(j == 0) {
				if(isRead1) {
					k = getInitBaseIndx1(indx);
				}
				else {
					k = getInitBaseIndx2(indx);
				}
				nucleIndxs[0] = k;
				readSeq[0] = nucleotides[k];
			}
			else {
				binIndx = j*binCount/n;
				preIndx = getIndexOfNucleotide(refSeq[j-1]);
				if(isRead1) {
					//k = getCondBaseIndx1(nucleIndxs[j-1], indx, binIndx);
					k = getCondBaseIndx1(preIndx, indx, binIndx);
				}
				else {
					//k = getCondBaseIndx2(nucleIndxs[j-1], indx, binIndx);
					k = getCondBaseIndx2(preIndx, indx, binIndx);
				}
				nucleIndxs[j] = k;
				readSeq[j] = nucleotides[k];
			}
			
			if(indx == k) {
				quality = getBaseQuality(k);
			}
			else {
				quality = getErrorBaseQuality(k);
			}
			baseQuality[j] = quality;
		}
		strncpy(&results[m], readSeq, n);
		strncpy(&results[m+n], baseQuality, n);
		m += 2*n;
	}
	results[m] = '\0';
}
