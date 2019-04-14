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
#include "psiFunc.h"
#include "Profile.h"

Profile::Profile() {
	subsDist1 = subsDist2 = NULL;
	subsCdf1 = subsCdf2 = NULL;
	qualityDist = qualityCdf = NULL;
}

Profile::~Profile() {
	int i, j;
	
	string bases = config.getStringPara("bases");
	int N = bases.length();
	
	if(subsDist1 != NULL) {
		for(i = 0; i < kmerCount; i++) {
			subsDist1[i].clear();
			subsCdf1[i].clear();
		}
		delete[] subsDist1;
		delete[] subsCdf1;
	}
	if(subsDist2 != NULL) {
		for(i = 0; i < kmerCount; i++) {
			subsDist2[i].clear();
			subsCdf2[i].clear();
		}
		delete[] subsDist2;
		delete[] subsCdf2;
	}
	if(qualityDist != NULL) {
		for(i = 0; i < N*N; i++) {
			qualityDist[i].clear();
			qualityCdf[i].clear();
		}
		delete[] qualityDist;
		delete[] qualityCdf;
	}
	
	if(kmers != NULL) {
		for(i = 0;i < kmerCount; i++) {
			delete[] kmers[i];
		}
		delete[] kmers;
	}
	
}

void Profile::initKmers() {
	string bases = config.getStringPara("bases");
	int kmer = config.getIntPara("kmer");
	int N = bases.length();
	int *tmp = new int[kmer];
	int i, j, k, n;
	
	kmerCount = 0;
	for(i = kmer-1; i >= 0; i--) {
		kmerCount += pow(N, kmer-i);
	}
	kmers = new char*[kmerCount];
	for(i = 0; i < kmerCount; i++) {
		kmers[i] = new char[kmer+1];
		kmers[i][kmer] = '\0';
	}
	
	k = 0;
	for(j = kmer-1; j >= 0; j--) {
		for(i = j; i < kmer; i++) {
			tmp[i] = 0;
		}
		while(tmp[j] < N) {
			KmerIndex *sIndx = &rIndex;
			for(i = 0; i < j; i++) {
				kmers[k][i] = 'X';
				sIndx = &(sIndx->nextIndexs['X']);
			}
			for(; i < kmer; i++) {
				char c = bases[tmp[i]];
				kmers[k][i] = c;
				sIndx = &(sIndx->nextIndexs[c]);
			}
			sIndx->index = k;
			k++;
			n = 1;
			for(i = kmer-1; i > j; i--) {
				if(n == 0) {
					break;
				}
				tmp[i] += n;
				if(tmp[i] == N) {
					tmp[i] = 0;
					n = 1;
				}
				else {
					n = 0;
				}
			}
			tmp[j] += n;
		}
	}
	
	delete[] tmp;
}

void Profile::init() {
	minBaseQuality = 33;
	maxBaseQuality = 126;
	
	indelRate = config.getRealPara("indelRate");
	
	int binCount = config.getIntPara("bins");
	
	initKmers();
	
	string bases = config.getStringPara("bases");
	int N = bases.length();
	
	int i, j;
	kmersDist.resize(binCount, kmerCount, 0);
	
	subsDist1 = new Matrix<double>[kmerCount];
	subsDist2 = new Matrix<double>[kmerCount];
	subsCdf1 = new Matrix<double>[kmerCount];
	subsCdf2 = new Matrix<double>[kmerCount];
	Matrix<double> tmp(binCount, N, 0);
	for(i = 0; i < kmerCount; i++) {
		subsDist1[i] = tmp;
		subsDist2[i] = tmp;
	}
	
	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	qualityDist = new Matrix<double>[N*N];
	qualityCdf = new Matrix<double>[N*N];
	Matrix<double> tmp1(binCount, baseQualtiyCount, 0);
	for(i = 0; i < N*N; i++) {
		qualityDist[i] = tmp1;
	}
	
	iSizeDist.resize(1, 1000, true);
	stdISize = 0;
}

int Profile::getKmerIndx(const char *s) {
	KmerIndex *sIndx = &rIndex;
	for(int i = 0; s[i] != '\0'; i++) {
		sIndx = &(sIndx->nextIndexs[s[i]]);
	}
	return sIndx->index;
}

int Profile::processRead(char* read) {
	if(read == NULL || read[0] == '\0') {
		return 0;
	}

	static int wxs = (genome.getInTargets().empty())? 0:1;
	
	static unsigned long readCount = 0;
	static unsigned long maxCount = 300000000;
	static vector<string>& chromosomes = genome.getChroms();
	vector<string>::iterator v_it;
	static int binCount = config.getIntPara("bins");
	
	int i, j, k, n;
	
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
	
	/***update subsDist***/
	int kmerIndx, refIndx, baseIndx, binIndx;
	static int kmer = config.getIntPara("kmer");
	n = strlen(refSeq);
	
	char seq[kmer+n];
	for(i = 0; i < kmer-1; i++) {
		seq[i] = 'X';
	}
	for(; i < kmer+n-1; i++) {
		if(altSeq[i-kmer+1] == readSeq[i-kmer+1]) {
			seq[i] = altSeq[i-kmer+1];
		}
		else {
			seq[i] = refSeq[i-kmer+1];
		}
	}
	for(i = 0; i < n; i++) {
		baseIndx = getIndexOfBase(readSeq[i]);
		binIndx = i*binCount/n;
		if(baseIndx != -1) {
			
			char c = seq[i+kmer];
			seq[i+kmer] = '\0';
			kmerIndx = getKmerIndx(&seq[i]);
			seq[i+kmer] = c;
			
			if(kmerIndx == -1) {
				continue;
			}
			
			if(isRead1) {
				k = subsDist1[kmerIndx].get(binIndx, baseIndx);
				subsDist1[kmerIndx].set(binIndx, baseIndx, k+1);
			}
			else {
				k = subsDist2[kmerIndx].get(binIndx, baseIndx);
				subsDist2[kmerIndx].set(binIndx, baseIndx, k+1);
			}
			
			k = kmersDist.get(binIndx, kmerIndx);
			kmersDist.set(binIndx, kmerIndx, k+1);
		}
	}
	
	/***update insert sizeset distribution***/
	//if(config.isPairedEnd() && tlen > 0) {
	if(tlen > 0) {
		if(tlen > iSizeDist.getCOLS()-1) {
			iSizeDist.resize(1, tlen+1, true);
		}
		iSizeDist.set(0, tlen, iSizeDist.get(0, tlen)+1);
	}
	
	int indx;
	static string bases = config.getStringPara("bases");
	int N = bases.length();
	/***update base quality distribution***/
	if(strlen(baseQuality) == strlen(readSeq)) {
		n = strlen(readSeq);
		for(i = 0; i < n; i++) {
			refIndx = getIndexOfBase(refSeq[i]);
			binIndx = i*binCount/n;
			baseIndx = getIndexOfBase(readSeq[i]);
			if(refIndx == -1 || baseIndx == -1) {
				continue;
			}
			if(altSeq[i] == readSeq[i]) {
				refIndx = getIndexOfBase(altSeq[i]);
			}
			indx = refIndx*N+baseIndx;
			
			j = baseQuality[i];
			//int j = baseQuality[i]-33;
			
			if(j >= minBaseQuality && j <= maxBaseQuality) {
				//qualityData[indx][binIndx].push_back(j);
				j -= minBaseQuality;
				k = qualityDist[indx].get(binIndx, j);
				qualityDist[indx].set(binIndx, j, k+1);
			}
		}
	}
	
	readCount++;
	
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
	
	if(wxs == 0 && readCount >= maxCount) {
		return 2;
	}
	/*
	else {
		return 0;
	}
	*/
	if(wxs == 1 && readCount >= 2*maxCount) {
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
	static double GC = -1;
	static int rc = 0;
	static int wxs = (genome.getInTargets().empty())? 0:1;
	
	static vector<string>& chromosomes = genome.getChroms(); 

	static int targetIndx = -1;
	
	char* refSeq;

	position -= 1;
	
	int i, k;
	if(chr.compare("X") == 0 || chr.compare("Y") == 0 || chr.compare("M") == 0)
	{
		return -1;
	}
	/*
	for(i = 0; i < chr.length(); i++) {
		k = chr[i]-'0';
		if(!(k >= 0 && k <= 9)) {
			break;
		}
	}
	if(i < chr.length()) {
		return -1;
	}
	*/
	
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
		if(GC > 0 && rc > 0) {
			if(wxs == 0) {
				gcs.push_back(GC);
				readCounts.push_back(rc);
			}
			else {
				int targetSize = rightPos-leftPos+1;
				rc = winSize*rc/targetSize;
				gcs.push_back(GC);
				readCounts.push_back(rc);
			}
		}
		if(wxs == 0) {
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
			GC = calculateGCContent(refSeq);
			delete[] refSeq;
		}
		else {
			GC = -1;
		}
		return rc;
	}
	else {
		if(preChr.compare("") != 0 && GC > 0 && rc > 0) {
			if(wxs == 0) {
				gcs.push_back(GC);
				readCounts.push_back(rc);
			}
			else {
				int targetSize = rightPos-leftPos+1;
				rc = winSize*rc/targetSize;
				gcs.push_back(GC);
				readCounts.push_back(rc);
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
			if(wxs == 0) {
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
			GC = calculateGCContent(refSeq);
			delete[] refSeq;
		}
		else {
			GC = -1;
		}
		return rc;
	}
	
}

void Profile::initGCParas() {
	int i;
	for(i = 0; i < 101; i++) {
		gcMeans[i] = 1;
	}
	gcStd = 1.0e-5;
}

void Profile::estimateGCParas(string outFile) {
	int i, j, k;
	
	int bins = 50;
	int *counts = new int[bins];
	memset(counts, 0, bins*sizeof(int));
	for(i = 0; i < gcs.size(); i++) {
		j = gcs[i]*bins;
		counts[j]++;
	}
	int expectCount = min(150000, (int) gcs.size())/bins;
	int *steps = new int[bins];
	for(i = 0; i < bins; i++) {
		steps[i] = max(1, counts[i]/expectCount);
	}
	delete[] counts;
	
	outFile = outFile + ".gc";
	ofstream ofs;
	ofs.open(outFile.c_str());
	vector<int> indxs;
	int *curCount = new int[bins];
	double med_rc = median(readCounts);
	for(i = 0; i < readCounts.size(); i++) {
		j = gcs[i]*bins;
		if(curCount[j]%steps[j] == 0) {
			//readCounts[i] = log(readCounts[i]/med_rc+ZERO_FINAL);
			readCounts[i] = readCounts[i]/(med_rc+ZERO_FINAL);
			if(readCounts[i] < 3) {
				ofs << readCounts[i] << '\t' << gcs[i] << endl;
				indxs.push_back(i);
			}
		}
		curCount[j]++;
	}
	delete[] curCount;
	delete[] steps;
	ofs.close();
	
	/*fitting locally weighted linear regression*/
	double tau = 5;
	double winSize = 0.03;
	Matrix<double> x(1, 2);
	x.set(0, 0, 1);
	int minGC = -1, maxGC = -1;
	for(k = 0; k <= 100; k++) {
		double gc = k/100.0;
		x.set(0, 1, gc);
		
		vector<double> gcInWin, rcInWin;
		for(i = 0; i < indxs.size(); i++) {
			j = indxs[i];
			if(fabs(gc-gcs[j]) <= winSize/2) {
				gcInWin.push_back(gcs[j]);
				rcInWin.push_back(readCounts[j]);
			}
			/*
			if(gcInWin.size() > 2000) {
				break;
			}
			*/
		}
		if(gcInWin.size() > 20) {
			if(minGC == -1) {
				minGC = k;
			}
			maxGC = k;
			Matrix<double> B(gcInWin.size(), 2);
			Matrix<double> y(gcInWin.size(), 1);
			for(i = 0; i < gcInWin.size(); i++) {
				B.set(i, 0, 1);
				B.set(i, 1, gcInWin[i]);
				y.set(i, 0, rcInWin[i]);
			}
			gcInWin.clear();
			rcInWin.clear();
			/*
			cerr << "B:" << endl;
			B.Print();
			cerr << "y:" << endl;
			y.Print();
			*/
			Matrix<double> W(B.getROWS(), B.getROWS(), 0);
			for(i = 0; i < B.getROWS(); i++) {
				double v = exp(-pow(B.get(i, 1)-gc, 2)/(2*tau));
				W.set(i, i, v);
			}
			//cerr << "W:" << endl;
			//W.Print();
			
			Matrix<double> beta = (B.transpose()*W*B).inverse()*B.transpose()*W*y;
			//cerr << "beta:" << endl;
			//beta.Print();
			Matrix<double> y_predict = x*beta;
			gcMeans[k] = max(0.0, y_predict.get(0, 0));
		}
		else {
			gcMeans[k] = 0;
		}
	}
	
	for(k = 0; k < minGC; k++) {
		gcMeans[k] = gcMeans[minGC]*k/minGC;
	}
	for(k = maxGC+1; k <= 100; k++) {
		gcMeans[k] = gcMeans[maxGC]-gcMeans[maxGC]*(k-maxGC)/(100-maxGC);
	}
	
	/*calculate standrad deviation*/
	gcStd = 0;
	for(i = 0; i < indxs.size(); i++) {
		j = indxs[i];
		k = gcs[j]*100;
		gcStd += pow(readCounts[j]-gcMeans[k], 2);
	}
	gcStd = sqrt(gcStd/indxs.size());
	cerr << "read counts std: " << gcStd << endl;
	
	gcs.clear();
	readCounts.clear();
}

void Profile::normParas(bool isLoaded) {
	int i, j, k;
	string bases = config.getStringPara("bases");
	int N = bases.length();
	int kmer = config.getIntPara("kmer");
	int binCount = config.getIntPara("bins");
	/***normalize probability matrix***/
	kmersDist.normalize(0);
	for(i = 0; i < kmerCount; i++) {
		subsDist1[i].normalize(0);
		int indx = getIndexOfBase(kmers[i][kmer-1]);
		Matrix<double> tmp = subsDist1[i].sumCols();
		for(j = 0; j < binCount; j++) {
			if(tmp.get(0, j) < ZERO_FINAL) {
				subsDist1[i].set(j, indx, 1);
			}
		}
		
		subsDist2[i].normalize(0);
		tmp = subsDist2[i].sumCols();
		for(j = 0; j < binCount; j++) {
			if(tmp.get(0, j) < ZERO_FINAL) {
				subsDist2[i].set(j, indx, 1);
			}
		}
	}
	
	for(i = 0; i < N*N; i++) {
		qualityDist[i].normalize(0);
	}
	
	if(!isLoaded) {
		int maxCount = 0;
		j = 0;
		for(i = 0; i < iSizeDist.getCOLS(); i++) {
			if(iSizeDist.get(0, i) > maxCount) {
				maxCount = iSizeDist.get(0, i);
				j = i;
			}
		}
		
		for(i = j*5; i < iSizeDist.getCOLS(); i++) {
			iSizeDist.set(0, i, 0);
		}

		iSizeDist.normalize(0);
		double meanTlen = 0;
		for(i = 0; i < iSizeDist.getCOLS(); i++) {
			meanTlen += iSizeDist.get(0, i)*i;
		}
		stdISize = 0;
		for(i = 0; i < iSizeDist.getCOLS(); i++) {
			stdISize += iSizeDist.get(0, i)*pow(i-meanTlen, 2);
		}
		stdISize = sqrt(stdISize);
		//cerr << "stdISize: " << stdISize << endl;
		//iSizeDist.clear();
	}
	else {
		baseAlphabet.resize(1, N, false);
		for(i = 0; i < N; i++) {
			baseAlphabet.set(0, i, i);
		}
		int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
		qualityAlphabet.resize(1, baseQualtiyCount, false);
		for(i = 0, j = minBaseQuality; i < baseQualtiyCount; i++, j++) {
			qualityAlphabet.set(0, i, j);
		}
	
		if(config.isPairedEnd() && stdISize > 0) {
			int meanInsertSize = config.getIntPara("insertSize")+1;
			int intervalLen = 6 * stdISize;
			int minInsertSize = max(meanInsertSize - intervalLen/2, config.getIntPara("readLength"));
			int maxInsertSize = 2*meanInsertSize - minInsertSize;
			//cerr << stdISize << endl;
			//cerr << minInsertSize << '\t' << maxInsertSize << '\t' << meanInsertSize << endl;
			iSizeAlphabet.resize(1, maxInsertSize-minInsertSize+1, false);
			for(i = 0; i < iSizeAlphabet.getCOLS(); i++) {
				iSizeAlphabet.set(0, i, minInsertSize++);
			}
		
			iSizeDist.resize(1, iSizeAlphabet.getCOLS(), false);
			for(i = 0; i < iSizeDist.getCOLS(); i++) {
				double pdf = normpdf(iSizeAlphabet.get(0, i), meanInsertSize, stdISize);
				iSizeDist.set(0, i, pdf);
			}
			iSizeDist.normalize(0);
		}
	}
}

void Profile::load(string proFile) {
	ifstream ifs;
	ifs.open(proFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open file " << proFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	
	string bases = "";
	int binCount = -1;
	int kmer = -1;
	
	string errMsg = "Error: malformed model file "+proFile+" @line ";
	
	// parse "coverage", "bases" and "binCount"
	while(getNextLine(ifs, line, lineNum)) {
		vector<string> fields = split(line, ':');
		if(fields.size() != 2) {
			cerr << errMsg << lineNum << "\n" << line;
			exit(1);
		}
		if(trim(fields[0]).compare("bases") == 0) {
			bases = trim(fields[1]);
			if(bases.empty()) {
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
		else if(trim(fields[0]).compare("kmer") == 0) {
			kmer = atoi(trim(fields[1]).c_str());
			if(kmer <= 0) {
				cerr << errMsg << lineNum << "\n" << line;
				exit(1);
			}
		}
		else {
			cerr << errMsg << lineNum << "\n" << line;
			exit(1);
		}
		if(!bases.empty() && binCount > 0 && kmer > 0) {
			break;
		}
	}
	if(bases.empty() || binCount <= 0 || kmer <= 0) {
		cerr << "Error: malformed model file " << proFile << endl;
		exit(1);
	}
	
	config.setStringPara("bases", bases);
	config.setIntPara("kmer", kmer);
	config.setIntPara("bins", binCount);
	
	init();

	int i, j, k, l;
	int N = bases.length();
	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	
	int paraLoadedCount = 0;
	vector<string> fields;
	while(getNextLine(ifs, line, lineNum)) {
		if(line.compare("[Kmer Distribution]") == 0) {
			for(j = 0; j < binCount; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != kmerCount) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					for(k = 0; k < kmerCount; k++) {
						double prob = atof(trim(fields[k]).c_str());
						kmersDist.set(j, k, prob);
					}
				}
				else {
					cerr << "Error: malformed profile file " << proFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Substitution Probs]") == 0) {
			for(i = 0; i < kmerCount; i++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, ':');
					if(fields.size() != 2 || trim(fields[0]).compare("kmer") != 0) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					fields[1] = trim(fields[1]);
					int kmerIndx = getKmerIndx(fields[1].c_str());
					if(kmerIndx == -1) {
						cerr << "Error: unrecognized kmer @line " << lineNum <<
								" in profile file " << proFile << endl;
						cerr << line << endl;
						exit(1);
					}
					for(j = 0; j < binCount*2; j++) {
						if(getNextLine(ifs, line, lineNum)) {
							fields = split(line, '\t');
							if(fields.size() != N) {
								cerr << errMsg << lineNum << "\n" << line << endl;
								exit(1);
							}
							for(k = 0; k < N; k++) {
								double prob = atof(trim(fields[k]).c_str());
								if(j < binCount) {
									subsDist1[kmerIndx].set(j, k, prob);
								}
								else {
									subsDist2[kmerIndx].set(j-binCount, k, prob);
								}
							}
						}
						else {
							cerr << "Error: malformed profile file " << proFile << endl;
							exit(1);
						}
					}
				}
				else {
					cerr << "Error: malformed profile file " << proFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Base Quality Distribution]") == 0) {
			for(i = 0; i < N*N; i++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, ':');
					if(fields.size() != 2 || trim(fields[0]).compare("basePairIndx") != 0) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					fields[1] = trim(fields[1]);
					int basePairIndx = atoi(fields[1].c_str());
					if(basePairIndx < 0 || basePairIndx > N*N-1) {
						cerr << "Error: unrecognized basePairIndx @line " << lineNum <<
								" in profile file " << proFile << endl;
						cerr << line << endl;
						exit(1);
					}
					for(j = 0; j < binCount; j++) {
						if(getNextLine(ifs, line, lineNum)) {
							fields = split(line, '\t');
							if(fields.size() != baseQualtiyCount) {
								cerr << errMsg << lineNum << "\n" << line << endl;
								exit(1);
							}
							for(k = 0; k < baseQualtiyCount; k++) {
								double prob = atof(trim(fields[k]).c_str());
								qualityDist[basePairIndx].set(j, k, prob);
							}
						}
						else {
							cerr << "Error: malformed profile file " << proFile << endl;
							exit(1);
						}
					}
				}
				else {
					cerr << "Error: malformed profile file " << proFile << endl;
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
				stdISize = atof(trim(line).c_str());
			}
			else {
				cerr << "Error: malformed model file " << proFile << endl;
				exit(1);
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Log Ratio Mean Value]") == 0) {
			for(j = 0; j < 101; j++) {
				if(getNextLine(ifs, line, lineNum)) {
					fields = split(line, '\t');
					if(fields.size() != 2) {
						cerr << errMsg << lineNum << "\n" << line << endl;
						exit(1);
					}
					int gc = atoi(fields[0].c_str());
					double mean = atof(fields[1].c_str());
					gcMeans[gc] = mean;
				}
				else {
					cerr << "Error: malformed model file " << proFile << endl;
					exit(1);
				}
			}
			paraLoadedCount++;
		}
		else if(line.compare("[Log Ratio Standard Deviation]") == 0) {
			if(getNextLine(ifs, line, lineNum)) {
				if(line.empty()) {
					cerr << errMsg << lineNum << "\n" << line << endl;
					exit(1);
				}
				gcStd = atof(trim(line).c_str());
			}
			else {
				cerr << "Error: malformed model file " << proFile << endl;
				exit(1);
			}
			paraLoadedCount++;
		}
		else {
			continue;
		}
	}
	ifs.close();
	
	if(paraLoadedCount < 6) {
		cerr << "Error: corrupted model file " << proFile <<
				", failed to load some parameters!" << endl;
		exit(1);
	}
	cerr << "profile was loaded from file " << proFile << endl;
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
	
	string bases = config.getStringPara("bases");
	int N = bases.length();
	int kmer = config.getIntPara("kmer");
	int binCount = config.getIntPara("bins");
	
	time_t timel;  
	time(&timel);     
	(*ost) << "#model created at " << asctime(gmtime(&timel));
	(*ost) << "#reads: " << bamFile << endl << endl;
	(*ost) << "bases: " << bases << endl;
	(*ost) << "binCount: " << binCount << endl;
	(*ost) << "kmer: " << kmer << endl << endl;
	
	int i, j, k, l;
	double *p;
	
	(*ost) << "\n[Kmer Distribution]" << endl;
	p = kmersDist.getEntrance();
	for(j = 0; j < binCount; j++) {
		for(k = 0; k < kmerCount; k++)  {
			if(k < kmerCount-1) {
				(*ost) << p[j*kmerCount+k] << '\t';
			}
			else {
				(*ost) << p[j*kmerCount+k] << endl;
			}
		}
	}
	
	(*ost) << "\n[Substitution Probs]" << endl;
	for(i = 0; i < kmerCount; i++) {
		(*ost) << "kmer: " << kmers[i] << endl;
		p = subsDist1[i].getEntrance();
		for(j = 0; j < binCount; j++) {
			for(k = 0; k < N; k++)  {
				if(k < N-1) {
					(*ost) << p[j*N+k] << '\t';
				}
				else {
					(*ost) << p[j*N+k] << endl;
				}
			}
		}
		p = subsDist2[i].getEntrance();
		for(j = 0; j < binCount; j++) {
			for(k = 0; k < N; k++)  {
				if(k < N-1) {
					(*ost) << p[j*N+k] << '\t';
				}
				else {
					(*ost) << p[j*N+k] << endl;
				}
			}
		}
	}
	
	(*ost) << "\n[Base Quality Distribution]" << endl;
	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	for(i = 0; i < N*N; i++) {
		(*ost) << "basePairIndx: " << i << endl;
		p = qualityDist[i].getEntrance();
		for(j = 0; j < binCount; j++) {
			for(k = 0; k < baseQualtiyCount; k++)  {
				if(k < baseQualtiyCount-1) {
					(*ost) << p[j*baseQualtiyCount+k] << '\t';
				}
				else {
					(*ost) << p[j*baseQualtiyCount+k] << endl;
				}
			}
		}
	}
	
	(*ost) << "\n[Insert Size Standard Deviation]" << endl;
	(*ost) << stdISize << endl;
	
	(*ost) << "\n[Log Ratio Mean Value]" << endl;
	for(i = 0; i < 101; i++) {
		(*ost) << i << '\t' << gcMeans[i] << endl;
	}
	(*ost) << "\n[Log Ratio Standard Deviation]" << endl;
	(*ost) << gcStd << endl;
	
	if(!outFile.empty()) {
		ofs.close();
	}
}

void Profile::initCDFs() {	
	unsigned int i, j, l;
	
	string bases = config.getStringPara("bases");
	int N = bases.length();
	int binCount = config.getIntPara("bins");
	
	//baseQuality
	int quality_th = 20;
	int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	for(i = 0; i < N*N; i++) {
		//qualityDist[i].normalize(0);
		/*
		if(i%(N+1) == 0) {
			for(j = 0; j < quality_th; j++) {
				qualityDist[i].setCol(j, 0.0);
			}
		}
		else {
			for(j = quality_th; j < baseQualtiyCount; j++) {
				qualityDist[i].setCol(j, 0.0);
			}
		}
		*/
		qualityDist[i].normalize(0);
		qualityCdf[i] = qualityDist[i].cumsum();
		qualityDist[i].clear();
	}
	
	//insertSize
	if(config.isPairedEnd() && stdISize > 0) {
		iSizeCdf = iSizeDist.cumsum();
		iSizeDist.clear();
	}
	
	//gc
	for(l = 0; l < 101; l++) {
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		normal_distribution<double> normal(gcMeans[l], gcStd);
		gc_generators.push_back(generator);
		gc_normDists.push_back(normal);
	}

	//subs
	for(i = 0; i < kmerCount; i++) {
		subsCdf1[i] = subsDist1[i].cumsum();
		if(config.isPairedEnd()) {
			if(stdISize > 0) {
				subsCdf2[i] = subsDist2[i].cumsum();
			}
			else {
				subsCdf2[i].clear();
			}
		}
		else {
			subsCdf2[i].clear();
		}
		subsDist1[i].clear();
		subsDist2[i].clear();
	}
}

void Profile::train(string proFile) {
	load(proFile);
	normParas(true);
	initCDFs();
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

	int wxs = (genome.getInTargets().empty())? 0:1;
	long minReadsRequired = 2000000;
	double med_rc = median(readCounts);
	// if((wxs == 0 && count < minReadsRequired) || (wxs == 1 && count < 2*minReadsRequired)) {
	if(med_rc < 5) {
		cerr << "\nWarning: no enough reads to evaluate GC-content effects!" << endl;		
		initGCParas();
	}
	else {	
		estimateGCParas(outFile);
	}
	normParas(false);
	saveResults(bamFile, outFile);
}

int Profile::yieldInsertSize() {
	if(iSizeAlphabet.getEntrance() == NULL) {
		return config.getIntPara("insertSize");
	}
	
	int k = randIndx(iSizeCdf.getEntrance(), iSizeAlphabet.getCOLS());
	return iSizeAlphabet.get(0, k);
}

double Profile::getStdISize() {
	return stdISize;
}

int Profile::getMaxInsertSize() {
	if(iSizeAlphabet.getEntrance() == NULL) {
		return config.getIntPara("insertSize");
	}
	int k = iSizeAlphabet.getCOLS();
	return iSizeAlphabet.get(0, k-1);
}

double Profile::getGCFactor(int gc) {
	if(gc < 0 || gc > 100) {
		return 0;
	}
	
	double v = gc_normDists[gc](gc_generators[gc]);
	while(v < 0) {
		v = gc_normDists[gc](gc_generators[gc]);
	}
	return v;
}

int Profile::getSubBaseIndx1(char *kmerSeq, int binIndx) {
	static int N = config.getStringPara("bases").length();
	static int kmer = config.getIntPara("kmer");
	int kmerIndx = getKmerIndx(kmerSeq);
	if(kmerIndx == -1) {
		return getIndexOfBase(kmerSeq[kmer-1]);
	}
	int k = randIndx(subsCdf1[kmerIndx].getEntrance()+binIndx*N, N);
	return k;
}

int Profile::getSubBaseIndx2(char *kmerSeq, int binIndx) {
	static int N = config.getStringPara("bases").length();
	static int kmer = config.getIntPara("kmer");
	int k;
	int kmerIndx = getKmerIndx(kmerSeq);
	if(kmerIndx == -1) {
		return getIndexOfBase(kmerSeq[kmer-1]);
	}
	if(subsCdf2[kmerIndx].getEntrance() == NULL) {
		int k = randIndx(subsCdf1[kmerIndx].getEntrance()+binIndx*N, N);
		return k;
	}
	else {
		int k = randIndx(subsCdf2[kmerIndx].getEntrance()+binIndx*N, N);
		return k;
	}
}

int Profile::getIndelSeq(vector<int> &baseIndxs) {
	static int N = config.getStringPara("bases").length();
	baseIndxs.clear();
	double p = threadPool->randomDouble(0, 1);
	if(p <= indelRate) {
		int i, k;
		int n = threadPool->randomInteger(1, 51);
		if(threadPool->randomDouble(0, 1) > 0.5) {
			for(i = 0; i < n; i++) {
				k = threadPool->randomInteger(0, N-1);
				baseIndxs.push_back(k);
			}
		}
		return n;
	}
	return 0;
}

int Profile::getBaseQuality(int basePairIndx, int binIndx) {
	static int baseQualtiyCount = maxBaseQuality-minBaseQuality+1;
	int i = randIndx(qualityCdf[basePairIndx].getEntrance()+binIndx*baseQualtiyCount, baseQualtiyCount);
	return qualityAlphabet.get(0, i);
}

int Profile::getRandBaseQuality() {
	return threadPool->randomInteger(minBaseQuality, minBaseQuality+20);
}

char* Profile::predict(char* refSeq, int isRead1) {
	if(refSeq == NULL || refSeq[0] == '\0') {
		return NULL;
	}
	
	int i, j, k;
	
	static string bases = config.getStringPara("bases");
	static int kmer = config.getIntPara("kmer");
	static int binCount = config.getIntPara("bins");
	static int N = bases.length();
	int n = strlen(refSeq);
	
	/*
	for(i = 0; i < n && refSeq[i] == 'N'; i++);
	for(j = i; j < n && refSeq[j] != 'N'; j++);
	
	refSeq = &refSeq[i];
	n = j-i;
	*/
	
	map<int, vector<int> > indelBaseIndxs;
	vector<int> indelLens;
	int indelLength = 0;
	for(j = 0; j < n;) {
		k = getIndelSeq(indelBaseIndxs[j]);
		if(indelBaseIndxs[j].empty() && k > 0) {
			k = min(n-j, k);
			indelLength -= k;
			indelLens.push_back(k);
			for(i = 1; i < k; i++) {
				indelLens.push_back(0);
			}
			j += k;
		}
		else {
			indelLength += k;
			j++;
			indelLens.push_back(k);
		}
	}
	if(n+indelLength < 50) {
		indelLength = 0;
		indelBaseIndxs.clear();
		indelLens.clear();
		for(j = 0; j < n; j++) {
			indelLens.push_back(0);
		}
	}
	
	char sourceSeq[n+indelLength+1];
	int m = 0;
	for(j = 0; j < n;) {
		if(indelBaseIndxs[j].empty() && indelLens[j] > 0) {	//del
			j += indelLens[j];
			continue;
		}
		else if(indelLens[j] == 0) {
			sourceSeq[m++] = refSeq[j];
			j++;
		}
		else {	//insert
			sourceSeq[m++] = refSeq[j];
			vector<int> &indxs = indelBaseIndxs[j];
			k = indxs.size();
			for(i = 0; i < k; i++) {
				sourceSeq[m++] = bases[indxs[i]];
			}
			j++;
		}
	}
	sourceSeq[n+indelLength] = '\0';
	n += indelLength;
	
	char seq[kmer+n];
	for(i = 0; i < kmer-1; i++) {
		seq[i] = 'X';
	}
	for(; i < kmer+n-1; i++) {
		seq[i] = sourceSeq[i-kmer+1];
	}
	
	int refIndx, basePairIndx, binIndx;
	char *results = new char[2*n+1];
	for(j = 0; j < n; j++) {
		refIndx = getIndexOfBase(sourceSeq[j]);
		binIndx = j*binCount/n;
		
		char c = seq[j+kmer];
		seq[j+kmer] = '\0';
		if(isRead1) {
			k = getSubBaseIndx1(&seq[j], binIndx);
		}
		else {
			k = getSubBaseIndx2(&seq[j], binIndx);
		}
		seq[j+kmer] = c;
		if(k == -1) {
			results[j] = 'N';
		}
		else {
			results[j] = bases[k];
		}
		
		if(k == -1) {
			results[j+n] = getRandBaseQuality();
		}
		else {
			basePairIndx = refIndx*N+k;
			results[j+n] = getBaseQuality(basePairIndx, binIndx);
		}
		
	}
	results[2*n] = '\0';
	return results;
}

void Profile::predict(char* refSeq, char* results, int num, int isRead1) {
	if(refSeq == NULL || refSeq[0] == '\0') {
		results[0] = '\0';
		return;
	}
	
	int i, j, k, quality;
	
	static string bases = config.getStringPara("bases");
	static int kmer = config.getIntPara("kmer");
	static int binCount = config.getIntPara("bins");
	static int N = bases.length();
	int n = strlen(refSeq);
	
	int refIndxs[n];
	char readSeq[n];
	char qualitySeq[n];
	
	char seq[kmer+n];
	for(i = 0; i < kmer-1; i++) {
		seq[i] = 'X';
	}
	for(; i < kmer+n-1; i++) {
		seq[i] = refSeq[i-kmer+1];
	}
	
	for(i = 0; i < n; i++) {
		refIndxs[i] = getIndexOfBase(refSeq[i]);
	}
	
	int basePairIndx, binIndx;
	
	int m = 0;
	for(i = 0; i < num; i++) {
		for(j = 0; j < n; j++) {
			binIndx = j*binCount/n;
			
			char c = seq[j+kmer];
			seq[j+kmer] = '\0';
			if(isRead1) {
				k = getSubBaseIndx1(&seq[j], binIndx);
			}
			else {
				k = getSubBaseIndx2(&seq[j], binIndx);
			}
			seq[j+kmer] = c;
			if(k == -1) {
				readSeq[j] = 'N';
			}
			else {
				readSeq[j] = bases[k];
			}
			
			if(k == -1) {
				qualitySeq[j] = minBaseQuality;
			}
			else {
				basePairIndx = refIndxs[j]*N+k;
				qualitySeq[j] = getBaseQuality(basePairIndx, binIndx);
			}
			
		}
		strncpy(&results[m], readSeq, n);
		strncpy(&results[m+n], qualitySeq, n);
		m += 2*n;
	}
	results[m] = '\0';
}
