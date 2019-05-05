// ***************************************************************************
// Segment.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "snp.h"
#include "MyDefine.h"
#include "Segment.h"

unsigned int Segment::segMaxSize = 1000000;
unsigned int Segment::fragSize = 1000;


Segment::Segment(int segIndx, string segChr, long segStartPos, long segEndPos) {
	assert(segIndx >= 0);
	assert(!segChr.empty());
	assert(segStartPos >= 0);
	assert(segEndPos >= segStartPos);
	this->segIndx = segIndx;
	this->segChr = segChr;
	this->segStartPos = segStartPos;
	this->segEndPos = segEndPos;
	this->CN = 2;
	this->mCN = 1;
	this->segSequences = NULL;
	
	readCount = 0;
	gcFactor = -1;

	initTargets();
	
	//threadFunc = &yieldReads;
	
	//yieldSegmentSequence();
}

Segment::Segment(int segIndx, string segChr, long segStartPos, long segEndPos, int CN, int mCN) {
	assert(segIndx >= 0);
	assert(!segChr.empty());
	assert(segStartPos >= 0);
	assert(segEndPos >= segStartPos);
	assert(mCN >= 0);
	assert(CN >= mCN);
	this->segIndx = segIndx;
	this->segChr = segChr;
	this->segStartPos = segStartPos;
	this->segEndPos = segEndPos;
	this->CN = CN;
	this->mCN = mCN;
	this->segSequences = NULL;
	
	readCount = 0;
	gcFactor = -1;

	initTargets();

	//threadFunc = &yieldReads;
	
	//yieldSegmentSequence();
}

void Segment::initTargets() {
	int i, j, k;
	if(!genome.getOutTargets().empty()) {
		vector<Target>& targetsOfChr = genome.getOutTargets(segChr);
		for(i = 0; i < targetsOfChr.size(); i++) {
			long spos = targetsOfChr[i].spos;
			long epos = targetsOfChr[i].epos;
			if((spos >= segStartPos && spos <= segEndPos) || (epos >= segStartPos && epos <= segEndPos) || (spos < segStartPos && epos > segEndPos)) {
				targetIndxs.push_back(i);
			}
		}
	}
}

char* Segment::getComplementSeq(char* sequence) {
	if(sequence == NULL) {
		return NULL;
	}
	int i, n = strlen(sequence);
	for(i = 0; i < n; i++) {
		char c;
		switch(sequence[i]){
			case 'A': c = 'T'; break;
			case 'T': c = 'A'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			case 'a': c = 't'; break;
			case 't': c = 'a'; break;
			case 'c': c = 'g'; break;
			case 'g': c = 'c'; break;
			case 'N': c = 'N'; break;
			default : c = 'N';
		}
		sequence[i] = c;
	}
	return sequence;
}

char* Segment::getSegSequence(int i) {
	assert(i >= 0 && i < config.getIntPara("ploidy"));
	return segSequences[i];
}

void Segment::clearSegSequences() {
	if(segSequences == NULL) {
		return;
	}
	int ploidy = config.getIntPara("ploidy");
	for(int i = 0; i < ploidy; i++) {
		if(segSequences[i] != NULL) {
			delete[] segSequences[i];
		}
	}
	delete[] segSequences;
	segSequences = NULL;
}

void Segment::generateSegSequences() {
	int ploidy = config.getIntPara("ploidy");
	vector<SNP>& snpsOfChr = genome.getSimuSNPs(segChr);
	vector<SNV>& snvsOfChr = genome.getSimuSNVs(genome.getCurrentPopu(), segChr);
	vector<Insert>& insertsOfChr = genome.getSimuInserts(genome.getCurrentPopu(), segChr);
	vector<Deletion>& delsOfChr = genome.getSimuDels(genome.getCurrentPopu(), segChr);
	
	if(CN == 0) {
		segSequences = NULL;
		return;
	}

	//char* refSeq = genome.getSubRefSequence(segChr, segStartPos-1, getRefSize());
	char* refSeq = genome.getSubSequence(segChr, segStartPos-1, getRefSize());
	if(refSeq == NULL) {
		segSequences = NULL;
		return;
	}

	transform(refSeq, refSeq+strlen(refSeq), refSeq, (int (*)(int))toupper);
	
	int i, j, k;
	unsigned int refSize = strlen(refSeq);
	vector<int>::iterator it;
	
	if(mIndx.empty()) {
		if(CN < ploidy) {
			for(i = 0; i < CN; i++) {
				while(1) {
					j = randomInteger(0, ploidy);
					it = find(seqReps.begin(), seqReps.end(), j);
					if(it == seqReps.end()) {
						seqReps.push_back(j);
						break;
					}
				}
			}
			for(i = 0; i < mCN; i++) {
				mIndx.push_back(seqReps[i]);
			}
		}
		else {
			for(i = 0; i < ploidy; i++) {
				seqReps.push_back(1);
			}
			int n = CN-ploidy;
			k = randomInteger(0, ploidy);
			for(i = n; i >= 0; i--) {
				if(seqReps[k]+i == mCN) {
					seqReps[k] += i;
					mIndx.push_back(k);
					break;
				}
				else if(seqReps[k]+i == CN-mCN) {
					seqReps[k] += i;
					for(j = 0; j < ploidy; j++) {
						if(j != k) {
							mIndx.push_back(j);
						}
					}
					break;
				}
			}
			if(i >= 0) {
				n -= i;
				while(n > 0) {
					j = randomInteger(0, ploidy);
					if(j != k) {
						seqReps[j]++;
						n--;
					}
				}
			}
			else {
				while(n > 0) {
					j = randomInteger(0, ploidy);
					seqReps[j]++;
					n--;
				}
				for(i = 0; i < ploidy; i++) {
					mIndx.push_back(i);
				}
			}
		}
	}
	
	vector<string> segSeqs;
	if(CN < ploidy) {
		for(i = 0; i < ploidy; i++) {
			it = find(seqReps.begin(), seqReps.end(), i);
			if(it != seqReps.end()) {
				segSeqs.push_back(refSeq);
			}
			else {
				segSeqs.push_back("");
			}
		}
	}
	else {
		for(i = 0; i < ploidy; i++) {
			string tmp = "";
			for(j = 0; j < seqReps[i]; j++) {
				tmp += refSeq;
			}
			segSeqs.push_back(tmp);
		}
	}

	delete[] refSeq;
	
	//simu SNP
	k = 0;
	for(i = 0; i < snpsOfChr.size(); i++) {
		SNP snp = snpsOfChr[i];
		long pos = snp.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			if(k == 0) {
				for(j = 0; j < mIndx.size(); j++) {
					string& segSeq = segSeqs[mIndx[j]];
					unsigned int segLen = segSeq.length();
					for(int t = 0; t < segLen/refSize; t++) {
						segSeq[sindx+t*refSize] = snp.getNucleotide();
					}
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if(it != mIndx.end()) {
						continue;
					}
					string& segSeq = segSeqs[j];
					unsigned int segLen = segSeq.length();
					for(int t = 0; t < segLen/refSize; t++) {
						segSeq[sindx+t*refSize] = snp.getNucleotide();
					}
				}
			}
			k = (k+1)%2;
		}
	}
	
	//simu SNV
	k = 0;
	for(i = 0; i < snvsOfChr.size(); i++) {
		SNV snv = snvsOfChr[i];
		long pos = snv.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			varType type = snv.getType();
			char alt = snv.getAlt();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					string& segSeq = segSeqs[j];
					unsigned int segLen = segSeq.length();
					for(int t = 0; t < segLen/refSize; t++) {
						segSeq[sindx+t*refSize] = alt;
					}
				}
			}
			else {
				if(k == 0) {
					for(j = 0; j < mIndx.size(); j++) {
						string& segSeq = segSeqs[mIndx[j]];
						unsigned int segLen = segSeq.length();
						for(int t = 0; t < segLen/refSize; t++) {
							segSeq[sindx+t*refSize] = alt;
						}
					}
				}
				else {
					for(j = 0; j < ploidy; j++) {
						it = find(mIndx.begin(), mIndx.end(), j);
						if(it != mIndx.end()) {
							continue;
						}
						string& segSeq = segSeqs[j];
						unsigned int segLen = segSeq.length();
						for(int t = 0; t < segLen/refSize; t++) {
							segSeq[sindx+t*refSize] = alt;
						}
					}
				}
				k = (k+1)%2;
			}
		}
	}
	
	//simu Insert
	map<int, map<int, int> > insertsPerhaploidy;
	map<int, int>::iterator m_it;
	int* insertLens = new int[ploidy];
	memset(insertLens, 0, ploidy*sizeof(int));
	k = 0;
	for(i = 0; i < insertsOfChr.size(); i++) {
		Insert insert = insertsOfChr[i];
		long pos = insert.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos+1-segStartPos; // insert before location pos+1
			varType type = insert.getType();
			string seq = insert.getSequence();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					string& segSeq = segSeqs[j];
					int n = segSeq.length()/(refSize+insertLens[j]);
					int len = seq.length();
					for(int t = 0; t < n; t++) {
						segSeq.insert(sindx+offset+t*(refSize+insertLens[j]+len), seq);
					}
					insertLens[j] += seq.length();
					insertedSeq.insert(make_pair(sindx, seq.length()));
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
						continue;
					}
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					string& segSeq = segSeqs[j];
					int n = segSeq.length()/(refSize+insertLens[j]);
					int len = seq.length();
					for(int t = 0; t < n; t++) {
						segSeq.insert(sindx+offset+t*(refSize+insertLens[j]+len), seq);
					}
					insertLens[j] += seq.length();
					insertedSeq.insert(make_pair(sindx, seq.length()));
				}
				k = (k+1)%2;
			}
		}
	}
	
	//simu Deletion
	map<int, map<int, int> > delsPerhaploidy;
	int* delLens = new int[ploidy];
	memset(delLens, 0, ploidy*sizeof(int));
	k = 0;
	for(i = 0; i < delsOfChr.size(); i++) {
		Deletion del = delsOfChr[i];
		long pos = del.getPosition();
		if(pos >= segStartPos && pos <= segEndPos) {
			int sindx = pos-segStartPos;
			varType type = del.getType();
			int delLen = del.getLength();
			if(type == homo) {
				for(j = 0; j < ploidy; j++) {
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					map<int, int>& delSeq = delsPerhaploidy[j];
					for(m_it = delSeq.begin(); m_it!= delSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset -= (*m_it).second;
						}
					}
					if(sindx+offset < 0) {
						continue;
					}
					string& segSeq = segSeqs[j];
					int n = segSeq.length()/(refSize+insertLens[j]-delLens[j]);
					for(int t = 0; t < n; t++) {
						segSeq.erase(sindx+offset+t*(refSize+insertLens[j]-delLens[j]-delLen), delLen);
					}
					delLens[j] += delLen;
					delSeq.insert(make_pair(sindx, delLen));
				}
			}
			else {
				for(j = 0; j < ploidy; j++) {
					it = find(mIndx.begin(), mIndx.end(), j);
					if((k == 0 && it == mIndx.end()) || (k == 1 && it != mIndx.end())) {
						continue;
					}
					int offset = 0;
					map<int, int>& insertedSeq = insertsPerhaploidy[j];
					for(m_it = insertedSeq.begin(); m_it!= insertedSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset += (*m_it).second;
						}
					}
					map<int, int>& delSeq = delsPerhaploidy[j];
					for(m_it = delSeq.begin(); m_it!= delSeq.end(); m_it++) {
						if((*m_it).first <= sindx) {
							offset -= (*m_it).second;
						}
					}
					if(sindx+offset < 0) {
						continue;
					}
					string& segSeq = segSeqs[j];
					int n = segSeq.length()/(refSize+insertLens[j]-delLens[j]);
					for(int t = 0; t < n; t++) {
						segSeq.erase(sindx+offset+t*(refSize+insertLens[j]-delLens[j]-delLen), delLen);
					}
					delLens[j] += delLen;
					delSeq.insert(make_pair(sindx, delLen));
				}
				k = (k+1)%2;
			}
		}
	}
	delete[] insertLens;
	delete[] delLens;
	
	segSequences = new char*[ploidy];
	for(i = 0; i < ploidy; i++) {
		if(segSeqs[i].length() == 0) {
			segSequences[i] = NULL;
		}
		else {
			segSequences[i] = new char[segSeqs[i].length()+1];
			strcpy(segSequences[i], segSeqs[i].c_str());
			transform(segSequences[i], segSequences[i]+strlen(segSequences[i]), segSequences[i], (int (*)(int))toupper);
		}
	}
	
}

void Segment::setReadCount(long readCount) {
	int i;
	double totalWL = getWeightedLength()+2.2204e-16;
	long sum = 0;
	fragRCs.clear();
	for(i = 0; i < fragWeights.size(); i++) {
		long rc = fragWeights[i]*readCount/totalWL;
		fragRCs.push_back(rc);
		sum += rc;
	}
	if(sum < readCount) {
		fragRCs[0] += readCount-sum;
	}
	this->readCount = readCount;
}

long Segment::getReadCount(int indx) {
	return fragRCs[indx];
}

long Segment::getReadCount() {
	return readCount;
}

/*
double Segment::getWeightedLength() {
	int i, j, k;
	double weightLen = 0;
	if(fragWeights.empty()) {
		unsigned int fragSize = Segment::fragSize;
		unsigned int segSize = getSeqSize()/CN;
		if(genome.getOutTargets().empty()) {
			k = segSize/fragSize;
			for(i = 0; i < k; i++) {
				long spos = i*fragSize;
				char* s = genome.getSubSequence(segChr, segStartPos+spos-1, fragSize);
				int gc = calculateGCPercent(s);
				double weight = profile.getGCFactor(gc)/fragSize;
				delete[] s;
				fragStartPos.push_back(spos);
				fragEndPos.push_back(spos+fragSize-1);
				fragWeights.push_back(weight);
				weightLen += weight;
			}
			if(k*fragSize < segSize) {
				long spos = k*fragSize;
				char* s = genome.getSubSequence(segChr, segStartPos+spos-1, segSize-k*fragSize);
				int gc = calculateGCPercent(s);
				double weight = profile.getGCFactor(gc)*(segSize-k*fragSize)/(fragSize*fragSize);
				delete[] s;
				fragStartPos.push_back(spos);
				fragEndPos.push_back(segSize-1);
				fragWeights.push_back(weight);
				weightLen += weight;
			}
		}
		else if(!targetIndxs.empty()) {
			vector<Target>& targetsOfChr = genome.getOutTargets(segChr);
			
			for(i = 0; i < targetIndxs.size(); i++) {
				j = targetIndxs[i];
				long spos = max(targetsOfChr[j].spos, segStartPos) - segStartPos;
				long epos = min(targetsOfChr[j].epos, segEndPos) - segStartPos;
				char* s = genome.getSubSequence(segChr, segStartPos+spos-1, epos-spos+1);
				int gc = calculateGCPercent(s);
				double weight = profile.getGCFactor(gc)*(epos-spos+1)/(fragSize*fragSize);
				delete[] s;
				fragStartPos.push_back(spos);
				fragEndPos.push_back(epos);
				fragWeights.push_back(weight);
				weightLen += weight;
			}
		}
		else {
			fragStartPos.push_back(0);
			fragEndPos.push_back(segSize-1);
			fragWeights.push_back(0);
		}
	}
	else {
		for(i = 0; i < fragWeights.size(); i++) {
			weightLen += fragWeights[i];
		}
	}
	return weightLen*CN/2;
}
*/

double Segment::getWeightedLength() {
	int i, j, k, n, m;
	int flag = 0;
	double weightLen = 0;
	int ploidy = config.getIntPara("ploidy");
	if(fragWeights.empty()) {
		if(segSequences == NULL) {
			flag = 1;
			generateSegSequences();
		}
		unsigned int fragSize = Segment::fragSize;
		unsigned int segSize = getSeqSize()/CN;
		if(genome.getOutTargets().empty()) {
			for(i = 0; i < ploidy; i++) {
				if(segSequences[i] == NULL) {
					continue;
				}
				char* p = segSequences[i];
				k = strlen(p)/fragSize;
				for(j = 0; j < k; j++) {
					long spos = j*fragSize;
					long epos = (j+1)*fragSize-1;
					char c = p[epos+1];
					p[epos+1] = '\0';
					int gc = calculateGCPercent(p+spos);
					p[epos+1] = c;
					double weight = profile.getGCFactor(gc)/fragSize;
					fragStartPos.push_back(spos);
					fragEndPos.push_back(epos);
					fragWeights.push_back(weight);
					hapIndxs.push_back(i);
					weightLen += weight;
				}
				if(k*fragSize < strlen(p)) {
					long spos = k*fragSize;
					int gc = calculateGCPercent(p+spos);
					double weight = profile.getGCFactor(gc)*(strlen(p)-spos)/(fragSize*fragSize);
					fragStartPos.push_back(spos);
					fragEndPos.push_back(strlen(p)-1);
					fragWeights.push_back(weight);
					hapIndxs.push_back(i);
					weightLen += weight;
				}
			}
		}
		else if(!targetIndxs.empty()) {
			vector<Target>& targetsOfChr = genome.getOutTargets(segChr);
			for(i = 0; i < ploidy; i++) {
				if(segSequences[i] == NULL) {
					continue;
				}
				char* p = segSequences[i];
				n = (seqReps.size() < ploidy)? 1:seqReps[i];
				long refLen = strlen(p)/n;
				for(k = 0; k < n; k++) {
					for(j = 0; j < targetIndxs.size(); j++) {
						m = targetIndxs[j];
						long spos = max(targetsOfChr[m].spos, segStartPos) - segStartPos;
						long epos = min(targetsOfChr[m].epos, segStartPos+refLen-1) - segStartPos;
						long spos_k = spos+k*strlen(p)/n;
						long epos_k = epos+k*strlen(p)/n;
						char c = p[epos_k+1];
						p[epos_k+1] = '\0';
						int gc = calculateGCPercent(p+spos_k);
						p[epos_k+1] = c;
						double weight = profile.getGCFactor(gc)*(epos_k-spos_k+1)/(fragSize*fragSize);
						fragStartPos.push_back(spos_k);
						fragEndPos.push_back(epos_k);
						fragWeights.push_back(weight);
						hapIndxs.push_back(i);
						weightLen += weight;
					}
				}
			}
		}
		else {
			fragStartPos.push_back(0);
			fragEndPos.push_back(0);
			fragWeights.push_back(0);
			hapIndxs.push_back(0);
		}
		if(flag == 1) {
			clearSegSequences();
		}
	}
	else {
		for(i = 0; i < fragWeights.size(); i++) {
			weightLen += fragWeights[i];
		}
	}
	return weightLen;
}

unsigned int Segment::getSeqSize() {
	if(this->segSequences == NULL) {
		return CN*getRefSize();
	}
	int ploidy = config.getIntPara("ploidy");
	unsigned int seqSize = 0;
	for(int i = 0; i < ploidy; i++) {
		if(segSequences[i] != NULL) {
			seqSize += strlen(segSequences[i]);
		}	
	}
	return seqSize;
}

unsigned int Segment::getTargetsSize() {
	if(targetIndxs.empty()) {
		return getRefSize();
	}
	vector<Target>& targetsOfChr = genome.getOutTargets(segChr);
	long sum = 0;
	int i, j;
	for(i = 0; i < targetIndxs.size(); i++) {
		j = targetIndxs[i];
		long spos = max(targetsOfChr[j].spos, segStartPos);
		long epos = min(targetsOfChr[j].epos, segEndPos);
		sum += epos-spos+1;
	}
	return sum;
}

void* Segment::yieldReads(const void* args) {
	Segment* seg = (Segment*) args;
	if(seg->getSegSequences() == NULL || seg->getReadCount() == 0) {
		return NULL;
	}

	int segIndx = seg->getSegIndx();
	string segChr = seg->getSegChr();
	long segStartPos = seg->getSegStartPos();
	long segEndPos = seg->getSegEndPos();
	int CN = seg->getCN();
	char** segSequences = seg->getSegSequences();
	vector<long>& fragStartPos = seg->getFragStartPos();
	vector<long>& fragEndPos = seg->getFragEndPos();
	vector<int>& hapIndxs = seg->getHapIndxs();

	bool paired = config.isPairedEnd();
	int seqLen, readLength = config.getIntPara("readLength");
	int ploidy = config.getIntPara("ploidy");
	
	//cerr << "Number of reads to sample on the segment: " << seg->getReadCount() << endl;
	
	unsigned long bufferSize = 50000000;
	char buf[10*readLength], buf1[10*readLength], buf2[10*readLength];
	char *results;
	char *outBuffer, *outBuffer1, *outBuffer2;
	unsigned long outIndx = 0, outIndx1 = 0, outIndx2 = 0;
	string outputDir = config.getStringPara("outputDir");
	if(paired) {
		outBuffer1 = new char[bufferSize];
		outBuffer2 = new char[bufferSize];
	}
	else {
		outBuffer = new char[bufferSize];
	}
	
	int i, j, k;
	
	unsigned int refSize = seg->getRefSize();
	unsigned int seqSize, segsize;
	seqSize = seg->getSeqSize();
	segsize = seqSize/CN;
	
	
	vector<int> segSeqIndxs;
	vector<int> cnStartCount;
	int count = 0;
	for(i = 0; i < ploidy; i++) {
		cnStartCount.push_back(count);
		if(segSequences[i] == NULL) {
			continue;
		}
		k = strlen(segSequences[i])/segsize;
		count += k;
		for(j = 0; j < k; j++) {
			segSeqIndxs.push_back(i);
		}
	}
	
	int fragCount = 0;
	int lastSegIndx, fragIndx;
	
	for(i = 0; i < fragStartPos.size(); i++) {
		long spos = fragStartPos[i];
		long epos = fragEndPos[i];
		int hapIndx = hapIndxs[i];
		long fragSize = epos-spos+1;
		int failCount = 0;
		int n = seg->getReadCount(i);
		while(n > 0) {
			long pos = threadPool->randomInteger(spos, epos+1);
			char* fragSeq;
			if(!paired) {
				fragSeq = getFragSequence(seg, hapIndx, pos, fragSize);
			}
			else {
				int insertSize = profile.yieldInsertSize();
				fragSeq = getFragSequence(seg, hapIndx, pos, insertSize);			
			}
			
			if(fragSeq == NULL || strlen(fragSeq) < readLength) {
				if(fragSeq != NULL) {
					delete[] fragSeq;
				}
				failCount++;
				if(failCount > 1000) {
					break;
				}
				continue;
			}			
			fragCount++;
			
			if(!paired) {
				k = threadPool->randomInteger(0, 2);
				if(k == 0) {
					char c = fragSeq[readLength];
					fragSeq[readLength] = '\0';
					results = profile.predict(fragSeq, 1);
					fragSeq[readLength] = c;
				}
				else {
					char* seq = &fragSeq[strlen(fragSeq)-readLength];
					seq = getComplementSeq(seq);
					reverse(seq, seq+readLength);
					results = profile.predict(seq, 1);
				}
				k = 0;
				sprintf(&buf[k], "@%s#%s#%ld#%d\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf);
				seqLen = strlen(results)/2;
				strncpy(&buf[k], &results[0], seqLen);
				strcpy(&buf[k+seqLen], "\n+\n");
				strncpy(&buf[k+seqLen+3], &results[seqLen], seqLen);
				buf[k+2*seqLen+3] = '\n';
				buf[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf)+outIndx < bufferSize) {
					strcpy(&outBuffer[outIndx], buf);
					outIndx += strlen(buf);
				}
				else {
					swp->write(outBuffer);
					strcpy(outBuffer, buf);
					outIndx = strlen(buf);
				}
				
				n--;
			}
			else {
				char c = fragSeq[readLength];
				fragSeq[readLength] = '\0';
				results = profile.predict(fragSeq, 1);
				fragSeq[readLength] = c;
				
				k = 0;
				sprintf(&buf1[k], "@%s#%s#%ld#%d/1\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf1);
				seqLen = strlen(results)/2;
				strncpy(&buf1[k], &results[0], seqLen);
				strcpy(&buf1[k+seqLen], "\n+\n");
				strncpy(&buf1[k+seqLen+3], &results[seqLen], seqLen);
				buf1[k+2*seqLen+3] = '\n';
				buf1[k+2*seqLen+4] = '\0';
				delete[] results;
				
				char* seq = &fragSeq[strlen(fragSeq)-readLength];
				seq = getComplementSeq(seq);
				reverse(seq, seq+readLength);
				results = profile.predict(seq, 0);
				k = 0;
				sprintf(&buf2[k], "@%s#%s#%ld#%d/2\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf2);
				seqLen = strlen(results)/2;
				strncpy(&buf2[k], &results[0], seqLen);
				strcpy(&buf2[k+seqLen], "\n+\n");
				strncpy(&buf2[k+seqLen+3], &results[seqLen], seqLen);
				buf2[k+2*seqLen+3] = '\n';
				buf2[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf1)+outIndx1 < bufferSize && strlen(buf2)+outIndx2 < bufferSize) {
					strcpy(&outBuffer1[outIndx1], buf1);
					strcpy(&outBuffer2[outIndx2], buf2);
					outIndx1 += strlen(buf1);
					outIndx2 += strlen(buf2);
				}
				else {
					swp->write(outBuffer1, outBuffer2);
					strcpy(outBuffer1, buf1);
					strcpy(outBuffer2, buf2);
					outIndx1 = strlen(buf1);
					outIndx2 = strlen(buf2);
				}
				
				n -= 2;
			}
			
			delete[] fragSeq;
		}
	}
	
	
	if(paired) {
		if(strlen(outBuffer1) > 0) {
			swp->write(outBuffer1, outBuffer2);
		}
		delete[] outBuffer1;
		delete[] outBuffer2;
	}
	else {
		if(strlen(outBuffer) > 0) {
			swp->write(outBuffer);
		}
		delete[] outBuffer;
	}
	
	return NULL;
}

/*
void* Segment::yieldReads(const void* args) {
	Segment* seg = (Segment*) args;
	if(seg->getSegSequences() == NULL || seg->getReadCount() == 0) {
		return NULL;
	}

	int segIndx = seg->getSegIndx();
	string segChr = seg->getSegChr();
	long segStartPos = seg->getSegStartPos();
	long segEndPos = seg->getSegEndPos();
	int CN = seg->getCN();
	char** segSequences = seg->getSegSequences();
	vector<long>& fragStartPos = seg->getFragStartPos();
	vector<long>& fragEndPos = seg->getFragEndPos();

	bool paired = config.isPairedEnd();
	int seqLen, readLength = config.getIntPara("readLength");
	int ploidy = config.getIntPara("ploidy");
	
	//cerr << "Number of reads to sample on the segment: " << seg->getReadCount() << endl;
	
	unsigned long bufferSize = 50000000;
	char buf[10*readLength], buf1[10*readLength], buf2[10*readLength];
	char *results;
	char *outBuffer, *outBuffer1, *outBuffer2;
	unsigned long outIndx = 0, outIndx1 = 0, outIndx2 = 0;
	string outputDir = config.getStringPara("outputDir");
	if(paired) {
		outBuffer1 = new char[bufferSize];
		outBuffer2 = new char[bufferSize];
	}
	else {
		outBuffer = new char[bufferSize];
	}
	
	int i, j, k;
	
	unsigned int refSize = seg->getRefSize();
	unsigned int seqSize, segsize;
	seqSize = seg->getSeqSize();
	segsize = seqSize/CN;
	
	
	vector<int> segSeqIndxs;
	vector<int> cnStartCount;
	int count = 0;
	for(i = 0; i < ploidy; i++) {
		cnStartCount.push_back(count);
		if(segSequences[i] == NULL) {
			continue;
		}
		k = strlen(segSequences[i])/segsize;
		count += k;
		for(j = 0; j < k; j++) {
			segSeqIndxs.push_back(i);
		}
	}
	
	int fragCount = 0;
	int lastSegIndx, fragIndx;
	
	for(i = 0; i < fragStartPos.size(); i++) {
		long spos = fragStartPos[i];
		long epos = fragEndPos[i];
		long fragSize = epos-spos+1;
		int failCount = 0;
		int n = seg->getReadCount(i);
		while(n > 0) {
			long pos = threadPool->randomInteger(0, fragSize*CN);
			k = pos/fragSize;
			pos = pos%fragSize+spos;
			j = segSeqIndxs[k];
			pos += (k-cnStartCount[j])*segsize;
			char* fragSeq;
			if(!paired) {
				fragSeq = getFragSequence(seg, j, pos, fragSize);
			}
			else {
				int insertSize = profile.yieldInsertSize();
				fragSeq = getFragSequence(seg, j, pos, insertSize);			
			}
			
			if(fragSeq == NULL || strlen(fragSeq) < readLength) {
				if(fragSeq != NULL) {
					delete[] fragSeq;
				}
				failCount++;
				if(failCount > 1000) {
					break;
				}
				continue;
			}			
			fragCount++;
			
			if(!paired) {
				k = threadPool->randomInteger(0, 2);
				if(k == 0) {
					char c = fragSeq[readLength];
					fragSeq[readLength] = '\0';
					results = profile.predict(fragSeq, 1);
					fragSeq[readLength] = c;
				}
				else {
					char* seq = &fragSeq[strlen(fragSeq)-readLength];
					seq = getComplementSeq(seq);
					reverse(seq, seq+readLength);
					results = profile.predict(seq, 1);
				}
				k = 0;
				sprintf(&buf[k], "@%s#%s#%ld#%d\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf);
				seqLen = strlen(results)/2;
				strncpy(&buf[k], &results[0], seqLen);
				strcpy(&buf[k+seqLen], "\n+\n");
				strncpy(&buf[k+seqLen+3], &results[seqLen], seqLen);
				buf[k+2*seqLen+3] = '\n';
				buf[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf)+outIndx < bufferSize) {
					strcpy(&outBuffer[outIndx], buf);
					outIndx += strlen(buf);
				}
				else {
					swp->write(outBuffer);
					strcpy(outBuffer, buf);
					outIndx = strlen(buf);
				}
				
				n--;
			}
			else {
				char c = fragSeq[readLength];
				fragSeq[readLength] = '\0';
				results = profile.predict(fragSeq, 1);
				fragSeq[readLength] = c;
				
				k = 0;
				sprintf(&buf1[k], "@%s#%s#%ld#%d/1\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf1);
				seqLen = strlen(results)/2;
				strncpy(&buf1[k], &results[0], seqLen);
				strcpy(&buf1[k+seqLen], "\n+\n");
				strncpy(&buf1[k+seqLen+3], &results[seqLen], seqLen);
				buf1[k+2*seqLen+3] = '\n';
				buf1[k+2*seqLen+4] = '\0';
				delete[] results;
				
				char* seq = &fragSeq[strlen(fragSeq)-readLength];
				seq = getComplementSeq(seq);
				reverse(seq, seq+readLength);
				results = profile.predict(seq, 0);
				k = 0;
				sprintf(&buf2[k], "@%s#%s#%ld#%d/2\n", genome.getCurrentPopu().c_str(), segChr.c_str(), pos%segsize, fragCount);
				k = strlen(buf2);
				seqLen = strlen(results)/2;
				strncpy(&buf2[k], &results[0], seqLen);
				strcpy(&buf2[k+seqLen], "\n+\n");
				strncpy(&buf2[k+seqLen+3], &results[seqLen], seqLen);
				buf2[k+2*seqLen+3] = '\n';
				buf2[k+2*seqLen+4] = '\0';
				delete[] results;
				
				if(strlen(buf1)+outIndx1 < bufferSize && strlen(buf2)+outIndx2 < bufferSize) {
					strcpy(&outBuffer1[outIndx1], buf1);
					strcpy(&outBuffer2[outIndx2], buf2);
					outIndx1 += strlen(buf1);
					outIndx2 += strlen(buf2);
				}
				else {
					swp->write(outBuffer1, outBuffer2);
					strcpy(outBuffer1, buf1);
					strcpy(outBuffer2, buf2);
					outIndx1 = strlen(buf1);
					outIndx2 = strlen(buf2);
				}
				
				n -= 2;
			}
			
			delete[] fragSeq;
		}
	}
	
	
	if(paired) {
		if(strlen(outBuffer1) > 0) {
			swp->write(outBuffer1, outBuffer2);
		}
		delete[] outBuffer1;
		delete[] outBuffer2;
	}
	else {
		if(strlen(outBuffer) > 0) {
			swp->write(outBuffer);
		}
		delete[] outBuffer;
	}
	
	return NULL;
}
*/

char* Segment::getFragSequence(Segment* seg, int segSeqIndx, long pos, int fragSize) {
	char** segSequences = seg->getSegSequences();
	char* s;
	if(strlen(segSequences[segSeqIndx])-pos >= fragSize) {
		s = new char[fragSize+1];
		strncpy(s, segSequences[segSeqIndx]+pos, fragSize);
		s[fragSize] = '\0';
	}
	else {
		int i = strlen(segSequences[segSeqIndx])-pos;
		int k = fragSize-i;
		char* p = genome.produceFragment(genome.getCurrentPopu(), seg->getSegChr(), seg->getSegIndx()+1, segSeqIndx, k);
		if(p != NULL) {
			s = new char[i+strlen(p)+1];
			strncpy(s, segSequences[segSeqIndx]+pos, i);
			strcpy(s+i, p);
			//s[i+strlen(p)] = '\0';
			delete[] p;
		}
		else {
			s = new char[i+1];
			strncpy(s, segSequences[segSeqIndx]+pos, i);
			s[i] = '\0';
		}
	}
	return s;
}
