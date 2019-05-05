// ***************************************************************************
// Genome.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _GENOME_H
#define _GENOME_H

#include <vector>
#include <map>
#include <string>

#include "Segment.h"
#include "snp.h"
#include "vcfparser.h"
#include "Fasta.h"

using namespace std;

class CNV {
	public:
		CNV() {spos = -1;}
		CNV(long spos, long epos, float CN, float mCN)
			: spos(spos), epos(epos), CN(CN), mCN(mCN) {}
		long getStartPos() {return spos;}
		long getEndPos() {return epos;}
		float getCopyNumber() {return CN;}
		float getMajorCopyNumber() {return mCN;}
		
		long spos;
		long epos;
		float CN;
		float mCN;
};

class Target {
	public:
		Target() {}
		Target(long spos, long epos)
			: spos(spos), epos(epos) {}
		long spos;
		long epos;
};

template <class dataType>
class PopuData : public map<string, vector<dataType> > {
	public:
		PopuData() {}
		~PopuData() {}
		//int getElementCount();
	private:
		int elemCount;
};
/*
template <class dataType>
int PopuData<dataType>::getElementCount() {
	map<string, vector<dataType> >::iterator it;
	int sum = 0;
	for(it = this->begin(); it != this->end(); it++) {
		sum += (it->second).size();
	}
	return sum;
}
*/

class Genome {
	private:
		vector<string> chromosomes;
		map<string, PopuData<CNV> > cnvs;
		map<string, PopuData<SNV> > snvs;
		map<string, PopuData<Insert> > inserts;
		map<string, PopuData<Deletion> > dels;
		
		map<string, PopuData<Segment> > segments;
		map<string, vector<Target> > inTargets;
		map<string, vector<Target> > outTargets;
		
		vector<vector<float> > mixProps;
		
		SNPOnChr sc;
		VcfParser vcfParser;
		FastaReference fr;
		
		
		//map<string, string> altSequence;
		string refSequence;
		string altSequence;
	
		string curPopu;
		string curChr;
		
		void loadAbers();
		void loadSNPs();
		void loadRefSeq();
		void loadTargets(map<string, vector<Target> >& targets, string targetFile);
		void loadAbundance();
		//void generateAltSequence();
		//void generateAltSequence(string chr);
		void generateChrSequence(string chr);
		void divideTargets(map<string, vector<Target> >& allTargets);
		
		void divideSegment(string popu, string chr, long segStartPos, long segEndPos, int CN, int mCN, int &segIndx);
		void calculateACNs(map<string, double>& ACNs);		
		
	public:
		Genome() {}
		~Genome() {}

		void loadData();
		void loadTrainData();

		vector<string>& getChroms() {return chromosomes;}
		
		vector<CNV>& getSimuCNVs(string popu, string chr);
		//vector<CNV>& getRealCNVs(string chr);
		vector<SNV>& getSimuSNVs(string popu, string chr);
		vector<SNV>& getRealSNVs(string chr) {return vcfParser.getSNVs(chr);}
		vector<SNP>& getSimuSNPs(string chr) {return sc[chr];}
		vector<Insert>& getSimuInserts(string popu, string chr);
		vector<Insert>& getRealInserts(string chr) {return vcfParser.getInserts(chr);}
		vector<Deletion>& getSimuDels(string popu, string chr);
		vector<Deletion>& getRealDels(string chr) {return vcfParser.getDels(chr);}
		vector<Segment>& getSegs(string popu, string chr);
		vector<Target>& getInTargets(string chr) {return inTargets[chr];}
		vector<Target>& getOutTargets(string chr) {return outTargets[chr];}

		map<string, vector<Target> >& getInTargets() {return inTargets;}
		map<string, vector<Target> >& getOutTargets() {return outTargets;}
		
		string getCurrentPopu() {return curPopu;}
		
		long getChromLen(string chr);
		long getGenomeLength();
		long getTargetLength();
		
		char* getSubSequence(string chr, int startPos, int length);
		char* getSubRefSequence(string chr, int startPos, int length);
		char* getSubAltSequence(string chr, int startPos, int length);
		
		int getSegReadCount(string popu, string chr, int segIndx);
		int getSegIndxWithOffset(string popu, string chr, int curSegIndx, int offset);
		void setReadCounts(string popu, long reads);
		
		void generateSegments();
		char* produceFragment(string popu, string chr, int startSegIndx, int segSeqIndx, int fragLen);
		void yieldReads();
};


#endif

