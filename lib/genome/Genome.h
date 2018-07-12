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
#include "Fasta.h"

using namespace std;

class CNV {
	public:
		CNV() {}
		CNV(long spos, long epos, float CN, float mCN)
			: spos(spos), epos(epos), CN(CN), mCN(mCN) {}
		~CNV() {}
		long spos;
		long epos;
		float CN;
		float mCN;
};

class SNV {
	public:
		SNV() {}
		SNV(long pos, char ref, char alt, string type)
			: pos(pos), ref(ref), alt(alt), type(type) {}
		~SNV() {}
		long pos;
		char ref;
		char alt;
		string type;
};

class Insert {
	public:
		Insert() {}
		Insert(long pos, string seq) : pos(pos), seq(seq) {}
		~Insert() {}
		long pos;
		string seq;
};

class Del {
	public:
		Del() {}
		Del(long pos, unsigned int len) : pos(pos), len(len) {}
		~Del() {}
		long pos;
		unsigned int len;
};

class Target {
	public:
		Target() {}
		Target(long spos, long epos)
			: spos(spos), epos(epos) {}
		~Target() {}
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
		map<string, PopuData<Del> > dels;
		
		map<string, PopuData<Segment> > segments;
		map<string, vector<Target> > inTargets;
		map<string, vector<Target> > outTargets;
		
		vector<vector<float> > mixProps;
		
		SNPOnChr scReal, scSim;	
		FastaReference fr;
		
		
		//map<string, string> altSequence;
		string altSequence;
	
		string curPopu;
		string curChr;
		
		void loadAbers();
		void loadSNPs(string snpFile, int n);
		void loadRefSeq(string refFile);
		void loadTargets(map<string, vector<Target> >& targets, string targetFile);
		void loadAbundance();
		//void generateAltSequence();
		void generateAltSequence(string chr);
		void divideTargets(map<string, vector<Target> >& allTargets);
		
		void divideSegment(string popu, string chr, long segStartPos, long segEndPos, int CN, int mCN, int &segIndx);
		void calculateACNs(map<string, double>& ACNs);		
		
	public:
		Genome() {}
		~Genome() {}

		void loadData();
		void loadData(string vcfFile, string refFile, string targetFile);

		vector<string>& getChroms() {return chromosomes;}
		
		vector<CNV>& getCNVs(string popu, string chr);
		vector<SNV>& getSNVs(string popu, string chr);
		vector<SNP>& getRealSNPs(string chr) {return scReal[chr];}
		vector<SNP>& getSimuSNPs(string chr) {return scSim[chr];}
		vector<Insert>& getInserts(string popu, string chr);
		vector<Del>& getDels(string popu, string chr);
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
		char* getSubAltSequence(string chr, int startPos, int length);
		
		int getSegReadCount(string popu, string chr, int segIndx);
		int getSegIndxWithOffset(string popu, string chr, int curSegIndx, int offset);
		void setReadCounts(string popu, long reads);
		
		void generateSegments();
		char* produceFragment(string popu, string chr, int startSegIndx, int segSeqIndx, int fragLen);
		void yieldReads();
};


#endif

