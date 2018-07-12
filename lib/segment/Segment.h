// ***************************************************************************
// Segment.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _SEGMENT_H
#define _SEGMENT_H

#include <cassert>
#include <vector>
#include <string>
#include <pthread.h>

#include "Matrix.h"

using namespace std;

class Segment {
	private:
		static unsigned int segMaxSize;
		
		int segIndx;
		string segChr;
		long segStartPos, segEndPos;
		int CN, mCN;
		
		//segment sequences
		char** segSequences;
		vector<int> mIndx;
		vector<int> seqReps;
		
		//target regions for exome sequencing
		vector<int> targetIndxs;
		
		long readCount;
		int replicates;
		double gcFactor;

		//map<long, double> fragWeights;
		//map<long, long> fragRCs;
		vector<long> fragStartPos;
		vector<long> fragEndPos;
		vector<long> fragRCs;
		vector<double> fragWeights;

		void initTargets();
		
		static char* getFragSequence(Segment* seg, int segSeqIndx, long pos, int fragSize);
		
	public:	
		static unsigned int fragSize;
	
		Segment(int segIndx, string segChr, long segStartPos, long segEndPos);
		Segment(int segIndx, string segChr, long segStartPos, long segEndPos, int CN, int mCN);
		~Segment() {}

		static char* getComplementSeq(char* sequence);

		static void setSegMaxSize(unsigned int size) {segMaxSize = size;}
		static unsigned int getSegMaxSize() {return segMaxSize;}
		static void setFragmentSize(unsigned int size) {fragSize = size;}
		static unsigned int getFragmentSize() {return fragSize;}
		
		int getSegIndx() {return segIndx;}
		string getSegChr() {return segChr;}
		long getSegStartPos() {return segStartPos;}
		long getSegEndPos() {return segEndPos;}
		int getCN() {return CN;}
		int getmCN() {return mCN;}
		
		vector<long>& getFragStartPos() {return fragStartPos;}
		vector<long>& getFragEndPos() {return fragEndPos;}
		vector<long>& getFragRCs() {return fragRCs;}

		unsigned int getSeqSize();
		unsigned int getRefSize() {return segEndPos-segStartPos+1;}
		unsigned int getTargetsSize();
		
		char* getSegSequence(int i);
		
		char** getSegSequences() {return segSequences;}

		void clearSegSequences();
		
		void setReadCount(long readCount);
		
		long getReadCount();
		long getReadCount(int indx);
		
		double getWeightedLength();
		
		void generateSegSequences();
		static void* yieldReads(void* args);
		
};


#endif

