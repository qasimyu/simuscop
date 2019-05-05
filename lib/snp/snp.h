// ***************************************************************************
// snp.h (c) 2013 Zhenhua Yu <yzh163@mail.ustc.edu.cn>
// HI Lab, University of Science and Technology of China
// All rights reserved.

#ifndef _SNP_H
#define _SNP_H

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctype.h>

using namespace std;

#ifndef _VARTYPE_
#define _VARTYPE_
enum varType{het, homo, unknown};
#endif

class SNP {

	public:
		SNP();
		SNP(string name, long long position, string observed, char strand, char ref);
		SNP(string name, long long position, char ref, char nucleotide);
		virtual ~SNP();
		string getName();
		long long getPosition();
		char getStrand();
		string getObserved();
		char getRef();
        char getNucleotide();
        float getFrequency();
		char getComplement(char nucleotide);
        varType getType();
        void setFrequency(float value);
        void setType(varType type);
	private:
		string name;
		//string chromosome;
		long long position;
		char strand;
		string observed;
		char ref;
		char nucleotide;
		float frequency;
		varType type;
};

class SNPOnChr : public map<string, vector<SNP> > {
	public:
		SNPOnChr();
        virtual ~SNPOnChr();
		string aberOfChr(string chromosome);
		void readSNPs(string fname);
		void readSNPsFromVCF(string fname);
		vector<string> getChroms();
		long SNPNumber();
	private:
		FILE * SNPFile;
		long snp_num;
};


#endif
