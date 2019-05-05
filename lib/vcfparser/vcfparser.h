// ***************************************************************************
// vcfparser.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _VCFPARSER_H
#define _VCFPARSER_H

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

class SNV {
	public:
		SNV() {position = -1;}
		SNV(long position, char ref, char alt, varType type) :
			position(position), ref(ref), alt(alt), type(type) {}
		long getPosition() {return position;}
		char getRef() {return ref;}
        char getAlt() {return alt;}
        varType getType() {return type;}
	private:
		long position;
		char ref;
		char alt;
		varType type;
};

class Insert {
	public:
		Insert() {position = -1;}
		Insert(long position, string sequence, varType type) :
			position(position), sequence(sequence), type(type) {}
		long getPosition() {return position;}
		string getSequence() {return sequence;}
		int getLength() {return sequence.length();}
		varType getType() {return type;}
	private:
		long position;
		string sequence;
		varType type;
};

class Deletion {
	public:
		Deletion() {position = -1;}
		Deletion(long position, int length, varType type) :
			position(position), length(length), type(type) {}
		long getPosition() {return position;}
		int getLength() {return length;}
		varType getType() {return type;}
	private:
		long position;
		int length;
		varType type;
};

class VcfParser {
	public:
		VcfParser() {vcfFile = "";}
		VcfParser(string vcfFile) : vcfFile(vcfFile) {}
		string aberOfChr(string chrName);
		void parse();
		void setVCF(string vcfFile) {this->vcfFile = vcfFile;}
		map<string, vector<SNV> >& getSNVs() {return snvs;}
		map<string, vector<Insert> >& getInserts() {return inserts;}
		map<string, vector<Deletion> >& getDels() {return deletions;}
		vector<SNV>& getSNVs(string chr) {return snvs[chr];}
		vector<Insert>& getInserts(string chr) {return inserts[chr];}
		vector<Deletion>& getDels(string chr) {return deletions[chr];}
		long getNumOfSNVs();
		long getNumOfInserts();
		long getNumOfDels();
	private:
		string vcfFile;
		map<string, vector<SNV> > snvs;
		map<string, vector<Insert> > inserts;
		map<string, vector<Deletion> > deletions;
};


#endif
