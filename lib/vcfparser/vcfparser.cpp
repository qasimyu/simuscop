// ***************************************************************************
// vcfparser.cpp (c) 2013 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstring>

#include "vcfparser.h"
#include "split.h"

string VcfParser::aberOfChr(string chrName) {
	size_t indx = chrName.find("chrom");
	if(indx == string::npos){
		indx = chrName.find("chr");
		if(indx != string::npos){
			chrName = chrName.substr(indx+3,chrName.size()-3);
		}
	}
	else{
		chrName = chrName.substr(indx+5,chrName.size()-5);
	}
	return chrName;
}

void VcfParser::parse() {
	if(vcfFile.empty()) {
		cerr << "Error: VCF file was not specified" << endl;
		return;
	}
	FILE *fp = fopen(vcfFile.c_str(),"r");
	if(!fp){
		cerr << "Error: cannot open VCF file " << vcfFile << endl;
		exit(-1);
	}
	varType type;
	
	float quality_th = 20;
	int depth_th = 10;
	char buf[20000];
	char * elems[10];
	int elemnum;
	long line_num = 0;
	int wrong_num = 0;
	while(fgets(buf, 20000, fp)){
		line_num++;
		if(buf[0] == '#') {
			continue;		
		}
		elemnum = split(buf, '\t', elems);
		if(elemnum < 10) {
			cerr << "Warning: malformed VCF file " << vcfFile <<
					", there should be at least 10 fields @line " << line_num << endl;
			cerr << buf << endl;
			wrong_num++;
			if(wrong_num > 10) {
				exit(1);
			}
			continue;
		}
		
		string info = elems[7];
		size_t indx = info.find("DP=");
		if(indx != string::npos) {
			size_t indx1 = info.find(";", indx);
			int depth = atoi(info.substr(indx+3, indx1-indx-3).c_str());
			//cerr << depth << endl;
			if(depth < depth_th) {
				continue;
			}
		}
		float quality = atof(elems[5]);
		if(quality < quality_th) {
			continue;
		}
		
		string chr = elems[0];
		chr = aberOfChr(chr);
		long pos = atol(elems[1]);
		char * fields[10];
		split(elems[9], ':', fields);
		if(strcmp(fields[0], "1/1") == 0) {
			type = het;
		}
		else {
			type = homo;
		}
		if(strlen(elems[3]) > 1) { //deletion
			Deletion del(pos+1, strlen(elems[3])-1, type);
			deletions[chr].push_back(del);
		}
		else if(strlen(elems[4]) > 1) { //insert
			Insert insert(pos, elems[4]+1, type);
			inserts[chr].push_back(insert);
		}
		else { //SNV
			SNV snv(pos, *elems[3], *elems[4], type);
			snvs[chr].push_back(snv);
		}
	}
	cerr << "total " << getNumOfSNVs() << " SNPs, " << getNumOfInserts() << " inserts and "
			<< getNumOfDels() << " deletions were loaded from file " << vcfFile << endl;
}

long VcfParser::getNumOfSNVs() {
	long snvNum = 0;
	map<string, vector<SNV> >::iterator it;
	for(it = snvs.begin(); it != snvs.end(); it++){
		snvNum += it->second.size();
	}
	return snvNum;
}

long VcfParser::getNumOfInserts() {
	long insertNum = 0;
	map<string, vector<Insert> >::iterator it;
	for(it = inserts.begin(); it != inserts.end(); it++){
		insertNum += it->second.size();
	}
	return insertNum;
}

long VcfParser::getNumOfDels() {
	long delNum = 0;
	map<string, vector<Deletion> >::iterator it;
	for(it = deletions.begin(); it != deletions.end(); it++){
		delNum += it->second.size();
	}
	return delNum;
}
