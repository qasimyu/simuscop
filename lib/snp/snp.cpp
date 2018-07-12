// ***************************************************************************
// snp.cpp (c) 2013 zhenhua yu <yzh163@mail.ustc.edu.cn>
// HI Lab, University of Science and Technology of China
// All rights reserved.

#include <cstring>

#include "snp.h"
#include "split.h"


SNP::SNP(string name, long long position, string observed, char strand, char ref) {
	this->name = name;
	this->position = position;
	this->observed = observed;
	this->strand = strand;
	this->ref = ref;
	vector<string> strs = split(observed, '/');
	if(strand == '-'){
		ref = getComplement(ref);
	}
	if(strs[0].at(0) == ref) {
	nucleotide = strs[1].at(0);
	} 
	else{
	nucleotide = strs[0].at(0);
	}

	if(strand == '-'){
	nucleotide = getComplement(nucleotide);
	}

	frequency = 0;
	status = unknow;
}

SNP::SNP(string name, long long position, char ref, char nucleotide) {
	this->name = name;
	this->position = position;
	this->ref = ref;
	this->nucleotide = nucleotide;
	
	frequency = 0;
	status = unknow;
}

SNP::SNP(void) {
}

SNP::~SNP(void) {
}

string SNP::getName() {
	return name;
}

long long SNP::getPosition() {
	return position;
}

char SNP::getStrand() {
	return strand;
}

string SNP::getObserved() {
	return observed;
}

char SNP::getRef() {
	return ref;
}

char SNP::getNucleotide() {
	return nucleotide;
}

float SNP::getFrequency() {
	return frequency;
}

snp_status SNP::getStatus() {
	return status;
}

char SNP::getComplement(char nucleotide) {
	switch(nucleotide){
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		case 'N': return 'N';
		default : return 'N';
	}
}

void SNP::setFrequency(float value) {
	this->frequency = value;
}

void SNP::setStatus(snp_status status) {
	this->status = status;
}


SNPOnChr::SNPOnChr(void) {
	snp_num = 0;
	SNPFile = NULL;
}

SNPOnChr::~SNPOnChr(void) {
	if(SNPFile != NULL) {
		fclose(SNPFile);
	}
}

vector<string> SNPOnChr::getChroms() {
	vector<string> chroms;
	map<string, vector<SNP> >::iterator it;
	for(it = this->begin(); it != this->end(); it++){
		chroms.push_back(it->first);
	}
	return chroms;
}

long SNPOnChr::SNPNumber() {
	return snp_num;
}

string SNPOnChr::aberOfChr(string chromosome) {
	size_t indx = chromosome.find("chrom");
	if(indx == string::npos){
		indx = chromosome.find("chr");
		if(indx != string::npos){
			chromosome = chromosome.substr(indx+3,chromosome.size()-3);
		}
	}
	else{
		chromosome = chromosome.substr(indx+5,chromosome.size()-5);
	}
	return chromosome;
}

void SNPOnChr::readSNPs(string fname) {
	SNPFile = fopen(fname.c_str(),"r");
	if(!(SNPFile = fopen(fname.c_str(),"r"))){
		cerr << "can not open SNP file " << fname << endl;
		exit(-1);
	}
	vector<string> chromosomes;
	vector<string>::iterator it;
	vector<vector<SNP> > AllSNP;
	vector<vector<SNP> >::iterator it2;
	snp_num = 0;
	char buf[1000];
	char * elems[10];
	int elemnum;
	long line_num = 0;
    	while(fgets(buf, 1000, SNPFile)){
		line_num++;
		// here assume the snp file is tab-delimited, every line being:
		// SNP id, chromosome, position, letters, strand, c_ref
		elemnum = split(buf, '\t', elems);
		if(elemnum == 6){
			snp_num++;
			//SNP snp(elems[0], atoll(elems[2]), elems[3], *elems[4], *elems[5]);
			SNP snp("", atoll(elems[2]), elems[3], *elems[4], *elems[5]);
			string chromosome = elems[1];
			chromosome = aberOfChr(chromosome);
			for(it = chromosomes.begin(); it < chromosomes.end(); it++){
				if(!chromosome.compare(*it)){
					break;
				}
			}
			if(it < chromosomes.end()){
				it2 = AllSNP.begin();
				it2 += (it-chromosomes.begin());
				(*it2).push_back(snp);
			}
			else{
				vector<SNP> v;
				v.push_back(snp);
				AllSNP.push_back(v);
				chromosomes.push_back(chromosome);
				//cerr << chromosome << endl;
			}
			
		}
		else{
			cerr << "Warning: malformed snp file " << fname <<
                    		", there should be 6 fields @line " << line_num << endl;
		    	cerr << buf << endl;
		    	//exit(1);
        	}
    	}
	for(size_t i = 0; i < chromosomes.size(); i++){
		this->insert(make_pair(chromosomes[i],AllSNP[i]));
	}
	
}

void SNPOnChr::readSNPsFromVCF(string fname) {
	SNPFile = fopen(fname.c_str(),"r");
	if(!(SNPFile = fopen(fname.c_str(),"r"))){
		cerr << "can not open VCF file " << fname << endl;
		exit(-1);
	}
	vector<string> chromosomes;
	vector<string>::iterator it;
	vector<vector<SNP> > AllSNP;
	vector<vector<SNP> >::iterator it2;
	
	float quality_th = 20;
	//float quality_th = 0;
	int depth_th = 10;
	snp_num = 0;
	char buf[20000];
	char * elems[10];
	int elemnum;
	long line_num = 0;
	int wrong_num = 0;
    	while(fgets(buf, 20000, SNPFile)){
		line_num++;
		if(buf[0] == '#') {
			continue;		
		}
        	elemnum = split(buf, '\t', elems);
		if(elemnum < 8){
			cerr << "Warning: malformed VCF file " << fname <<
                    		", there should be at least 8 fields @line " << line_num << endl;
            		cerr << buf << endl;
			wrong_num++;
			if(wrong_num > 10) {
				exit(1);
			}
			continue;
			//exit(1);
		}
		string info = elems[7];
		if(info.find("INDEL") != string::npos) {
			continue;		
		}
		if(strlen(elems[3]) > 1 || strlen(elems[4]) > 1) {
			continue;
		}
		size_t indx = info.find("DP=");
		if(indx != string::npos) {
			size_t indx1 = info.find(";", indx);
			int depth = atoi(info.substr(indx+3, indx1-indx-3).c_str());
			//cerr << depth << endl;
			if(depth < depth_th)
				continue;
		}
		float quality = atof(elems[5]);
        	if(quality < quality_th) {
			continue;
		}
		snp_num++;
		//SNP snp(elems[2], atoll(elems[1]), *elems[3], *elems[4]);
		SNP snp("", atoll(elems[1]), *elems[3], *elems[4]);
		string chromosome = elems[0];
		chromosome = aberOfChr(chromosome);
		for(it = chromosomes.begin(); it < chromosomes.end(); it++){
			if(!chromosome.compare(*it)){
				break;
			}
		}
		if(it < chromosomes.end()){
			it2 = AllSNP.begin();
			it2 += (it-chromosomes.begin());
			(*it2).push_back(snp);
		}
		else{
			vector<SNP> v;
			v.push_back(snp);
			AllSNP.push_back(v);
			chromosomes.push_back(chromosome);
			//cerr << chromosome << endl;
		}
    	}
	for(size_t i = 0; i < chromosomes.size(); i++){
		this->insert(make_pair(chromosomes[i],AllSNP[i]));
	}
	
}
