// ***************************************************************************
// Genome.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <unistd.h>
#include <algorithm>

#include "split.h"
#include "MyDefine.h"
#include "Genome.h"

void Genome::loadData() {
	loadAbers();
	/*
	loadCNVs();
	loadSNVs();
	*/
	loadSNPs(config.getStringPara("snp"), 1);
	loadRefSeq(config.getStringPara("ref"));
	loadTargets(outTargets, config.getStringPara("target"));
	divideTargets(outTargets);
	loadAbundance();
	
	curChr = "nodefined";
}

void Genome::loadData(string vcfFile, string refFile, string targetFile) {
	loadSNPs(vcfFile, 0);
	loadRefSeq(refFile);
	if(!targetFile.empty()) {
		loadTargets(inTargets, targetFile);
		divideTargets(inTargets);
	}
	curChr = "nodefined";
}

void Genome::loadAbers() {
	string aberFile = config.getStringPara("variation");
	if(aberFile.empty()) {
		return;
	}
	vector<string>& popuNames = config.getPopuNames();
	vector<string>::iterator v_it;
	map<string, vector<CNV> >::iterator m1_it;
	map<string, vector<SNV> >::iterator m2_it;
	map<string, vector<Insert> >::iterator m3_it;
	map<string, vector<Del> >::iterator m4_it;
	
	ifstream ifs;
	ifs.open(aberFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open file " << aberFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	int cnvCount = 0, snvCount = 0, insertCount = 0, delCount = 0;
	while(getline(ifs,line)) {
		lineNum++;
		if(line.empty()) {
			continue;
		}
		//cerr << line << "\n\t";
		vector<string> fields = split(line,'\t');
		
		string aberType = fields[0];
		//CNV
		if(aberType.compare("c") == 0) {
			if(fields.size() != 7) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string popuName = fields[1];
			v_it = find(popuNames.begin(), popuNames.end(), popuName);
			if(v_it == popuNames.end()) {
				cerr << "ERROR: unrecognized population identifier at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[2]);
			CNV cnv;
			cnv.spos = atol(fields[3].c_str());// 1-based
			cnv.epos = atol(fields[4].c_str());// 1-based
			cnv.CN = atof(fields[5].c_str());
			cnv.mCN = atof(fields[6].c_str());
			if(cnv.CN < cnv.mCN) {
				cerr << "ERROR: total copy number should be not lower than major copy number at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			if(cnv.CN-cnv.mCN > cnv.mCN) {
				cnv.mCN = cnv.CN-cnv.mCN;
			}
			PopuData<CNV>& cnvsOfPopu = cnvs[popuName];
			for(m1_it = cnvsOfPopu.begin(); m1_it != cnvsOfPopu.end(); m1_it++) {
				if((*m1_it).first.compare(chr) == 0) {
					(*m1_it).second.push_back(cnv);
					break;
				}
			}
			if(m1_it == cnvsOfPopu.end()) {
				vector<CNV> temp;
				temp.push_back(cnv);
				cnvsOfPopu.insert(make_pair(chr, temp));
			}
			cnvCount++;
		}
		//SNV
		else if(aberType.compare("s") == 0) {
			if(fields.size() != 7) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string popuName = fields[1];
			v_it = find(popuNames.begin(), popuNames.end(), popuName);
			if(v_it == popuNames.end()) {
				cerr << "ERROR: unrecognized population identifier at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[2]);
			SNV snv;
			snv.pos = atol(fields[3].c_str());// 1-based
			snv.ref = fields[4].at(0);
			snv.alt = fields[5].at(0);
			if(snv.ref == snv.alt) {
				cerr << "ERROR: the mutated allele should be not same as the reference allele at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			snv.type = fields[6];
			if(snv.type.compare("homo") != 0 && snv.type.compare("het") != 0) {
				cerr << "ERROR: unrecognized SNV type at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			PopuData<SNV>& snvsOfPopu = snvs[popuName];
			for(m2_it = snvsOfPopu.begin(); m2_it != snvsOfPopu.end(); m2_it++) {
				if((*m2_it).first.compare(chr) == 0) {
					(*m2_it).second.push_back(snv);
					break;
				}
			}
			if(m2_it == snvsOfPopu.end()) {
				vector<SNV> temp;
				temp.push_back(snv);
				snvsOfPopu.insert(make_pair(chr, temp));
			}
			snvCount++;
		}
		//Insert
		else if(aberType.compare("i") == 0) {
			if(fields.size() != 5) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string popuName = fields[1];
			v_it = find(popuNames.begin(), popuNames.end(), popuName);
			if(v_it == popuNames.end()) {
				cerr << "ERROR: unrecognized population identifier at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[2]);
			Insert insert;
			insert.pos = atol(fields[3].c_str());// 1-based
			insert.seq = fields[4];
			PopuData<Insert>& insertsOfPopu = inserts[popuName];
			for(m3_it = insertsOfPopu.begin(); m3_it != insertsOfPopu.end(); m3_it++) {
				if((*m3_it).first.compare(chr) == 0) {
					(*m3_it).second.push_back(insert);
					break;
				}
			}
			if(m3_it == insertsOfPopu.end()) {
				vector<Insert> temp;
				temp.push_back(insert);
				insertsOfPopu.insert(make_pair(chr, temp));
			}
			insertCount++;
		}
		//Del
		else if(aberType.compare("d") == 0) {
			if(fields.size() != 5) {
				cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string popuName = fields[1];
			v_it = find(popuNames.begin(), popuNames.end(), popuName);
			if(v_it == popuNames.end()) {
				cerr << "ERROR: unrecognized population identifier at line " << lineNum << " in file " << aberFile << endl;
				cerr << line << endl;
				exit(1);
			}
			string chr = abbrOfChr(fields[2]);
			Del del;
			del.pos = atol(fields[3].c_str());// 1-based
			del.len = atoi(fields[4].c_str());
			PopuData<Del>& delsOfPopu = dels[popuName];
			for(m4_it = delsOfPopu.begin(); m4_it != delsOfPopu.end(); m4_it++) {
				if((*m4_it).first.compare(chr) == 0) {
					(*m4_it).second.push_back(del);
					break;
				}
			}
			if(m4_it == delsOfPopu.end()) {
				vector<Del> temp;
				temp.push_back(del);
				delsOfPopu.insert(make_pair(chr, temp));
			}
			delCount++;
		}
		else {
			cerr << "ERROR: unrecognized aberraton type at line " << lineNum << " in file " << aberFile << endl;
			cerr << line << endl;
			exit(1);
		}
	}
	ifs.close();
	cerr << "\nDetails of the aberrations loaded from file " << aberFile << " are as follows:" << endl;
	cerr << "CNV: " << cnvCount << endl;
	cerr << "SNV: " << snvCount << endl;
	cerr << "Insert: " << insertCount << endl;
	cerr << "Del: " << delCount << endl;
}

void Genome::loadSNPs(string snpFile, int n) {
	if(n == 0) {
		scReal.readSNPsFromVCF(snpFile);
		cerr << "\n" << scReal.SNPNumber() << " filtered germline SNPs were loaded from file " << snpFile << endl;
	}
	else {
		if(snpFile.empty()) {
			return;
		}
		scSim.readSNPs(snpFile);
		cerr << "\n" << scSim.SNPNumber() << " SNPs to simulate were loaded from file " << snpFile << endl;
	}
}

void Genome::loadRefSeq(string refFile) {
	string tmp = refFile;
	if(tmp.empty()) {
		cerr << "genome sequence file not specified!" << endl;
		exit(1);
	}
	if(strcmp(tmp.substr(tmp.length()-3).c_str(), ".gz") == 0) {
		string cmd = "gzip -cd " + tmp + " > "+tmp.substr(0, tmp.length()-3);
		system(cmd.c_str());
		tmp = tmp.substr(0, refFile.length()-3);
	}
	fr.open(tmp, false);
	chromosomes = (*(fr.index)).sequenceNames;
	if(chromosomes.empty()) {
		cerr << "ERROR: reference sequence cannot be empty!" << endl;
		exit(1);
	}
	cerr << "\nReference sequence was loaded from file " << refFile << endl;
}

void Genome::loadTargets(map<string, vector<Target> >& targets, string targetFile) {
	if(targetFile.empty()) {
		return;
	}
	map<string, vector<Target> >::iterator m_it;
	
	ifstream ifs;
	ifs.open(targetFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open target file " << targetFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	int targetNum = 0;
	while(getline(ifs,line)) {
		lineNum++;
		//cerr << line << "\n\t";
		vector<string> fields = split(line,'\t');
		if(fields.size() < 3) {
			cerr << "ERROR: line " << lineNum << " should have at least 3 fields in file " << targetFile << endl;
			cerr << line << endl;
			exit(1);
		}
		
        string chr = abbrOfChr(fields[0]);
		long chrLen = getChromLen(chr);
		if(chrLen <= 0) {
			continue;
		}
		Target target;
		target.spos = max((long) 1, atol(fields[1].c_str())-50+1); // 0-based 
		//target.epos = min(chrLen, atol(fields[2].c_str())+50); // 1-based
		long tmp;
		if(atol(fields[2].c_str()) <= 0) {
			tmp = chrLen-(-atol(fields[2].c_str()))%chrLen;
		}
		else {
			tmp = atol(fields[2].c_str());
		}
		target.epos = min(chrLen, tmp+50); // 1-based
		for(m_it = targets.begin(); m_it != targets.end(); m_it++) {
			if((*m_it).first.compare(chr) == 0) {
				(*m_it).second.push_back(target);
				break;
			}
		}
		if(m_it == targets.end()) {
			vector<Target> temp;
			temp.push_back(target);
			targets.insert(make_pair(chr, temp));
		}
		targetNum++;
	}
	ifs.close();
	cerr << "\ntotal " << targetNum << " targets were loaded from file " << targetFile << endl;
}

void Genome::loadAbundance() {
	string abundanceFile = config.getStringPara("abundance");
	if(abundanceFile.empty()) {
		return;
	}
	
	vector<string>& popuNames = config.getPopuNames();
	
	ifstream ifs;
	ifs.open(abundanceFile.c_str());
	if(!ifs.is_open()) {
		cerr << "can not open abundance file " << abundanceFile << endl;
		exit(-1);
	}
	
	string line;
	int lineNum = 0;
	while(getline(ifs,line)) {
		lineNum++;
		//cerr << line << "\n\t";
		vector<string> fields = split(line,'\t');
		if(fields.size() != popuNames.size()) {
			cerr << "ERROR: line " << lineNum << " has wrong number of fields in file " << abundanceFile << endl;
			cerr << line << endl;
			exit(1);
		}
        vector<float> props;
		float sum = 0;
		for(int i = 0; i < fields.size(); i++) {
			float prop = atof(fields[i].c_str());
			sum += prop;
			props.push_back(prop);
		}
		if(fabs(1-sum) > 0.001) {
			cerr << "ERROR: the sum of abundances is not equal to one at line " << lineNum << " in file " << abundanceFile << endl;
			cerr << line << endl;
			exit(1);
		}
		mixProps.push_back(props);
	}
	ifs.close();
	cerr << "\ntotal " << mixProps.size() << " population combinations were loaded from file " << abundanceFile << endl;
}

/*
void Genome::generateAltSequence() {
	int i, j;
	for(i = 0; i < chromosomes.size(); i++) {
		string chr = chromosomes[i];
		vector<SNP>& snpsOfChr = getSNPs(chr);
		string refSeq = fr.getSubSequence(chr, 0, getChromLen(chr));
		transform(refSeq.begin(), refSeq.end(), refSeq.begin(), (int (*)(int))toupper);
		for(j = 0; j < snpsOfChr.size(); j++) {
			SNP snp = snpsOfChr[j];
			long pos = snp.getPosition();
			refSeq[pos-1] = snp.getNucleotide();
		}
		altSequence.insert(make_pair(chr, refSeq));
	}
}
*/

vector<CNV>& Genome::getCNVs(string popu, string chr) {
	PopuData<CNV>& cnvsOfPopu = cnvs[popu];
	return cnvsOfPopu[chr];
}

vector<SNV>& Genome::getSNVs(string popu, string chr) {
	PopuData<SNV>& snvsOfPopu = snvs[popu];
	return snvsOfPopu[chr];
}

vector<Insert>& Genome::getInserts(string popu, string chr) {
	PopuData<Insert>& insertsOfPopu = inserts[popu];
	return insertsOfPopu[chr];
}

vector<Del>& Genome::getDels(string popu, string chr) {
	PopuData<Del>& delsOfPopu = dels[popu];
	return delsOfPopu[chr];
}

vector<Segment>& Genome::getSegs(string popu, string chr) {
	PopuData<Segment>& segsOfPopu = segments[popu];
	return segsOfPopu[chr];
}

long Genome::getChromLen(string chr) {
	FastaIndexEntry fie;
	int i;
	for(i = 0; i < chromosomes.size(); i++) {
		if(chromosomes[i].compare(chr) == 0) {
			break;
		}
	}
	if(i == chromosomes.size()) {
		return 0;
	}
	fie = fr.index->entry(chr);
	return fie.length;
}
long Genome::getGenomeLength() {
	long n = 0;
	for(int i = 0; i < chromosomes.size(); i++) {
		n += getChromLen(chromosomes[i]);
	}
	return n;
}

long Genome::getTargetLength() {
	if(!outTargets.empty())  {
		long n = 0;
		map<string, vector<Target> >::iterator it;
		for(it = outTargets.begin(); it != outTargets.end(); it++) {
			vector<Target>& targets = it->second;
			for(int j = 0; j < targets.size(); j++) {
				n += targets[j].epos-targets[j].spos+1;
			}
		}
		return n;
	}
	else {
		return getGenomeLength();
	}
}

char* Genome::getSubSequence(string chr, int startPos, int length) {
	string seq = fr.getSubSequence(chr, startPos, length);
	transform(seq.begin(), seq.end(), seq.begin(), (int (*)(int))toupper);
	char* s = new char[seq.length()+1];
	strcpy(s, seq.c_str());
	return s;
}

void Genome::generateAltSequence(string chr) {	
	int i, j;
	vector<string>::iterator it = find(chromosomes.begin(), chromosomes.end(), chr);
	if(it == chromosomes.end()) {
		cerr << "ERROR: unrecognized chromosome identifier \"" << chr 
			 << "\" when generating alternative sequence!" << endl;
		exit(1);
	}
	
	curChr = chr;
	vector<SNP>& snpsOfChr = getRealSNPs(chr);
	altSequence = fr.getSubSequence(chr, 0, getChromLen(chr));
	for(j = 0; j < snpsOfChr.size(); j++) {
		SNP snp = snpsOfChr[j];
		long pos = snp.getPosition();
		altSequence[pos-1] = snp.getNucleotide();
	}
	transform(altSequence.begin(), altSequence.end(), altSequence.begin(), (int (*)(int))toupper);
}


char* Genome::getSubAltSequence(string chr, int startPos, int length) {
	if(curChr.compare(chr) != 0) {
		generateAltSequence(chr);
	}
	//string seq = altSequence[chr].substr(startPos, length);
	string seq = altSequence.substr(startPos, length);
	char* s = new char[seq.length()+1];
	strcpy(s, seq.c_str());
	return s;
}


int Genome::getSegReadCount(string popu, string chr, int segIndx) {
	PopuData<Segment>& segsOfPopu = segments[popu];
	vector<Segment>& chrSegs = segsOfPopu[chr];
	if(segIndx < 0 || segIndx > chrSegs.size()-1) {
		return 0;
	}
	return chrSegs[segIndx].getReadCount();
}

int Genome::getSegIndxWithOffset(string popu, string chr, int curSegIndx, int offset) {
	PopuData<Segment>& segsOfPopu = segments[popu];
	vector<Segment>& chrSegs = segsOfPopu[chr];
	int i;
	if(offset == 0) {
		return curSegIndx;
	}
	if((offset > 0) && (curSegIndx < chrSegs.size()-1)) {
		unsigned int n = 0;
		for(i = curSegIndx+1; i < chrSegs.size(); i++) {
			if(offset <= chrSegs[i].getSeqSize()+n) {
				break;
			}
			n += chrSegs[i].getSeqSize();
		}
		if(i < chrSegs.size()) {
			return i;
		}
	}
	if((offset < 0) && (curSegIndx > 0)) {
		unsigned int n = 0;
		offset = -offset;
		for(i = curSegIndx-1; i >= 0; i--) {
			if(offset <= chrSegs[i].getSeqSize()+n) {
				break;
			}
			n += chrSegs[i].getSeqSize();
		}
		if(i >= 0) {
			return i;
		}
	}
	return -1;
}

char* Genome::produceFragment(string popu, string chr, int startSegIndx, int segSeqIndx, int fragLen) {
	PopuData<Segment>& segsOfPopu = segments[popu];
	vector<Segment>& chrSegs = segsOfPopu[chr];
	if(fragLen <= 0) {
		return NULL;
	}
	
	int i, len = 0;
	char* fragSeq = new char[fragLen+1];
	char* segSeq;
	for(i = startSegIndx; i < chrSegs.size(); i++) {
		segSeq = chrSegs[i].getSegSequence(segSeqIndx);
		if(segSeq != NULL) {
			if(len+strlen(segSeq) >= fragLen) {
				strncpy(fragSeq+len, segSeq, fragLen-len);
				fragSeq[fragLen] = '\0';
				break;
			}
			else {
				strncpy(fragSeq+len, segSeq, strlen(segSeq));
				len += strlen(segSeq);
			}
		}
	}
	if(i == chrSegs.size()) {
		char *p = fragSeq;
		p[fragLen] = '\0';
		fragSeq = new char[len+1];
		strncpy(fragSeq, p, len);
		fragSeq[len] = '\0';
		delete[] p;
	}
	return fragSeq;
}

void Genome::generateSegments() {
	int i, j, k;
	vector<string>& popuNames = config.getPopuNames();
	int ploidy = config.getIntPara("ploidy");
	int mCN = (int) ceil((float) ploidy/2);
	
	if(!outTargets.empty()) {
		map<string, vector<Target> >::iterator it;
		chromosomes.clear();
		for(it = outTargets.begin(); it != outTargets.end(); it++) {
			string chr = it->first;
			for(i = 0; i < chromosomes.size(); i++) {
				if(chr.compare(chromosomes[i]) == 0) {
					break;
				}		
			}
			if(i == chromosomes.size()) {
				chromosomes.push_back(chr);
			}
		}
	}
	for(i = 0; i < popuNames.size(); i++) {
		string popu = popuNames[i];
		for(j = 0; j < chromosomes.size(); j++) {
			string chr = chromosomes[j];
			int segIndx = 0;
			vector<CNV>& cnvsOfChr = getCNVs(popu, chr);
			long segStartPos = 1, segEndPos = -1;
			
			for(k = 0; k < cnvsOfChr.size(); k++) {		
				if(segStartPos > getChromLen(chr)) {
					break;
				}
				CNV cnv = cnvsOfChr[k];
				cnv.epos = min(cnv.epos, getChromLen(chr));
				if(segStartPos < cnv.spos) {
					divideSegment(popu, chr, segStartPos, cnv.spos-1, ploidy, mCN, segIndx);
				}
				divideSegment(popu, chr, cnv.spos, cnv.epos, cnv.CN, cnv.mCN, segIndx);
				segStartPos = cnv.epos+1;
			}
			if(segStartPos <= getChromLen(chr)) {
				segEndPos = getChromLen(chr);
				divideSegment(popu, chr, segStartPos, getChromLen(chr), ploidy, mCN, segIndx);
			}
		}
	}
	
}

void Genome::divideTargets(map<string, vector<Target> >& allTargets) {
	int i, j, k;
	unsigned int targetMaxSize = Segment::fragSize;
	map<string, vector<Target> > newTargets;
	map<string, vector<Target> >::iterator it, m_it;
	for(it = allTargets.begin(); it != allTargets.end(); it++) {
		string chr = it->first;
		vector<Target>& targets = it->second;
		for(j = 0; j < targets.size(); j++) {
			Target target = targets[j];
			long spos = target.spos;
			long tsize = target.epos-target.spos+1;
			int k = tsize/targetMaxSize;
			for(i = 0; i < k; i++) {
				Target newTarget;
				newTarget.spos = spos;
				if(i == k-1) {
					newTarget.epos = target.epos;
				}
				else {
					newTarget.epos = spos+targetMaxSize-1;
				}
				spos = newTarget.epos+1;
				for(m_it = newTargets.begin(); m_it != newTargets.end(); m_it++) {
					if((*m_it).first.compare(chr) == 0) {
						(*m_it).second.push_back(newTarget);
						break;
					}
				}
				if(m_it == newTargets.end()) {
					vector<Target> temp;
					temp.push_back(newTarget);
					newTargets.insert(make_pair(chr, temp));
				}
			}
			if(spos <= target.epos) {
				Target newTarget;
				newTarget.spos = spos;
				newTarget.epos = target.epos;
				for(m_it = newTargets.begin(); m_it != newTargets.end(); m_it++) {
					if((*m_it).first.compare(chr) == 0) {
						(*m_it).second.push_back(newTarget);
						break;
					}
				}
				if(m_it == newTargets.end()) {
					vector<Target> temp;
					temp.push_back(newTarget);
					newTargets.insert(make_pair(chr, temp));
				}	
			}
		}
	}
	allTargets = newTargets;
	
}

void Genome::divideSegment(string popu, string chr, long segStartPos, long segEndPos, int CN, int mCN, int &segIndx) {
	unsigned int segMaxSize = Segment::getSegMaxSize();
	PopuData<Segment>& segsOfPopu = segments[popu];
	long segSize = segEndPos-segStartPos+1;
	int n = segSize/segMaxSize;
	unsigned int m = segSize-n*segMaxSize;
	for(int i = 0; i < n; i++) {
		if(i == n-1 && m < segMaxSize/2) {
			Segment seg(segIndx++, chr, segStartPos, segEndPos, CN, mCN);
			segsOfPopu[chr].push_back(seg);
			segStartPos = segEndPos+1;
		}
		else {
			Segment seg(segIndx++, chr, segStartPos, segStartPos+segMaxSize-1, CN, mCN);
			segsOfPopu[chr].push_back(seg);
			segStartPos += segMaxSize;
		}
	}
	if(segStartPos <= segEndPos) {
		Segment seg(segIndx++, chr, segStartPos, segEndPos, CN, mCN);
		segsOfPopu[chr].push_back(seg);
	}
}

void Genome::calculateACNs(map<string, double>& ACNs) {
	int i;
	map<string, PopuData<Segment> >::iterator it1;
	map<string, vector<Segment> >::iterator it2;
	for(it1 = segments.begin(); it1 != segments.end(); it1++) {
		PopuData<Segment>& segsOfPopu = it1->second;
		long sum = 0;
		for(it2 = segsOfPopu.begin(); it2 != segsOfPopu.end(); it2++) {
			vector<Segment>& chrSegs = it2->second;
			for(i = 0; i < chrSegs.size(); i++) {
				sum += chrSegs[i].getSeqSize();
			}
		}
		double acn = (double) sum/getGenomeLength();
		ACNs.insert(make_pair(it1->first, acn));
	}
}

void Genome::setReadCounts(string popu, long reads) {
	int i, j, k;
	PopuData<Segment>& segsOfPopu = segments[popu];
	map<string, double> chrWLens;
	double WL = 0;
	for(i = 0; i < chromosomes.size(); i++) {
		string chr = chromosomes[i];
		vector<Segment>& chrSegs = segsOfPopu[chr];
		double chrWL = 0;
		for(j = 0; j < chrSegs.size();j++) {
			chrWL += chrSegs[j].getWeightedLength();
		}
		WL += chrWL;
		chrWLens.insert(make_pair(chr, chrWL));
	}
	
	long curReads = 0;
	long chrReads;
	for(i = 0; i < chromosomes.size(); i++) {
		string chr = chromosomes[i];
		vector<Segment>& chrSegs = segsOfPopu[chr];
		double chrWL = chrWLens[chr];
		//cerr << chr << ": " << chrWL << endl;
		if(i < chromosomes.size()-1) {
			chrReads = reads*(chrWL/WL);
		}
		else {
			chrReads = reads-curReads;
		}
		long sum = 0;
		for(j = 0; j < chrSegs.size();j++) {
			if(j < chrSegs.size()-1) {
				long segReadCount = (chrSegs[j].getWeightedLength()/chrWL)*chrReads;
				chrSegs[j].setReadCount(segReadCount);
				sum += segReadCount;
			}
			else {
				chrSegs[j].setReadCount(chrReads-sum);
			}
		}
		curReads += chrReads;
	}
}

void Genome::yieldReads() {	
	int i, j, k, m, n = 0;
	vector<string>& popuNames = config.getPopuNames();
	
	long reads = getTargetLength()*config.getIntPara("coverage")/config.getIntPara("readLength");
	if(config.isVerbose()) {	
		cerr << "\nNumber of reads to sample: " << reads << endl;
	}
	
	map<string, double> ACNs;
	map<string, double>::iterator it;
	calculateACNs(ACNs);
	if(config.isVerbose()) {
		if(ACNs.size() > 1) {
			cerr << "\nAverage copy number of populations: " << endl;
		}
		else {
			cerr << "\nAverage copy number: " << endl;
		}
		for(it = ACNs.begin(); it != ACNs.end(); it++) {
			cerr << it->first << ": " << it->second << endl;
		}
	}
	
	cerr << "\n*****Generating samples*****" << endl;
	srand(time(0));
	
	if(mixProps.empty()) {
		curPopu = popuNames[0];
		
		string fqFilePrefix = config.getStringPara("outputDir")+"/";
		if(config.isPairedEnd()) {
			string outFile1 = fqFilePrefix+popuNames[0]+"_1.fq";
			string outFile2 = fqFilePrefix+popuNames[0]+"_2.fq";
			swp = new SeqWriter(outFile1, outFile2);
		}
		else {
			string outFile = fqFilePrefix+popuNames[0]+".fastq";
			swp = new SeqWriter(outFile);
		}
		
		setReadCounts(popuNames[0], reads);
		PopuData<Segment>& segsOfPopu = segments[popuNames[0]];
		for(j = 0; j < chromosomes.size(); j++) {
			string chr = chromosomes[j];
			curChr = chr;
			
			vector<Segment>& chrSegs = segsOfPopu[chr];
			
			for(k = 0; k < chrSegs.size(); k++) {
				chrSegs[k].generateSegSequences(); 
			}
			for(k = 0; k < chrSegs.size(); k++) {
				//Segment::yieldReads(&chrSegs[k]);
				threadPool->pool_add_work(&Segment::yieldReads, &chrSegs[k], n++); 
			}
			threadPool->wait();
			for(k = 0; k < chrSegs.size(); k++) {
				chrSegs[k].clearSegSequences();
			}
		}
		delete swp;
	}
	else {
		for(m = 0; m < mixProps.size(); m++) {
			vector<float> props = mixProps[m];
			double w_acn = 0;
			char buf[1000];
			string fn = "";
			if(config.isVerbose()) {
				cerr << "\ncurrent mixed proportions: " << endl;
			}
			for(i = 0; i < popuNames.size(); i++) {
				if(config.isVerbose()) {
					cerr << props[i] << '\t';
				}
				string popu = popuNames[i];
				float prop = props[i];
				double acn = ACNs[popu];
				w_acn += prop*acn;
				if(i == 0) {
					sprintf(buf, "%s_%.3f", popu.c_str(), prop);
				}
				else {
					sprintf(buf, "+%s_%.3f", popu.c_str(), prop);
				}
				fn += buf;
			}
			if(config.isVerbose()) {
				cerr << "\nfile name:" << endl;
				cerr << fn << endl;
			}
			
			string fqFilePrefix = config.getStringPara("outputDir")+"/";
			if(config.isPairedEnd()) {
				string outFile1 = fqFilePrefix+fn+"_pe_1.fastq";
				string outFile2 = fqFilePrefix+fn+"_pe_2.fastq";
				swp = new SeqWriter(outFile1, outFile2);
			}
			else {
				string outFile = fqFilePrefix+fn+".fastq";
				swp = new SeqWriter(outFile);
			}
			
			for(i = 0; i < popuNames.size(); i++) {
				string popu = popuNames[i];
				curPopu = popu;
				
				long popuReads = reads*props[i]*ACNs[popu]/w_acn;
				setReadCounts(popu, popuReads);
				PopuData<Segment>& segsOfPopu = segments[popu];
				for(j = 0; j < chromosomes.size(); j++) {
					string chr = chromosomes[j];
					curChr = chr;
					
					vector<Segment>& chrSegs = segsOfPopu[chr];
					
					for(k = 0; k < chrSegs.size(); k++) {
						chrSegs[k].generateSegSequences();
					}
					for(k = 0; k < chrSegs.size(); k++) {
						//Segment::yieldReads(&chrSegs[k]);
						threadPool->pool_add_work(&Segment::yieldReads, &chrSegs[k], n++); 
					}
					threadPool->wait();
					for(k = 0; k < chrSegs.size(); k++) {
						chrSegs[k].clearSegSequences();
					}
				}
			}
			delete swp;
		}
	}
}


