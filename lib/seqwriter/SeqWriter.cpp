// ***************************************************************************
// SeqWriter.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdlib>
#include <iostream>

#include "SeqWriter.h"


SeqWriter::SeqWriter(string seFile) {
	se_ofs.open(seFile.c_str());
	if(!se_ofs.is_open()) {
		cerr << "Error: can not open fastq file to save results:\n" << seFile << endl;
		exit(-1);
	}
	pe = false;
	pthread_mutex_init(&pm, NULL);
}

SeqWriter::SeqWriter(string peFile1, string peFile2) {
	pe_ofs1.open(peFile1.c_str());
	if(!pe_ofs1.is_open()) {
		cerr << "Error: can not open fastq file to save results:\n" << peFile1 << endl;
		exit(-1);
	}
	pe_ofs2.open(peFile2.c_str());
	if(!pe_ofs2.is_open()) {
		cerr << "Error: can not open fastq file to save results:\n" << peFile2 << endl;
		exit(-1);
	}
	pe = true;
	pthread_mutex_init(&pm, NULL);
}

SeqWriter::~SeqWriter() {
	close();
}

void SeqWriter::write(char* str) {
	//printf("thread 0x%x go in\n", pthread_self());
	pthread_mutex_lock(&pm);
	se_ofs << str;
	pthread_mutex_unlock(&pm);
	//printf("thread 0x%x go out\n", pthread_self());
}

void SeqWriter::write(char* str1, char* str2) {
	pthread_mutex_lock(&pm);
	pe_ofs1 << str1;
	pe_ofs2 << str2;
	pthread_mutex_unlock(&pm);
}

void SeqWriter::close() {
	if(pe){
		pe_ofs1.close();
		pe_ofs2.close();
	}
	else {
		se_ofs.close();
	}
}
