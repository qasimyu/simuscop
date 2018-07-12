// ***************************************************************************
// SeqWriter.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _SEQWRITER_H
#define _SEQWRITER_H

#include <fstream>
#include <string>
#include <pthread.h>

using namespace std;

class SeqWriter {
	private:
		ofstream se_ofs;
		ofstream pe_ofs1, pe_ofs2;
		
		bool pe;
		pthread_mutex_t pm;
		
	public:		
		SeqWriter(string seFile);
		SeqWriter(string peFile1, string peFile2);
		~SeqWriter();
		
		void write(char* str);
		void write(char* str1, char* str2);
		
		void close();
		
};


#endif

