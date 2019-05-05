// ***************************************************************************
// seqToProfile.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <ctime>

#include "MyDefine.h"
#include "Genome.h"
#include "Profile.h"


void usage(const char* app);

int main(int argc, char *argv[]) {
	string bamFile = "", targetFile = "";
	string vcfFile = "", refFile = "";
	string outFile = "", samtools = "";
	
	int binCount = 50;
	int kmer = 3;

	/*** record elapsed time ***/
	time_t start_t, end_t;
	long time_used;
	start_t = time(NULL);

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"bam", required_argument, 0, 'b'},
		{"target", required_argument, 0, 't'},
		{"vcf", required_argument, 0, 'v'},
		{"ref", required_argument, 0, 'r'},
		{"output", required_argument, 0, 'o'},
		{"samtools", required_argument, 0, 's'},
		{"kmer", required_argument, 0, 'k'},
		{"bins", required_argument, 0, 'B'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hb:t:v:r:o:s:k:B:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage(argv[0]);
				exit(0);
			case 'b':
				bamFile = optarg;
				break;
			case 't':
				targetFile = optarg;
				break;
			case 'v':
				vcfFile = optarg;
				break;
			case 'r':
				refFile = optarg;
				break;
			case 'o':
				outFile = optarg;
				break;
			case 's':
				samtools = optarg;
				break;
			case 'k':
				kmer = atoi(optarg);
				break;
			case 'B':
				binCount = atoi(optarg);
				break;
			default :
				usage(argv[0]);
				exit(1);
		}
	}

	if(bamFile.empty()){
		cerr << "Use --bam to specify a normal BAM file." << endl;
		usage(argv[0]);
		exit(1);
	}

	/*
	if(targetFile.empty()){
		cerr << "\nsequencing strategy: WGS!" << endl;
	}
	else {
		cerr << "\nsequencing strategy: WES!" << endl;
	}
	*/
	
	if(vcfFile.empty()){
		cerr << "Use --vcf to specify the VCF file generated from the normal BAM." << endl;
		usage(argv[0]);
		exit(1);
	}

	if(refFile.empty()){
		cerr << "Use --ref to specify the reference file(.fasta) to which the reads are aligned." << endl;
		usage(argv[0]);
		exit(1);
	}

	if(samtools.empty()) {
		cerr << "\nWarning: the path of samtools not specified!" << endl;
		cerr << "Assume the tool has been installed and included in the system PATH!" << endl;
	}
	
	if(kmer < 1 || kmer > 5) {
		cerr << "Error: parameter \"kmer\" should be a positive integer with maximum value of 5!" << endl;
		exit(1);
	}
	
	if(binCount < 10) {
		cerr << "Error: parameter \"bins\" should be a positive integer with minimum value of 10!" << endl;
		exit(1);
	}
	
	config.setIntPara("kmer", kmer);
	config.setIntPara("bins", binCount);
	config.setStringPara("bam", bamFile);
	config.setStringPara("samtools", samtools);
	config.setStringPara("output", outFile);
	config.setStringPara("vcf", vcfFile);
	config.setStringPara("target", targetFile);
	config.setStringPara("ref", refFile);
	
	/*** load data ***/
	genome.loadTrainData();

	/*** profile learning ***/
	profile.init();
	profile.train();
	
	end_t = time(NULL);
	time_used = end_t-start_t;
	int minutes = time_used/60;
	int seconds = time_used - minutes*60;
	cerr << "\nElapsed time: " << minutes << " minutes and " << seconds << " seconds!\n" << endl;
	
	return 0;
}

void usage(const char* app) {
	cerr << "\nVersion: " << current_version << endl << endl;
	cerr << "Usage: " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -b, --bam <string>              normal BAM file" << endl
		<< "    -t, --target <string>           exome target file (.bed) for whole-exome sequencing[default:null]" << endl
		<< "    -v, --vcf <string>              the VCF file generated from the normal BAM" << endl
		<< "    -r, --ref <string>              genome reference file (.fasta) to which the reads were aligned" << endl
		<< "    -o, --output <string>           output file" << endl
		<< "    -s, --samtools <string>         the path of samtools [default:samtools]" << endl
		<< "    -k, --kmer <int>                the length of kmer sequence [default:3]" << endl
		<< "    -B, --bins <int>                the number of bins into which bases of read are grouped [default:50]" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -b /path/to/normal.bam -t /path/to/normal.bed -v /path/to/normal.vcf -r /path/to/ref.fa -k 4 > /path/to/results.model" << endl
		<< endl
		<< "    " << app << " -b /path/to/normal.bam -v /path/to/normal.vcf -r /path/to/ref.fa -o /path/to/results.model -s /path/to/samtools" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

