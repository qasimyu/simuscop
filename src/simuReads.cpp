// ***************************************************************************
// simuRead.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <ctime>
#include <sys/stat.h>

#include "MyDefine.h"
#include "Config.h"
#include "Genome.h"
#include "Profile.h"
#include "ThreadPool.h"
#include "split.h"

void usage(const char* app);

int main(int argc, char *argv[]) {
	if(argc == 1) {
		cerr << "Error: configuration file is required!" << endl;
		usage(argv[0]);
		exit(1);
	}
	if(argc > 2) {
		cerr << "Error: too many input arguments!" << endl;
		usage(argv[0]);
        exit(1);
	}
	
	/*** record elapsed time ***/
	time_t start_t, end_t;
	long time_used;
	start_t = time(NULL);
	
	/*** fetch absolute path of the application ***/
	string binPath = argv[0];
	size_t indx = binPath.find_last_of('/');
	binPath = binPath.substr(0, indx);
	
	/*** config file ***/
	string configFile = argv[1];
	
	/*** read configuration file ***/
	config.setConfigFile(configFile);
	config.loadConfig();
	config.binPath = binPath;
	
	/*** create genome data ***/
	genome.loadData();

	/*** create output directory ***/
	string outputDir = config.getStringPara("output");
	string cmd = "test ! -e "+outputDir+" && mkdir -m 755 -p "+outputDir;
	system(cmd.c_str());
	
	/*** create thread pool ***/
	threadPool = new ThreadPool(config.getIntPara("threads"));
	threadPool->pool_init();
	//cerr << "\nThreadPool init done!" << endl;
	
	/*** load model ***/
	profile.train(config.getStringPara("profile"));
	
	/*** sampling reads ***/
	genome.generateSegments();
	//cerr << "\nsegments generation done!" << endl;
	genome.yieldReads();
	cerr << "\nReads generation done!" << endl;
	
	/*** remove temporary files ***/
	//cmd = "rm -R "+outputDir+"/tmp";
	//system(cmd.c_str());
	
	end_t = time(NULL);
	time_used = end_t-start_t;
	int minutes = time_used/60;
	int seconds = time_used - minutes*60;
	cerr << "\nElapsed time: " << minutes << " minutes and " << seconds << " seconds!\n" << endl;
	
	return 0;
}

void usage(const char* app) {
	cerr << "\nVersion: " << current_version << endl << endl;	
	cerr << "Usage: " << app << " <configuration file>" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " /path/to/config.txt" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

