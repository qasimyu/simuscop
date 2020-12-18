// ***************************************************************************
// Config.cpp (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.
#include <random>
#include <iostream>
#include <fstream>
#include <cassert>

#include "MyDefine.h"
#include "split.h"
#include "Config.h"

Config::Config() {
	configFile = "";
	string strParaNames[] = {"bam", "profile", "ref", "variation", "snp", "vcf", "target",
							"bases", "output", "abundance", "layout", "samtools"};

	/*---start default configuration---*/

	int n = sizeof(strParaNames)/sizeof(string);
	for(int i = 0; i < n; i++) {
		if(strParaNames[i].compare("layout") == 0) {
			stringParas.insert(make_pair(strParaNames[i], "SE"));
		}
		else if(strParaNames[i].compare("bases") == 0) {
			stringParas.insert(make_pair(strParaNames[i], "ACTG"));
		}
		else {
			stringParas.insert(make_pair(strParaNames[i], ""));
		}
	}
	intParas.insert(make_pair("kmer", 0));
	intParas.insert(make_pair("bins", 0));
	intParas.insert(make_pair("threads", 1));
	intParas.insert(make_pair("verbose", 1));
	intParas.insert(make_pair("readLength", 0));
	intParas.insert(make_pair("coverage", 0));
	intParas.insert(make_pair("ploidy", 2));
	intParas.insert(make_pair("insertSize", 350));

	realParas.insert(make_pair("indelRate", 0.00025));
	/*---end default configuration---*/
}

void Config::loadConfig() {
	if(configFile.empty()) {
		cerr << "Error: configuration file not specified!" << endl;
		exit(-1);
	}
	ifstream ifs;
	ifs.open(configFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: can not open configuration file" << configFile << endl;
		exit(-1);
	}

	string line;
	int lineNum = 0, indx;
	while(getline(ifs, line)) {
		lineNum++;
		line = trim(line);
		if(line[0] == '#' || line.empty()) {
			continue;
		}
		indx = line.find('=');
		if(indx == string::npos) {
			cerr << "ERROR: line " << lineNum << " is incorrectly formatted in file " << configFile << endl;
			cerr << line << endl;
			exit(1);
		}
		string key = trim(line.substr(0, indx));
		string value = trim(line.substr(indx+1));

		if(stringParas.find(key) != stringParas.end()) {
			stringParas[key] = value;
		}
		else if(intParas.find(key) != intParas.end()) {
			intParas[key] = atoi(value.c_str());
		}
		else if(realParas.find(key) != realParas.end()) {
			realParas[key] = atof(value.c_str());
		}
		else if(key.compare("name") == 0) {
			popuNames = split(value, ',');
			for(int i = 0; i < popuNames.size(); i++) {
				popuNames[i] = trim(popuNames[i]);
			}
		}
		else {
			cerr << "ERROR: unrecognized item \"" << key << "\" @line " << lineNum << " in file " << configFile << endl;
			cerr << line << endl;
			exit(1);
		}
	}
	ifs.close();

	checkParas();
}

void Config::checkParas() {
	if(stringParas["profile"].empty()) {
		cerr << "Error: sequencing profile must be specified!" << endl;
		exit(1);
	}
	if(stringParas["snp"].empty()) {
		cerr << "Warning: SNP file not specified!" << endl;
		cerr << "No SNPs will be inserted into the genome." << endl;
	}
	if(stringParas["variation"].empty()) {
		cerr << "Warning: variation file not specified!" << endl;
		cerr << "No variations will be inserted into the genome." << endl;
	}
	if(stringParas["ref"].empty()) {
		cerr << "Error: reference file not specified!" << endl;
		exit(1);
	}
	if(popuNames.empty()) {
		cerr << "Error: population names not specified!" << endl;
		exit(1);
	}
	if(popuNames.size() > 1 && stringParas["abundance"].empty()) {
		cerr << "Error: abundance file not specified!" << endl;
		exit(1);
	}
	if(stringParas["output"].empty()) {
		cerr << "Error: output directory not specified!" << endl;
		exit(1);
	}

	if(stringParas["layout"].empty()) {
		cerr << "Warning: sequence layout not specified!" << endl;
		cerr << "use the default value: \"single end\"" << endl;
		stringParas["layout"] = "SE";
	}
	else if(stringParas["layout"].compare("SE") != 0 && stringParas["layout"].compare("PE") != 0) {
		cerr << "Error: sequence layout incorrectly specified!" << endl;
		cerr << "should be SE or PE" << endl;
		exit(1);
	}

	if(intParas["threads"] < 1) {
		cerr << "Error: number of threads should be a positive integer!" << endl;
		exit(1);
	}
	/*
	if(intParas["readLength"] < 1 || intParas["readLength"] > 100000) {
		cerr << "Error: read length should be a positive integer with maximum value of 100000!" << endl;
		exit(1);
	}
	*/
	if(intParas["coverage"] < 1) {
		cerr << "Error: sequence coverage should be a positive integer!" << endl;
		exit(1);
	}
	if(intParas["ploidy"] < 1) {
		cerr << "Error: genome ploidy should be a positive integer!" << endl;
		exit(1);
	}
	if(stringParas["layout"].compare("PE") == 0 && intParas["insertSize"] < intParas["readLength"]) {
		cerr << "Error: insert size should be not smaller than read length!" << endl;
		exit(1);
	}

	if(realParas["indelRate"] < 0 || realParas["indelRate"] > 0.001) {
		cerr << "Error: indel error rate should be a value between 0 to 0.001!" << endl;
		exit(1);
	}

	string tmp = stringParas["outputDir"];
	if(tmp[tmp.length()-1] == '/') {
		tmp.erase(tmp.end()-1);
		stringParas["outputDir"] = tmp;
	}
}

string Config::getStringPara(string paraName) {
	if(stringParas.find(paraName) != stringParas.end()) {
		return stringParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
	return "";
}

void Config::setStringPara(string paraName, string value) {
	if(stringParas.find(paraName) != stringParas.end()) {
		stringParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

int Config::getIntPara(string paraName) {
	if(intParas.find(paraName) != intParas.end()) {
		return intParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setIntPara(string paraName, int value) {
	if(intParas.find(paraName) != intParas.end()) {
		intParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

double Config::getRealPara(string paraName) {
	if(realParas.find(paraName) != realParas.end()) {
		return realParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setRealPara(string paraName, double value) {
	if(realParas.find(paraName) != realParas.end()) {
		realParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}
