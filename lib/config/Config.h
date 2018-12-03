// ***************************************************************************
// Config.h (c) 2018 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _CONFIG_H
#define _CONFIG_H

#include <vector>
#include <map>
#include <cstdlib>
#include <string>

using namespace std;

class Config {
	private:
		string configFile;
		
		map<string, string> stringParas;
		map<string, int> intParas;
		map<string, double> realParas;
		vector<string> popuNames;
	public:
		Config();
		Config(string configFile) : configFile(configFile){}
		~Config() {}
		
		string binPath;
		
		void setConfigFile(string configFile) {
			this->configFile = configFile;
		}
		
		bool isVerbose() {return !(intParas["verbose"] == 0);}
		bool isPairedEnd() {return stringParas["layout"].compare("PE") == 0;}
		void loadConfig();
		void checkParas();
		
		string getStringPara(string paraName);
		void setStringPara(string paraName, string value);
		int getIntPara(string paraName);
		void setIntPara(string paraName, int value);
		double getRealPara(string paraName);
		void setRealPara(string paraName, double value);
		
		vector<string>& getPopuNames() {return popuNames;}
		
};


#endif

