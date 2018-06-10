/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        EGN2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/


#include "functions.h"
#include <map>
#include <string>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <list>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <cfloat>
#include <cstring>
#include <vector>
#include <iterator>
#include <set>

using namespace std;

// Get EGN2 main options from main.config file
void getMainOptions(string configDir, long double& evalueLimit, float& pidentLimit, unsigned short int& minCov, unsigned short int& mutual, unsigned int& nCpu, unsigned short int& verbose)
{
	string mainConfigFile = configDir + "/main.config"; 
	ifstream lines(mainConfigFile.c_str());
	string line;
	while(getline(lines,line))
	{
		if(line.find("#") != std::string::npos || line.empty())
		{
			//Do nothing
		}
		else
		{
			stringstream lineInfo(line);
                        // Split line and store values into lineValues vector
                        vector<string> lineValues;
                        copy(istream_iterator<string>(lineInfo),istream_iterator<string>(),back_inserter<vector<string> >(lineValues));
			// Get options
			if(lineValues[0] == "evalueLimit")
			{
				evalueLimit = atof(lineValues[1].c_str());
			}
			else if(lineValues[0] == "pidentLimit")
			{
				pidentLimit = atof(lineValues[1].c_str());
			}
			else if(lineValues[0] == "minCov")
			{
				minCov = atoi(lineValues[1].c_str());
			}
			else if(lineValues[0] == "mutual")
			{
				mutual = atoi(lineValues[1].c_str());
			}
			else if(lineValues[0] == "nCpu")
			{
				nCpu = atoi(lineValues[1].c_str());
			}
			else if(lineValues[0] == "verbose")
			{
				verbose = atoi(lineValues[1].c_str());
			}
		}
	}
}
//END
