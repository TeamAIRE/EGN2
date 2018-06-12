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


void renameFasta(string fileIn, string outDir, map<string, unsigned long long int>& seqNewName, map<unsigned long long int, string>& seqRealName)
{
	string exonerateDir = getenv("EGN2EXONERATE");
	string fileRenamed = outDir + "/" + fileIn + ".renamed";
	ofstream outputRenamed(fileRenamed.c_str());
	unsigned long long int i = 0;
	ifstream lines(fileIn.c_str());
	string line;
	while(getline(lines,line))
	{
		if(line.find(">") != std::string::npos)
		{
			int r=0;
			i++;
			line.erase(line.begin());
			outputRenamed << ">" << i << endl;
			seqNewName[line]=i;
			seqRealName[i]=line;
		}
		else
		{
			outputRenamed << line << endl;
		}
	}
	system((exonerateDir + "/fastareformat -f " + fileRenamed + " > " + fileRenamed + ".tmp").c_str());
	system(("rm " + fileRenamed).c_str());
	system(("mv " + fileRenamed + ".tmp " + fileRenamed).c_str());
	cout << "Nb of Sequences: " << i << endl; 
}
//END
