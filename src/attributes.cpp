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

void getAttributes(string attrFile, map<unsigned long long int, unsigned long long int>& gene2genome, map<string, unsigned long long int>& seqNewName, map<unsigned long long int, genomeInfo>& genomes, map<string, unsigned long long int>& genomeNewName)
{
	unsigned long long int i = 0;
	ifstream lines(attrFile.c_str());
	set<string> tmpGenomes;
	string line;
	while(getline(lines,line))
	{
		//
		replace(line.begin(), line.end(), ' ', '_');
		// Get hit line
                istringstream lineInfo(line);
                // Split line and store values into hitValues vector
                vector<string> lineValues;
                copy(istream_iterator<string>(lineInfo),istream_iterator<string>(),back_inserter<vector<string> >(lineValues));
		if(tmpGenomes.find(lineValues[1]) != tmpGenomes.end())
		{
			gene2genome[seqNewName[lineValues[0]]]=genomeNewName[lineValues[1]];
			genomes[genomeNewName[lineValues[1]]].addSeq();
		}
		else
		{
			tmpGenomes.insert(lineValues[1]);
			i++;
			genomes[i]=genomeInfo(1,lineValues[1]);
			genomeNewName[lineValues[1]]=i;
			gene2genome[seqNewName[lineValues[0]]]=genomeNewName[lineValues[1]];
		}
	}
	//cout << "Nb of Genomes: " << i << endl; 
}


// get attribute when user have already computed the similarity

void getAttributes2(string attrFile, map<string, unsigned long long int>& gene2genome, map<unsigned long long int, genomeInfo>& genomes, map<string, unsigned long long int>& genomeNewName)
{
        unsigned long long int i = 0;
        ifstream lines(attrFile.c_str());
        set<string> tmpGenomes;
        string line;
        while(getline(lines,line))
        {
		//
		replace(line.begin(), line.end(), ' ', '_');
                // Get hit line
                istringstream lineInfo(line);
                // Split line and store values into hitValues vector
                vector<string> lineValues;
                copy(istream_iterator<string>(lineInfo),istream_iterator<string>(),back_inserter<vector<string> >(lineValues));
                if(tmpGenomes.find(lineValues[1]) != tmpGenomes.end())
                {
                        gene2genome[lineValues[0]]=genomeNewName[lineValues[1]];
                        genomes[genomeNewName[lineValues[1]]].addSeq();
                }
                else
                {
                        tmpGenomes.insert(lineValues[1]);
                        i++;
                        genomes[i]=genomeInfo(1,lineValues[1]);
                        genomeNewName[lineValues[1]]=i;
                        gene2genome[lineValues[0]]=genomeNewName[lineValues[1]];
                }
        }
        //cout << "Nb of Genomes: " << i << endl;
}
//END
