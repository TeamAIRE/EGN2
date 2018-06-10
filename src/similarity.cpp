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
#include <algorithm>
#include <sys/stat.h>
#include <sys/resource.h>
#include <math.h>
#include <set>
#include <vector>
#include <thread>

using namespace std;

// Compute similarity
void runBlast(string seqType, string subFile, string dbFile, string simToolOptions, unsigned int l)
{
	// Get similarity tool parameters
	string blastDir = getenv("EGN2BLAST");
	string simCmd;
	if(seqType == "prot")
	{
		simCmd = blastDir + "/blastp -query " + subFile + " -db " + dbFile + " " + simToolOptions + " -out " + subFile + ".out";
	}
	else if(seqType == "nucl")
	{
		simCmd = blastDir + "/blastn -query " + subFile + " -db " + dbFile + " " + simToolOptions + " -out " + subFile + ".out";
	}
	system((simCmd).c_str());
}


void runDiamond(string seqType, string subFile, string dbFile, string simToolOptions, unsigned int l)
{
        // Get similarity tool parameters
	string diamondDir = getenv("EGN2DIAMOND");
	string simCmd;
        if(seqType == "prot")
        {
                simCmd = diamondDir + "/diamond blastp --query " + subFile + " --db " + dbFile + " " + simToolOptions + " --out " + subFile + ".out --quiet";
        }
        else if(seqType == "nucl")
        {
                simCmd = diamondDir + "/diamond blastn --query " + subFile + " --db " + dbFile + " " + simToolOptions + " --out " + subFile + ".out --quiet";
        }
        system((simCmd).c_str());
}

// Load Network (MultiThread)

void runSimilarity( string fileIn, string simTool, string seqType, string configDir, string outDir, unsigned int nCpu, string timeInfo)
{
	//
	string blastDir = getenv("EGN2BLAST");
	string diamondDir = getenv("EGN2DIAMOND");
	string exonerateDir = getenv("EGN2EXONERATE");
	string simDir = outDir + "/similarityRes";
	system(("mkdir " + simDir).c_str());
        // split network
	string tmpDir = outDir + "/tmp";
	system(("mkdir " + tmpDir).c_str());
	string basename = fileIn.substr(fileIn.find_last_of("\\/")+1);
        if( nCpu > 1 )
        {
                // Split file into N parts where N = nCpu
                string splitCmd = exonerateDir + "/fastasplit -f " + fileIn + " -o " + tmpDir + " -c " + to_string(nCpu);
                system((splitCmd).c_str());
		//cout << splitCmd << endl;
		//exit(1);
                // Rename splited files
                string modNameCmd = "for i in $(ls " + tmpDir + "/); do N=$(($N+1)); mv " + tmpDir + "/$i "+ tmpDir + "/$N." + simTool + ".fasta_" + timeInfo + "; done";
                //
                system((modNameCmd).c_str());
        }
        else
        {
                // create link
                system(("ln -s " + fileIn + " " + tmpDir + "/1." + simTool + ".fasta_" + timeInfo).c_str());
        }
	//exit(1);
	// Get similarity tool options
	string simToolConfig;
	if(simTool == "blast")
	{
		if(seqType == "prot")
		{
			// Get config file
			simToolConfig = configDir + "/blastp.config";
			// Create db
			system((blastDir + "/makeblastdb -in " + fileIn + " -dbtype prot &> /dev/null").c_str());
		}
		else if(seqType == "nucl")
		{
			// Get config file
                        simToolConfig = configDir + "/blastn.config";
			// Create db
			system((blastDir + "/makeblastdb -in " + fileIn + " -dbtype nucl &> /dev/null").c_str());
		}
	}
	else if(simTool == "diamond")
	{
		// Get config file
		simToolConfig = configDir + "/diamond.config";
		// Create db
		system((diamondDir + "/diamond makedb --in " + fileIn + " --db " + fileIn + ".diamondDB &> /dev/null").c_str());
	}
	string simToolOptions;
        ifstream parameters(simToolConfig.c_str());
        string parameter;
        while(getline(parameters,parameter))
        {
                if(parameter.find("#") != std::string::npos)
                {       
                        // do nothing
                }
                else
                {
                        simToolOptions = simToolOptions + " " + parameter;   
                }
        }
	//cout << simToolOptions << endl;
	//exit(1);
        // Multithreading
        vector<thread> threads;
        for(unsigned int i = 1; i < nCpu+1; ++i)
        {
                string subFasta = tmpDir + "/" + to_string(i) + "." + simTool + ".fasta_" + timeInfo;
		if(simTool == "blast")
		{
                	threads.push_back(thread(runBlast, seqType, subFasta, fileIn, simToolOptions, i));
		}
		else if(simTool == "diamond")
		{
			threads.push_back(thread(runDiamond, seqType, subFasta, fileIn, simToolOptions, i));
		}
        }
        for(auto& thread : threads)
        {
                thread.join();
        }
        // Merge similarity tool output
	for( unsigned int i = 1; i < nCpu+1; i++)
        {
		string fileToMerge = tmpDir + "/" + to_string(i) + "." + simTool + ".fasta_" + timeInfo + ".out";
                system(("cat " + fileToMerge + " >> " + simDir + "/" + basename + ".out").c_str());
                system(("rm " + fileToMerge).c_str());
	}
        system(("rm -r " + tmpDir).c_str());
}
//END
