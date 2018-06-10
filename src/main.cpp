/*

	Written by Jananan PATHMANATHAN, 2017-2018
	
	This file is part of EGN2.
	
	EGN2 is shared under Creative commons licence: 
	
	Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
	
	See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "parseHeader.h"
#include "functions.h"
#include "getOptions.h"
#include "renameFasta.h"
#include "attributes.h"
#include "similarity.h"
#include "networks.h"
#include "createNetwork.h"
#include "loadNetwork.h"
#include "getTime.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <list>
#include <iterator>
#include <unistd.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>
#include <thread>
#include <cmath>
#include <ctime>
#include <set>
#include <vector>

using namespace std;

int main(int argc, char * argv[])
{
	// ALL TIME (start)
	time_t Astart = time(NULL);
	// *************** Default Values ***************
	// Blastp file header
	// Clean Blastp file header
	string header = "qseqid sseqid qstart qend qlen sstart send slen pident evalue";
	// Evalue threshold
	long double evalueLimit;
	long double& refEvalueLimit = evalueLimit; 
	// Percentage of identity threshold
	float pidentLimit;
	float& refPidentLimit = pidentLimit;
	// Minimum coverage between two sequences
	unsigned short int minCov;
	unsigned short int& refMinCov = minCov;
	// Mutual Coverage
	unsigned short int mutual;
	unsigned short int& refMutual = mutual;
	// Number of CPU(s)
	unsigned int nCpu;
	unsigned int& refNcpu = nCpu;
	// Verbose mode
	unsigned short int verbose;
	unsigned short int& refVerbose = verbose;
	// Check cmd line 
	if(argc < 1 )
	{
		//help();
		cout << endl;
                cout << "----- Usage of ENG2 -----" << endl << endl;
		cout << "Case 1: You do not have computed the similarity between sequences." << endl;
		cout << "./runEGN2 -f file.fasta -s [nucl or prot] -m [blast or diamond] -a attributes.file" << endl << endl;
		cout << "Case 2: You have computed the similarity between sequences." << endl;
                cout << "./runEGN2 -n similarity_file -a attributes.file -h header" << endl << endl;
                cout << "-- optional -- " << endl;
		cout << "-h header : " << endl;
                cout << "-d workDir" << endl << endl;
		exit(1);
	}
	// Time
	time_t start, end;
	int c;
	// sequence fasta file
	string fileIn;
	// FusedTriplets output file
	string fileOut;
	// Working Directory
	string workDir = get_current_dir_name();
	// Config file Directory
	string configDir = getenv("EGN2CONFIGDIR");
	// Attribute file
	string attrFile;
	// blast, blat or diamond
	string simTool;
	// Sequences type: nucl (nucleotide) or prot (protein)
	string seqType;
	// Similarity file
	string simFile;
	// Get arguments
	while( (c = getopt (argc, argv, "f:s:m:a:d:n:h:")) != -1)
	{
		switch(c)
		{
			case 'f':
				fileIn = optarg;
				break;
			case 's':
				seqType = optarg;
				break;
			case 'm':
				simTool = optarg;
				break;
			case 'a':
				attrFile = optarg;
				break;
			case 'd':
				workDir = optarg;
				break;
			case 'n':
                                simFile = optarg;
                                break;
			case 'h':
				header = optarg;
				break;
			default:
				cout << endl;
		                cout << "----- Usage of ENG2 -----" << endl << endl;
        		        cout << "Case 1: You do not have computed the similarity between sequences." << endl;
                		cout << "./runEGN2 -f file.fasta -s [nucl or prot] -m [blast or diamond] -a attributes.file" << endl << endl;
                		cout << "Case 2: You have computed the similarity between sequences." << endl;
                		cout << "./runEGN2 -n similarity_file -a attributes.file -h header" << endl << endl;
                		cout << "-- optional -- " << endl;
                		cout << "-h header"<< endl;
                		cout << "-d workDir" << endl << endl;
				exit(1);
				break;
		}
	}
	//
	unsigned short int STOP = 0;
	if( fileIn == "" && simFile == "" )
	{
		cout << "Please provide a fasta file [-f] or a similarity file [-n] !" << endl << endl;
		cout << "----- Usage of ENG2 -----" << endl << endl;
                cout << "Case 1: You do not have computed the similarity between sequences." << endl;
                cout << "./runEGN2 -f file.fasta -s [nucl or prot] -m [blast or diamond] -a attributes.file" << endl << endl;
                cout << "Case 2: You have computed the similarity between sequences." << endl;
                cout << "./runEGN2 -n similarity_file -a attributes.file -h header" << endl << endl;
                cout << "-- optional -- " << endl;
                cout << "-h header" << endl;
                cout << "-d workDir" << endl << endl;
                exit(1);
	}
	if( seqType == "" && simFile == "" )
        {
        	cout << "Sequence type [-s] value is empty !" << endl;
                cout << "Please select one of those : nucl (nucleotide) or prot (protein)" << endl;
		STOP = 1;
        }
        else if( seqType != "nucl" && seqType != "prot" && simFile == "" )
        {
            	cout << "Sequence type [-s] " << seqType << " unknown." << endl;
               	cout << "Please select one of those : nucl (nucleotide) or prot (protein)" << endl;
		STOP = 1;
        }
        if( simTool == "" && simFile == "" )
        {
		cout << "Similarity tool [-m] value is empty !" << endl;
		cout << "Please select one of those : blast or diamond" << endl;
		STOP = 1;
        }
        else if( simTool != "blast" && simTool != "diamond" && simFile == "" )
        {
		cout << "Similarity tool [-m] value " << simTool << " unknown." << endl;
		cout << "Please select one of those : blast, blat or diamond" << endl;
		STOP = 1;
        }
	if(attrFile == "")
	{
		cout << "Please provide the attributes file [-a] !" << endl;
		STOP = 1;
	}
	if(STOP == 1)
	{
		cout << endl;
		cout << "----- Usage of ENG2 -----" << endl << endl;
                cout << "Case 1: You do not have computed the similarity between sequences." << endl;
                cout << "./runEGN2 -f file.fasta -s [nucl or prot] -m [blast or diamond] -a attributes.file" << endl << endl;
                cout << "Case 2: You have computed the similarity between sequences." << endl;
                cout << "./runEGN2 -n similarity_file -a attributes.file -h header" << endl << endl;
                cout << "-- optional -- " << endl;
                cout << "-h header" << endl;
                cout << "-d workDir" << endl << endl;
                exit(1);
	}
	// input file basename
	string basename;
	if(fileIn != "" && simFile == "")
	{
		basename = fileIn.substr(fileIn.find_last_of("\\/")+1);
	}
	else if(fileIn == "" && simFile != "")
        {
                basename = simFile.substr(fileIn.find_last_of("\\/")+1);
        }
	string filename = basename;
	replace(basename.begin(), basename.end(), '.', '_');
	string folder = "egn2_" + basename;
	// Output directory
	string folderTime = asctime(localtime(&Astart));
	replace(folderTime.begin(), folderTime.end(), ' ', '_');
	replace(folderTime.begin(), folderTime.end(), ':', '_');
	folderTime.pop_back();
	string outputDir = workDir + "/" + folder + "_" + folderTime;
	system(("mkdir " + outputDir).c_str());
        // Ouput file basename
        fileOut = outputDir + "/" + basename;
	// Get EGN2 main parameters
	getMainOptions(configDir, refEvalueLimit, refPidentLimit, refMinCov, refMutual, refNcpu, refVerbose);
	// copy config file 
	system(("mkdir " + outputDir + "/configfiles").c_str());
	system(("cp " + configDir + "/*.config " + outputDir + "/configfiles/").c_str());
	// Print running command line
	if( verbose == 1 && simFile == "")
	{
		cout << "---------- PARAMETERS ----------" << endl;
		cout << "Input         : " << fileIn << endl;
		cout << "Ouput         : " << outputDir << endl;
		cout << "Sequence      : " << seqType << endl;
		cout << "SimTool       : " << simTool << endl;
		cout << "E-value       : " << evalueLimit << endl;
		cout << "Pident        : " << pidentLimit << endl;
		cout << "MinCov        : " << minCov << endl;
		cout << "Mutual        : " << mutual << endl;
		cout << "Nb CPUs       : " << nCpu << endl;
		cout << "Verbose       : " << verbose << endl;
		cout << "--------------------------------" << endl << endl;
	}
	else if(verbose == 1 && simFile != "")
	{
		cout << "---------- PARAMETERS ----------" << endl;
                cout << "Input         : " << simFile << endl;
                cout << "Ouput         : " << outputDir << endl;
                cout << "E-value       : " << evalueLimit << endl;
                cout << "Pident        : " << pidentLimit << endl;
                cout << "MinCov        : " << minCov << endl;
                cout << "Mutual        : " << mutual << endl;
                cout << "Verbose       : " << verbose << endl;
                cout << "--------------------------------" << endl << endl;
	}
	// Create all maps and their reference
	// edges[query,subject]:qstart, qenid, sstart, send, evalue and pident
	map<pair<unsigned long long int, unsigned long long int>, edgeValues> edges;
	map<pair<unsigned long long int, unsigned long long int>, edgeValues>& refEdges = edges;
	// all_neighbors[node]: list of all neighbors
	map<unsigned long long int, geneInfo > genes;
	map<unsigned long long int, geneInfo >& refGenes = genes;
	// Sequences new name
	map<unsigned long long int, string> seqRealName;
        map<unsigned long long int, string>& refSeqRealName = seqRealName;
	map<string, unsigned long long int> seqNewName;
        map<string, unsigned long long int>& refSeqNewName = seqNewName;
	// Check header
	map<string, unsigned short int> positionsList;
	map<string, unsigned short int>& refPositionsList = positionsList;
	positionsList = getPositions(header);
	// gene to genome dico
	map<unsigned long long int, unsigned long long int> gene2genome;
	map<unsigned long long int, unsigned long long int>& refGene2genome = gene2genome;
	// gene to genome dico2 if sim file given by user
	map<string, unsigned long long int> gene2genome2;
        map<string, unsigned long long int>& refGene2genome2 = gene2genome2;
	// genome info
	map<string, unsigned long long int > genomeNewName;
        map<string, unsigned long long int >& refGenomeNewName = genomeNewName;
	map<unsigned long long int, genomeInfo> genomes;
	map<unsigned long long int, genomeInfo>& refGenomes = genomes;
	map<unsigned long long int, string> genomeRealName;
	map<unsigned long long int, string>& refGenomeRealName = genomeRealName;
	set<unsigned long long int> genomeList;
	set<unsigned long long int>& refGenomeList = genomeList;
	// Connected Compnents
	map<unsigned long long int, ccInfo> CCs;
	map<unsigned long long int, ccInfo>& refCCs = CCs;
	// Step 1 : Rename Sequences
	if(simFile == "")
	{
		start = time(NULL);
		if( verbose == 1 )
       		{
       	        	cout << "----- RENAME SEQUENCES -----" << endl << endl;
                	cout << "START - " << asctime(localtime(&start)) << endl;
        	}
		renameFasta(fileIn, outputDir, refSeqNewName, refSeqRealName);
		end = time(NULL);
		if( verbose == 1 )
        	{
                	cout << "END   - " << asctime(localtime(&end)) << endl;
        	}
		duration( start, end, verbose );
        	cout << endl;
		// Step 2 : Get Attributes
		getAttributes(attrFile, refGene2genome, refSeqNewName, refGenomes, refGenomeNewName);
		// Step 3 : Compute Similarity
		start = time(NULL);
		if( verbose == 1 )
		{
			cout << "----- COMPUTE SIMILARITY -----" << endl << endl;
			cout << "START - " << asctime(localtime(&start)) << endl;
		}
		//
		string renamedFasta = outputDir + "/" + filename + ".renamed";
		runSimilarity(renamedFasta, simTool, seqType, configDir, outputDir, nCpu, folderTime);
		end = time(NULL);
		if( verbose == 1 )
		{
			cout << "END   - " << asctime(localtime(&end)) << endl;
		}
		duration( start, end, verbose );
		cout << endl;
		// Step 4 : Create Network
		start = time(NULL);
        	if( verbose == 1 )
        	{
        	        cout << "----- NETWORK CREATION -----" << endl << endl;
               	 	cout << "START - " << asctime(localtime(&start)) << endl;
                	cout << endl;
        	}
		system(("mkdir " + outputDir + "/cleanNetwork").c_str());
		string fileOutCN = outputDir + "/cleanNetwork/" + filename + ".renamed.cleanNetwork";
		string simFile = outputDir + "/similarityRes/" + filename + ".renamed.out";
        	createNetwork(simFile, fileOutCN, pidentLimit, evalueLimit, refPositionsList, refEdges, refGenes, refSeqRealName, refGene2genome, refGenomeList, verbose);
        	end = time(NULL);
        	if( verbose == 1 )
        	{
        	        cout << "END   - " << asctime(localtime(&end)) << endl;
        	}
        	// Get time
        	duration( start, end, verbose );
        	cout << endl;
	}
	else
	{
		getAttributes2(attrFile, refGene2genome2, refGenomes, refGenomeNewName);
		start = time(NULL);
		if( verbose == 1 )
                {
                        cout << "----- LOADING NETWORK -----" << endl << endl;
                        cout << "START - " << asctime(localtime(&start)) << endl;
                        cout << endl;
                }
                loadNetwork(simFile, pidentLimit, evalueLimit, refPositionsList, refEdges, refGenes, refSeqNewName, refSeqRealName, refGene2genome2, refGenomeList, verbose);
                end = time(NULL);
                if( verbose == 1 )
                {
                        cout << "END   - " << asctime(localtime(&end)) << endl;
                }
                // Get time
                duration( start, end, verbose );
                cout << endl;
	}
	// Step 5: Connecnted Components
	string outputDirNet = outputDir + "/networks";
        system(("mkdir " + outputDirNet).c_str());
	start = time(NULL);
        if( verbose == 1 )
        {
                cout << "----- COMPUTE CONNECTED COMPONENTS -----" << endl << endl;
                cout << "START - " << asctime(localtime(&start)) << endl;
        }
        computeCC(refGenes, refEdges, minCov, mutual, refCCs, refGenomes, outputDirNet, refGenomeList, verbose);
        end = time(NULL);
        if( verbose == 1 )
        {
                cout << "END   - " << asctime(localtime(&end)) << endl;
        }
        duration( start, end, verbose );
        cout << endl;
	// Step 6: Genome and CC networks
	start = time(NULL);
        if( verbose == 1 )
        {
                cout << "----- COMPUTE GENOME AND CC NETWORKS -----" << endl << endl;
                cout << "START - " << asctime(localtime(&start)) << endl;
        }
	computeGenomeAndCCnetwork(refGenes, refEdges, refCCs, refGenomes, outputDirNet, verbose);
        end = time(NULL);
        if( verbose == 1 )
        {
                cout << "END   - " << asctime(localtime(&end)) << endl;
        }
        duration( start, end, verbose );
        cout << endl;
	// ALL TIME (end)
	time_t Aend = time(NULL);
	// REAL TIME
	if( verbose == 1 )
	{
		cout << "----- TOTAL COMPUTING TIME -----" << endl;
	}
	duration( Astart, Aend, verbose);
	cout << endl;
}
//END
