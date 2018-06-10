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
#include <algorithm>
#include <sys/stat.h>
#include <sys/resource.h>
#include <math.h>
#include <set>
#include <vector>

using namespace std;

void loadNetwork(string fileIn, float pidentLimit, long double evalueLimit, map<string, unsigned short int>& positionsList, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, map<unsigned long long int, geneInfo>& genes, map<string, unsigned long long int>& seqNewName, map<unsigned long long int, string>& seqRealName, map<string, unsigned long long int>& gene2genome, set<unsigned long long int>& genomeList, unsigned short int verbose)
{
	unsigned short int numbering = 1;
	// Get headers position value
	// qseqid position
	unsigned short int qseqid_p = positionsList["qseqid"];
	// sseqid position
	unsigned short int sseqid_p = positionsList["sseqid"];
	// pident position
	unsigned short int pident_p = positionsList["pident"];
	// evalue position
	unsigned short int evalue_p = positionsList["evalue"];
	// qstart position
	unsigned short int qstart_p = positionsList["qstart"];
	// qend position
	unsigned short int qend_p = positionsList["qend"];
	// qlen position
	unsigned short int qlen_p = positionsList["qlen"];
	// sstart position
	unsigned short int sstart_p = positionsList["sstart"];
	// send position
	unsigned short int send_p = positionsList["send"];
	// slen position
	unsigned short int slen_p = positionsList["slen"];
	// numbering
	//map<string,unsigned long long int> newName;
	unsigned long long int num = 0;
	// Read input file
	ifstream blastp(fileIn.c_str());
	// Check file
	if(!blastp)
	{
		cout << "Can't open file " << fileIn << endl;
	}
	else
	{
		// Read blastp out and create the network
		string hit;
		// Start reading blastp out
		while(getline(blastp,hit))
		{
			// Get hit line
			istringstream hitInfo(hit);
			// Split line and store values into hitValues vector
			vector<string> hitValues;
			copy(istream_iterator<string>(hitInfo),istream_iterator<string>(),back_inserter<vector<string> >(hitValues));
			// Values to keep
			unsigned long long int qseqid, sseqid;
			unsigned short int qstart, qend, qlen, sstart, send, slen;
			float pident;
			long double evalue;
			// Get only the necessary values from the headers
			// Do not keep self hits and hits out of the pident thresholds
			if(numbering != 0)
			{
				if(seqNewName.find(hitValues[qseqid_p]) != seqNewName.end())
				{
					qseqid=seqNewName[hitValues[qseqid_p]];
				}
				else
				{
					num++;
					seqNewName[hitValues[qseqid_p]]=num;
					seqRealName[num]=hitValues[qseqid_p];
					qseqid=num;
				}
				if(seqNewName.find(hitValues[sseqid_p]) != seqNewName.end())
				{
					sseqid=seqNewName[hitValues[sseqid_p]];
				}
				else
				{
					num++;
					seqNewName[hitValues[sseqid_p]]=num;
					seqRealName[num]=hitValues[qseqid_p];
					sseqid=num;
				}
			}
			else
			{
				qseqid=atoll(hitValues[qseqid_p].c_str());
				sseqid=atoll(hitValues[sseqid_p].c_str());
			}
			if(qseqid != sseqid and atof(hitValues[pident_p].c_str()) >= pidentLimit)
			{
				// format hits in order to have the smallest sequence id first
				if(qseqid < sseqid)
				{
					qstart = atoi(hitValues[qstart_p].c_str());
					qend = atoi(hitValues[qend_p].c_str());
					qlen = atoi(hitValues[qlen_p].c_str());
					sstart = atoi(hitValues[sstart_p].c_str());
					send = atoi(hitValues[send_p].c_str());
					slen = atoi(hitValues[slen_p].c_str());
					pident = atof(hitValues[pident_p].c_str());
					evalue = atof(hitValues[evalue_p].c_str());
				}
				else if(qseqid > sseqid)
				{
					unsigned long long int tmp = qseqid;
					qseqid = sseqid;
					sseqid = tmp;
                               		qstart = atoi(hitValues[sstart_p].c_str());
                                	qend = atoi(hitValues[send_p].c_str());
					qlen = atoi(hitValues[slen_p].c_str());
                                	sstart = atoi(hitValues[qstart_p].c_str());
                                	send = atoi(hitValues[qend_p].c_str());
					slen = atoi(hitValues[qlen_p].c_str());
                                	pident = atof(hitValues[pident_p].c_str());
                               		evalue = atof(hitValues[evalue_p].c_str());
				}
				//cout << evalue << endl;
				// Free hitValues memory
				hitValues.clear();
				// Fill up edges map
				if( edges.find(make_pair(qseqid,sseqid)) != edges.end() )
				{
					// Get edge values for qseqid-sseqid hit
					edgeValues oldValues = edges[make_pair(qseqid,sseqid)];
					// Check if current e-value is better than old one.
					if( oldValues.getEvalue() > evalue )
					{
						// Update edgeValues for qseqid-sseqid hit [ keep best hit ]
						edges[make_pair(qseqid,sseqid)] = edgeValues(qstart, qend, sstart, send, pident, evalue);
					}
				}
				else
				{
					// Add new edge
					edges[make_pair(qseqid,sseqid)] = edgeValues(qstart, qend, sstart, send, pident, evalue);
					genes[qseqid]=geneInfo(qlen,0,seqRealName[qseqid],gene2genome[seqRealName[qseqid]],false);
                                       	genes[sseqid]=geneInfo(slen,0,seqRealName[sseqid],gene2genome[seqRealName[qseqid]],false);
                                       	genes[qseqid].insertNeighbor(sseqid);
                                       	genes[sseqid].insertNeighbor(qseqid);
                                       	genomeList.insert(gene2genome[seqRealName[qseqid]]);
                                       	genomeList.insert(gene2genome[seqRealName[qseqid]]);
				}
			}
		}
	}
	// Print Network information
	if( verbose == 1 )
	{
		cout << "---------- NETWORK INFO ----------" << endl << endl;
		cout << "Number of nodes        :\t" << genes.size() << endl << endl;
		cout << "Number of edges        :\t" << edges.size() << endl << endl;
		cout << "Number of genomes      :\t" << genomeList.size() << endl << endl;
		cout << "----------------------------------" << endl << endl;
	}
	seqRealName.clear();
        gene2genome.clear();
}
//end
