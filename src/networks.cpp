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
#include <set>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

// Connected component detection using DFS algo

void computeCC( map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, unsigned short int minCov, unsigned short int mutual, map<unsigned long long int, ccInfo>& CCs,  map<unsigned long long int, genomeInfo>& genomes, string outputDir, set<unsigned long long int>& genomeList, unsigned short int verbose)
{
	// Number of connected components
	unsigned long long int nCC = 0;
	// cc info file
	string fileOutCCinfo = outputDir + "/CC.info";
	ofstream outputCCinfo(fileOutCCinfo.c_str());
	outputCCinfo << "#CCid\tnodes\tedges\tconnectivity\tmeanPident\tsdPident\tNbGenomes" << endl;
	// nodes file
	string fileOutCCnodes = outputDir + "/CC.nodes";
	ofstream outputCCnodes(fileOutCCnodes.c_str());
	// Gephi tmp node file
	string tmpGDFnodes = outputDir + "/tmp_nodes.gdf";
	ofstream outputGDFnodes(tmpGDFnodes.c_str());
	outputGDFnodes << "nodedef>name VARCHAR,label VARCHAR,genome VARCHAR,cc VARCHAR" << endl;
	// edges file
	string fileOutCCedges = outputDir + "/CC.edges";
	ofstream outputCCedges(fileOutCCedges.c_str());
	// Gephi tmp edge file
	string tmpGDFedges = outputDir + "/tmp_edges.gdf";
        ofstream outputGDFedges(tmpGDFedges.c_str());
        outputGDFedges << "edgedef>node1 VARCHAR,node2 VARCHAR,pident DOUBLE" << endl;
	// CC genome presence absence
        string fileOutCCinfop = outputDir + "/CC.genomeInfo";
        ofstream outputCCinfop(fileOutCCinfop.c_str());
        set<unsigned long long int>::iterator fc;
        outputCCinfop << "#nCC";
        for(fc=genomeList.begin();fc!=genomeList.end();++fc)
        {
        	outputCCinfop << "\t" << genomes[*fc].getGenomeRealName();
        }
        outputCCinfop << endl;
        // Compute connected components
	map<unsigned long long int, geneInfo>::iterator it;
	for( it = genes.begin(); it != genes.end(); it++)
	{
		// Check if node has been visited, if true go to next one
		if(genes[it->first].isVisited() ) continue;
		// Start
		// cc family number
		nCC++;
		// cc info
		CCs[nCC]=ccInfo(0);
		outputCCnodes << ">CC" << nCC << endl;
		outputCCedges << ">CC" << nCC << endl;
		// count cc edges
		unsigned long long int NbEdges = 0;
		// count cc nodes
		unsigned long long int NbNodes = 0;
		// Percentage of identity info
		float sumPident = 0.0;
		// vetor of pident for standard deviatin calculation
		vector<float> pidentList;
		vector<float>& refPidentList = pidentList;
		// Temporary list containing nodes belonging to CC
		set<unsigned long long int> tmpNodes;
		set<unsigned long long int>::iterator val;
		// Insert CC starting node
		tmpNodes.insert(it->first);
		// Get CC until tmpNodes is not empty
		while(!tmpNodes.empty())
		{
			// Get node
			val = tmpNodes.begin();
			unsigned long long int node = *val;
			// Delete node from tmpNodes
			tmpNodes.erase(val);
			// If not visited check node
			unsigned long long int i;
			if( !genes[node].isVisited() )
			{
				// Change node status to visited
				genes[node].setVisited(true);
				NbNodes++;
				genes[node].setCcID(nCC);
				CCs[nCC].insertGenome(genes[node].getGenomeID());
				outputCCnodes << node << endl;
				outputGDFnodes << node << "," << genes[node].getRealName() << "," << genomes[genes[node].getGenomeID()].getGenomeRealName() << "," << nCC<< endl;
				// Current node's neighnors
				vector<unsigned long long int> currentNodeNeighbors = genes[node].getNeighbors();
				while(!currentNodeNeighbors.empty())
				{
					// Get neighbor
					unsigned long long int neighbor = currentNodeNeighbors.back();
					currentNodeNeighbors.pop_back();
					// Check neighbor's status
					if(!genes[neighbor].isVisited())
					{
						if(checkCoverage(node, neighbor, minCov, genes, edges, mutual))
						{
							// Add neighbor to tmpNodes
							tmpNodes.insert(neighbor);
							// count edge
							NbEdges++;
							float pident = getEdgeValues(edges,node,neighbor).getPident();
							pidentList.push_back(pident);
							sumPident+=pident;
							outputCCedges << node << "\t" << neighbor << "\t" << pident << endl;
							outputGDFedges << node << "," << neighbor << "," << pident << endl;
						}
					}
				}
			}
		}
		// output info
		float connectivity = (2.0*(float)NbEdges)/((float)NbNodes*((float)NbNodes-1.0));
		float meanPident = sumPident/(float)NbEdges;
		float standardDev = computeSD(meanPident,refPidentList);
		CCs[nCC].setSize(NbNodes);
		if(NbNodes == 1)
		{
			outputCCinfo << "CC" << nCC << "\t" << NbNodes << "\t0\t1\t0\t0\t" << CCs[nCC].getNbGenomes() << endl;
		}
		else if(NbNodes == 2)
		{
			 outputCCinfo << "CC" << nCC << "\t" << NbNodes << "\t" << NbEdges << "\t1\t" << setprecision (2) << fixed << sumPident << "\t0\t" << CCs[nCC].getNbGenomes() << endl;
		}
		else
		{
                	outputCCinfo << "CC" << nCC << "\t" << NbNodes << "\t" << NbEdges << "\t" << setprecision (2) << fixed << connectivity << "\t" << setprecision (2) << fixed << meanPident << "\t" << setprecision (2) << fixed << standardDev << "\t" << CCs[nCC].getNbGenomes() << endl;
		}
		pidentList.clear();
		// CC genome presence absence
		set<unsigned long long int>::iterator fc;
		outputCCinfop << "CC" << nCC;
		set<unsigned long long int> tmp = CCs[nCC].getGenomes();
		for(fc=genomeList.begin();fc!=genomeList.end();++fc)
		{
			if(tmp.find(*fc) != tmp.end())
			{
				outputCCinfop << "\t1";
			}
			else
			{
				outputCCinfop << "\t0";
			}
		}
		outputCCinfop << endl;
	}
	// Clear isVisted map container
	if(verbose == 1)
	{
		cout << "Nb of CC : " << nCC << endl;
	}
	// Create final GDF file
	string finalGDF = outputDir + "/ConnetedComponents.gdf";
	system(("cat " + tmpGDFnodes + " >> " + finalGDF).c_str());
	system(("cat " + tmpGDFedges + " >> " + finalGDF).c_str());
	system(("rm " + tmpGDFnodes + " " + tmpGDFedges).c_str());
}


// Compute CC and Genome network
void computeGenomeAndCCnetwork( map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, map<unsigned long long int, ccInfo>& CCs, map<unsigned long long int, genomeInfo>& genomes, string outputDir, unsigned short int verbose)
{
	// cc network gdf
        string ccNetGDF = outputDir + "/ccNetwork.gdf";
        ofstream ccNetGDFoutput(ccNetGDF.c_str());
	// cc id
	set<unsigned long long int> ccNetNodes;
	// cc net edges
	map<pair<unsigned long long int, unsigned long long int>,unsigned long long int> ccNetEdges;
	// genome network gdf
	string genomeNetGDF = outputDir + "/GenomeSeqNetwork.gdf";
        ofstream genomeNetGDFoutput(genomeNetGDF.c_str());
	// genome id
	set<unsigned long long int> genomeNetNodes;
	// genome net edges
	map<pair<unsigned long long int, unsigned long long int>,unsigned long long int> genomeNetEdges;
	//
	map<pair<unsigned long long int,unsigned long long int>, edgeValues>::iterator it;
        for(it=edges.begin(); it!=edges.end(); ++it)
	{
		// gene net nodes
		unsigned long long int e1 = it->first.first;
		unsigned long long int e2 = it->first.second;
		// CC net
		unsigned long long int cc1 = genes[e1].getCcID();
		unsigned long long int cc2 = genes[e2].getCcID();
		if(cc1 != cc2)
		{
			ccNetNodes.insert(cc1);
			ccNetNodes.insert(cc2);
			if(cc1 < cc2)
			{
				ccNetEdges[make_pair(cc1,cc2)]+=1;
			}
			else
			{
				ccNetEdges[make_pair(cc2,cc1)]+=1;
			}
		}
		else
		{
			ccNetNodes.insert(cc1);
		}
		// Genome net
		unsigned long long int genome1 = genes[e1].getGenomeID();
		unsigned long long int genome2 = genes[e2].getGenomeID();
		//cout << "n1-" << e1 << "\tn2-" << e2 << " --- " << "g1-" << genome1 << "\tg2-" << genome2 << endl; 
		if(genome1 != genome2)
		{
			genomeNetNodes.insert(genome1);
			genomeNetNodes.insert(genome2);
			if(genome1 < genome2)
			{
				genomeNetEdges[make_pair(genome1,genome2)]+=1;
			}
			else
			{
				genomeNetEdges[make_pair(genome2,genome1)]+=1;
			}
		}
		else
		{
			genomeNetNodes.insert(genome1);
		}
	}	
	// Output cc net
	ccNetGDFoutput << "nodedef>name VARCHAR,label VARCHAR,nseq VARCHAR,ngenome VARCHAR" << endl;
	std::set<unsigned long long int>::iterator cc;
	for(cc=ccNetNodes.begin(); cc!=ccNetNodes.end(); ++cc)
	{
		ccNetGDFoutput << *cc << ",CC" << *cc << "," << CCs[*cc].getSize() << "," << CCs[*cc].getNbGenomes() << endl;
	}
	ccNetGDFoutput << "edgedef>node1 VARCHAR,node2 VARCHAR,nHit DOUBLE" << endl;
	map<pair<unsigned long long int, unsigned long long int>,unsigned long long int>::iterator cce;
	for(cce=ccNetEdges.begin(); cce!=ccNetEdges.end(); ++cce)
	{
		//cc net nodes
                unsigned long long int n1 = cce->first.first;
                unsigned long long int n2 = cce->first.second;
		unsigned long long int nHit = cce->second;
		ccNetGDFoutput << n1 << "," << n2 << "," << nHit << endl;
	}
	// Output genome net
	genomeNetGDFoutput << "nodedef>name VARCHAR,label VARCHAR,nseq VARCHAR" << endl;
	std::set<unsigned long long int>::iterator g;
	for(g=genomeNetNodes.begin(); g!=genomeNetNodes.end(); ++g)
	{
		genomeNetGDFoutput << *g << "," << genomes[*g].getGenomeRealName() << "," << genomes[*g].getSize() << endl;
	}
	genomeNetGDFoutput << "edgedef>node1 VARCHAR,node2 VARCHAR,nHit DOUBLE" << endl;
	map<pair<unsigned long long int, unsigned long long int>,unsigned long long int>::iterator ge;
	for(ge=genomeNetEdges.begin(); ge!=genomeNetEdges.end(); ++ge)
	{
		//genome net nodes
                unsigned long long int g1 = ge->first.first;
                unsigned long long int g2 = ge->first.second;
                unsigned long long int nHit = ge->second;
                genomeNetGDFoutput << g1 << "," << g2 << "," << nHit << endl;
	}
	// genome CC net
	// genome cc network gdf
        string gcNetGDF = outputDir + "/GenomeCCNetwork.gdf";
        ofstream gcNetGDFoutput(gcNetGDF.c_str());
        // genome id
        set<unsigned long long int> gcNetNodes;
        // cc net edges
        map<pair<unsigned long long int, unsigned long long int>,unsigned long long int> gcNetEdges;
	//
	map<unsigned long long int,ccInfo>::iterator c;
	for(c=CCs.begin(); c!=CCs.end(); ++c)
	{
		set<unsigned long long int> tmpGenomes = CCs[c->first].getGenomes();
		if(tmpGenomes.size() == 1)
		{
			gcNetNodes.insert(*tmpGenomes.begin());
		}
		else
		{
			while(!tmpGenomes.empty())
			{
				unsigned long long int g1 = *tmpGenomes.begin();
				tmpGenomes.erase(tmpGenomes.begin());
				set<unsigned long long int>::iterator i;
				for(i=tmpGenomes.begin(); i!=tmpGenomes.end(); ++i)
				{
					unsigned long long int g2 = *i;
					if(g1 != g2)
			                {
                			        gcNetNodes.insert(g1);
			                        gcNetNodes.insert(g2);
		                	        if(g1 < g2)
		            	        	{
	                		                gcNetEdges[make_pair(g1,g2)]+=1;
	                		        }
	                       			else
			                        {
	                		                gcNetEdges[make_pair(g2,g1)]+=1;
	                       			}
	                		}
	               	 		else
	                		{
	                        		gcNetNodes.insert(g1);
	                		}
				}
			}
		}

	}
	// Output
	gcNetGDFoutput << "nodedef>name VARCHAR,label VARCHAR" << endl;
        std::set<unsigned long long int>::iterator gn;
        for(gn=gcNetNodes.begin(); gn!=gcNetNodes.end(); ++gn)
        {
                gcNetGDFoutput << *gn << "," << genomes[*gn].getGenomeRealName() << endl;
        }
        gcNetGDFoutput << "edgedef>node1 VARCHAR,node2 VARCHAR,sharedcc DOUBLE" << endl;
        map<pair<unsigned long long int, unsigned long long int>,unsigned long long int>::iterator gne;
        for(gne=gcNetEdges.begin(); gne!=gcNetEdges.end(); ++gne)
        {
                //genome net nodes
                unsigned long long int g1 = gne->first.first;
                unsigned long long int g2 = gne->first.second;
                unsigned long long int scc = gne->second;
                gcNetGDFoutput << g1 << "," << g2 << "," << scc << endl;
        }
}
//end
