/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        EGN2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/
#include "functions.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

// Classe geneInfo
geneInfo::geneInfo() : m_length(0), m_ccID(0), m_realName(""), m_genomeID(0),m_visited(false)
{

}

geneInfo::geneInfo(unsigned int length, unsigned long long int ccID, string realName, unsigned long long int genomeID, bool visited) : m_length(length), m_ccID(ccID), m_realName(realName), m_genomeID(genomeID), m_visited(visited)
{

}
// Set or insert values
void geneInfo::setLength(unsigned int length)
{
        m_length = length;
}

void geneInfo::setCcID(unsigned long long int ccID)
{
        m_ccID = ccID;
}

void geneInfo::setRealName(string realName)
{
        m_realName = realName;
}

void geneInfo::setGenomeID(unsigned long long int genomeID)
{
        m_genomeID = genomeID;
}

void geneInfo::setVisited(bool visited)
{
        m_visited = visited;
}

void geneInfo::insertNeighbor(unsigned long long int neighbor)
{
	m_neighbors.push_back(neighbor);
}

// get values

unsigned int geneInfo::getLength()
{
	return m_length;
}

unsigned long long int geneInfo::getCcID()
{
        return m_ccID;
}

string geneInfo::getRealName()
{
        return m_realName;
}

unsigned long long int geneInfo::getGenomeID()
{
        return m_genomeID;
}

bool geneInfo::isVisited()
{
        return m_visited;
}

vector<unsigned long long int> geneInfo::getNeighbors()
{
	return m_neighbors;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Class edgeValues
edgeValues::edgeValues() : m_qstart(0), m_qend(0), m_sstart(0), m_send(0), m_pident(0), m_evalue(0)
{

}

edgeValues::edgeValues(unsigned short int qstart, unsigned short int qend, unsigned short int sstart, unsigned short int send, float pident, long double evalue) : m_qstart(qstart), m_qend(qend), m_sstart(sstart), m_send(send), m_pident(pident), m_evalue(evalue)
{

}

void edgeValues::setQstart(unsigned short int qstart)
{
	m_qstart = qstart;
}

void edgeValues::setQend(unsigned short int qend)
{
	m_qend = qend;
}

void edgeValues::setSstart(unsigned short int sstart)
{
	m_sstart = sstart;
}

void edgeValues::setSend(unsigned short int send)
{
	m_send = send;
}

void edgeValues::setPident(float pident)
{
        m_pident = pident;
}

void edgeValues::setEvalue(long double evalue)
{
	m_evalue = evalue;
}

unsigned short int edgeValues::getQstart()
{
	return m_qstart;
}

unsigned short int edgeValues::getQend()
{
	return m_qend;
}

unsigned short int edgeValues::getSstart()
{
	return m_sstart;
}

unsigned short int edgeValues::getSend()
{
	return m_send;
}

float edgeValues::getPident()
{
        return m_pident;
}

long double edgeValues::getEvalue()
{
	return m_evalue;
}


// CCinfo
ccInfo::ccInfo() : m_size(0)
{

}

ccInfo::ccInfo(unsigned long long int size) : m_size(size)
{

}

// set and insert info
void ccInfo::setSize(unsigned long long int size)
{
        m_size = size;
}

void ccInfo::insertGenome(unsigned long long int ID)
{
        m_genomes.insert(ID);
}

// get info
unsigned long long int ccInfo::getSize()
{
        return m_size;
}

set<unsigned long long int> ccInfo::getGenomes()
{
        return m_genomes;
}

unsigned long long int ccInfo::getNbGenomes()
{
        return m_genomes.size();
}

// genome info
genomeInfo::genomeInfo() : m_size(0), m_genomeRealName("")
{

}

genomeInfo::genomeInfo(unsigned long long int size, string genomeRealName) : m_size(size), m_genomeRealName(genomeRealName)
{

}

// set and insert info
void genomeInfo::setSize(unsigned long long int size)
{
        m_size = size;
}


void genomeInfo::setGenomeRealName(string genomeRealName)
{
        m_genomeRealName = genomeRealName;
}


void genomeInfo::addSeq()
{
        m_size+=1;
}

// get info
unsigned long long int genomeInfo::getSize()
{
        return m_size;
}
string genomeInfo::getGenomeRealName()
{
        return m_genomeRealName;
}


// Get edge values for a given hit
edgeValues getEdgeValues(map<pair<unsigned long long int,unsigned long long int>, edgeValues> &edges, unsigned long long int a, unsigned long long int b)
{
        //if ( edges.find(make_pair(a,b)) == edges.end() )
	if( a > b )
        {
                edgeValues toReverse = edges[make_pair(b,a)];
                return reverseEdgeValues(toReverse);
                //return reverseEdgeValues(edges[make_pair(b,a)]);
        }
        else
        {
                return edges[make_pair(a,b)];
        }
}


// Function reversing edge values
edgeValues reverseEdgeValues(edgeValues edge)
{
        unsigned short int qstart, qend, qlen, sstart, send, slen;
        // Initial values
        qstart = edge.getSstart();
        qend = edge.getSend();
        sstart = edge.getQstart();
        send = edge.getQend();
        // Set reverse values
        edge.setQstart(qstart);
        edge.setQend(qend);
        edge.setSstart(sstart);
        edge.setSend(send);
        // Return
        return edge;
}

// Dichotomic search to check if a value is present in a list
bool FoundIn(unsigned long long int value, vector<unsigned long long int>& List)
{
	bool found;
	unsigned long long int start;
	unsigned long long int end;
	unsigned long long int middle;
	// initialization
	found = false;
	start = 0;
	end = List.size();
	// start dichotomic search
	while(!found && ((end - start) > 1))
	{
		middle = (start + end)/2;
		found = (List.at(middle) == value);
		if(List.at(middle) > value)
		{
			end = middle;
		}
		else
		{
			start = middle;
		}
	}
	if(List.at(start) == value)
	{
		 return true;
	}
	else
	{
		return false;
	}
}


bool checkCoverage(unsigned long long int node1, unsigned long long int node2, unsigned short int minCov, map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, unsigned short int mutual)
{
        // Edge value
        edgeValues edge = getEdgeValues(edges, node1, node2);
        // Query coverage
        unsigned short int qcov = floor((((float)edge.getQend()-(float)edge.getQstart()+1.0)*100.0)/(float)genes[node1].getLength());
        // Subject coverage
        unsigned short int scov = floor((((float)edge.getSend()-(float)edge.getSstart()+1.0)*100.0)/(float)genes[node2].getLength());
        // Return true if both cov are higher or equal to the minimum coverage value 
        // fixed by the user or the default one (80%)
        if( (qcov >= minCov and scov >= minCov) and mutual == 1 )
        {
                return true;
        }
	else if((qcov >= minCov or scov >= minCov) and mutual == 0 )
	{
		return true;
	}
        else
        {
                return false;
        }
}


// Standard deviatio funtion
float computeSD(float mean, vector<float>values)
{
	float val = 0.0;
	unsigned long long int N = values.size();
	unsigned long long int i;
	for(i = 0; i < N; i++)
	{
		val += pow(values[i] - mean, 2);
	}
    	return sqrt(val / N);
}

//END
