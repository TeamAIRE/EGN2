/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        EGN2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/
#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

// Genes
class geneInfo
{
	public:
		geneInfo();
		geneInfo(unsigned int length, unsigned long long int ccID, std::string realName, unsigned long long int genomeID, bool visited);
		// Set or insert values
		void setLength(unsigned int);
		void setCcID(unsigned long long int);
		void setRealName(std::string);
		void setGenomeID(unsigned long long int);
		void setVisited(bool);
		void insertNeighbor(unsigned long long int);
		// Get values
		unsigned int getLength();
		unsigned long long int getCcID();
		std::string getRealName();
		unsigned long long int getGenomeID();
		bool isVisited();
		std::vector<unsigned long long int> getNeighbors();
	private:
		unsigned int m_length;
		unsigned long long int m_ccID;
		std::string m_realName;
		unsigned long long int m_genomeID;
		bool m_visited;
		std::vector<unsigned long long int> m_neighbors;
};

// EDGES	

class edgeValues
{
	public:
		edgeValues();
		edgeValues(unsigned short int qstart, unsigned short int qend, unsigned short int sstart, unsigned short int send, float pident, long double evalue);
		// Set values for the edge
		void setQstart(unsigned short int);
		void setQend(unsigned short int);
		void setSstart(unsigned short int);
		void setSend(unsigned short int);
		void setPident(float);
		void setEvalue(long double);
		// Get values 
		unsigned short int getQstart();
		unsigned short int getQend();
		unsigned short int getSstart();
		unsigned short int getSend();
		float getPident();
		long double getEvalue();
	private:
		// Query
		unsigned short int m_qstart;
		unsigned short int m_qend;
		// Subject
		unsigned short int m_sstart;
		unsigned short int m_send;
		// Hits values
		float m_pident;
		long double m_evalue;
};


// Connected components
class ccInfo
{
        public:
                ccInfo();
                ccInfo(unsigned long long int size);
                // Set or insert values
                void setSize(unsigned long long int);
                void insertGenome(unsigned long long int);
                // Get values
                unsigned long long int getSize();
                std::set<unsigned long long int> getGenomes();
		unsigned long long int getNbGenomes();
        private:
                unsigned long long int m_size;
                std::set<unsigned long long int> m_genomes;
};

// Genomes
class genomeInfo
{
        public:
                genomeInfo();
                genomeInfo(unsigned long long int size, std::string genomeRealName);
                // Set or insert values
                void setSize(unsigned long long int);
		void setGenomeRealName(std::string);
                // Get values
                unsigned long long int getSize();
		std::string getGenomeRealName();
		void addSeq();
        private:
                unsigned long long int m_size;
                std::string m_genomeRealName;
};



edgeValues getEdgeValues(std::map<std::pair<unsigned long long int,unsigned long long int>, edgeValues> &edges, unsigned long long int a, unsigned long long int b);

edgeValues reverseEdgeValues(edgeValues edge);



bool FoundIn(unsigned long long int value, std::vector<unsigned long long int>& List);

bool checkCoverage(unsigned long long int node1, unsigned long long int node2, unsigned short int minCov, std::map<unsigned long long int, geneInfo>& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, unsigned short int mutual);

float computeSD(float mean, std::vector<float>values);

#endif // FUNCTIONS_H_INCLUDED
