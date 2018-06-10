/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of FamilyDetector.
        
        MultiTwin is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/
#ifndef NETWORKS_H_INCLUDED
#define NETWORKS_H_INCLUDED

#include "functions.h"
#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>


void computeCC(std::map<unsigned long long int, geneInfo>& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, unsigned short int minCov, unsigned short int mutual, std::map<unsigned long long int, ccInfo>& CCs, std::map<unsigned long long int, genomeInfo>& genomes, std::string outputDir, std::set<unsigned long long int>& genomeList, unsigned short int verbose);

void computeGenomeAndCCnetwork(std::map<unsigned long long int, geneInfo>& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, std::map<unsigned long long int, ccInfo>& CCs, std::map<unsigned long long int, genomeInfo>& genomes, std::string outputDir, unsigned short int verbose);

#endif // NETWORKS_H_INCLUDED
