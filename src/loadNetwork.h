/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CleanBlast.
        
        MultiTwin is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef LOADNETWORK_H_INCLUDED
#define LOADNETWORK_H_INCLUDED

#include "functions.h"
#include <map>
#include <list>
#include <set>
#include <vector>
#include <string>
#include <sys/stat.h>


void loadNetwork(std::string fileIn, float pidentLimit, long double evalueLimit, std::map<std::string, unsigned short int>& positionsList, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, std::map<unsigned long long int, geneInfo>& genes, std::map<std::string, unsigned long long int>& seqNewName, std::map<unsigned long long int, std::string>& seqRealName, std::map<std::string, unsigned long long int>& gene2genome, std::set<unsigned long long int>& genomeList, unsigned short int verbose);

#endif // LOADNETWORK_H_INCLUDED
