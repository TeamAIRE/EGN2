/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        EGN2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/
#ifndef ATTRIBUTES_H_INCLUDED
#define ATTRIBUTES_H_INCLUDED

#include <map>
#include <string>
#include "functions.h"

void getAttributes(std::string attrFile, std::map<unsigned long long int, unsigned long long int>& gene2genome, std::map<std::string, unsigned long long int>& seqNewName, std::map<unsigned long long int, genomeInfo>& genomes, std::map<std::string, unsigned long long int>& genomeNewName);

void getAttributes2(std::string attrFile, std::map<std::string, unsigned long long int>& gene2genome, std::map<unsigned long long int, genomeInfo>& genomes, std::map<std::string, unsigned long long int>& genomeNewName);

#endif // ATTRIBUTES_H_INCLUDED
