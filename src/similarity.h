/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        ENG2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef SIMILARITY_H_INCLUDED
#define SIMILARITY_H_INCLUDED

#include <map>
#include <list>
#include <set>
#include <vector>
#include <string>
#include <sys/stat.h>


void runBlastp(std::string seqType, std::string subFile, std::string dbFile, std::string simToolConfig, unsigned int l);
void runDiamond(std::string seqType, std::string subFile, std::string dbFile, std::string simToolOptions, unsigned int l);
void runSimilarity( std::string fileIn, std::string simTool, std::string seqType, std::string configDir, std::string outDir, unsigned int nCpu, std::string timeInfo);

#endif // SIMILARITY_H_INCLUDED
