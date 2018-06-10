/*

        Written by Jananan PATHMANATHAN, 2017-2018
        
        This file is part of EGN2.
        
        EGN2 is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/
#ifndef RENAMEFASTA_H_INCLUDED
#define RENAMEFASTA_H_INCLUDED

#include <map>
#include <string>

void renameFasta(std::string fileIn, std::string outDir, std::map<std::string, unsigned long long int>& seqNewName, std::map<unsigned long long int, std::string>& seqRealName);
#endif // RENAMEFASTA_H_INCLUDED
