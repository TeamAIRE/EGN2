/*

        Written by Jananan PATHMANATHAN, 2017-2018

        This file is part of EGN2.

        ENG2 is shared under Creative commons licence:

        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)

        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef GETOPTIONS_H_INCLUDED
#define GETOPTIONS_H_INCLUDED

#include <string>


void getMainOptions(std::string configDir, long double& evalueLimit, float& pidentLimit, unsigned short int& minCov, unsigned short int& mutual, unsigned int& nCpu, unsigned short int& verbose);

#endif // GETOPTIONS_H_INCLUDED
