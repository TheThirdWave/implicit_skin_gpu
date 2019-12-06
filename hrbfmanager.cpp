#include "hrbfmanager.h"

HRBFManager::HRBFManager()
{
    numHRBFS = 0;
    recalc = true;
}

bool HRBFManager::getNeedRecalc()
{
    recalc = false;
    for(int i = 0; i < numHRBFS; i++)
    {
        if(hrbfs[i].getNeedRecalc() == true) recalc = true;
    }
    return recalc;
}


