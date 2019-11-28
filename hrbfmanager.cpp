#include "hrbfmanager.h"

HRBFManager::HRBFManager()
{
    recalc = true;
}

bool HRBFManager::getNeedsRecalc()
{
    return recalc;
}


