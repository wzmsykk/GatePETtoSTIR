#ifndef UTCHAIN_H
#define UTCHAIN_H
#include <utype.h>
#include <string>
#include "TChain.h"
class uTChain{
public:
    uTChain();
    ~uTChain();
    int GetEntry(unsigned long int i);
    int SetBranchAddress(std::string name,void *p);

    bool branchAddressAllSet=false;
    Int_t *compton1, *compton2;
    Int_t *runID,
    *eventID1, *eventID2;
    Int_t *crystalID1, *crystalID2;
    Int_t *submoduleID1,
    *submoduleID2, *moduleID1, *moduleID2, *rsectorID1, *rsectorID2; Int_t
    *rotAngleID1, *rotAngleID2; Float_t *energy_1, *energy_2; Float_t *globalPosX1,
    *globalPosX2; Int_t *comptonPhantom1, *comptonPhantom2; Int_t *sourceID1,
    *sourceID2; Float_t *sourcePosX1, *sourcePosY1, *sourcePosZ1, *sourcePosX2,
    *sourcePosY2, *sourcePosZ2;

};
#endif // UTCHAIN_H
