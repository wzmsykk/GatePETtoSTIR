#include <iostream>
#include<stdio.h>

#include "deQTConfig.h"
#include "rootgen.h"
using namespace std;

int main()
{
    rootGen myWork=rootGen();
    cout << "Console ROOT to michealGram generator." << endl;
    fprintf(stdout,"Version %d.%d\n",
               deQT_VERSION_MAJOR,
               deQT_VERSION_MINOR);
    cout << "Started" << endl;
    myWork.stirTemplateGen();
    cout << "stirTemplateGenOK" << endl;
    myWork.createEmptyMichelogram();
    cout << "EmptyMichelogramGenerating" << endl;
    myWork.createROOTMichelogram();
    cout << "MichelogramGenerated" << endl;
    myWork.saveMichelogram();
    cout << "MichelogramSaved" << endl;
    std::cin.get();
    return 0;
}
