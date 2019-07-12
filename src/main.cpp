#include <iostream>
#include<stdio.h>

#include "deQTConfig.h"
#include "rootgen.h"
using namespace std;

int main()
{
    rootGen* myWork=new rootGen();
    myWork->load_test=false;
    cout << "Console ROOT to michealGram generator." << endl;
    fprintf(stdout,"Version %d.%d\n",
               deQT_VERSION_MAJOR,
               deQT_VERSION_MINOR);
    cout<<"input .root filename should be file.root"<<endl;
    cout << "Started" << endl;
    myWork->processMacData();
    myWork->stirTemplateGen();
    cout << "stirTemplateGenOK" << endl;
    if(!myWork->load_test){

    myWork->createEmptyMichelogram();

    myWork->clearMichelogram();
    cout << "EmptyMichelogramGenerated" << endl;
    }
    cout << "loadROOTfiles" << endl;
    myWork->loadROOTfiles();

    myWork->createROOTMichelogram();
    myWork->showFirstaa();
    cout << "MichelogramGenerated" << endl;
    if(!myWork->load_test){
    //myWork->saveMichelogram();
    cout << "press any key to exit";
    cout << "MichelogramSaved" << endl;}
    delete myWork;
    std::cin.get();
    return 0;
}
