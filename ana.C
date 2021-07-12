#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include "RHAna.h"

/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=1)
{
  const char *hstfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events");
  //Declare an instance of our code class
  RHAna m_selec;
  
  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("VLLM100.root");
    //Set name of Output file
    hstfilename = "VLLM100_Modified.root";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  
  
  std::cout<<"Output file is "<<hstfilename<<endl;
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);
}

