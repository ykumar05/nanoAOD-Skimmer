#define RHAna_cxx

#include "RHAna.h"
#include <TH2.h>
#include <TStyle.h>
#include<fstream>
#include<iomanip>
using namespace std;

// The begin functions are called before the process function
// Things that should happen before looping over all events
// go here.
void RHAna::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   myTree = new TTree("EventTree","Modified_Tree");
   //Connecting the branches of the input root file to the newly defined variables.
   myTree ->Branch("nMuon",      &nMu_);
   myTree ->Branch("Muon_pt",   &muPt_);
   myTree ->Branch("Muon_eta", &muEta_);
}
void RHAna::SlaveBegin(TTree * /*tree*/)   
{
   TString option = GetOption();
   nEvtTotal = 0;
   nEvtRan = 0;
   //Create the Root file
   _HstFile = new TFile(_HstFileName,"recreate");
   cout<<"The file that is being genrated is: "<<_HstFileName<<endl;
}

void RHAna::SlaveTerminate()
{
  //cout<<"Printing the tree...."<<endl;
  //t->Print();
  cout<<"Writing the skimmed tree......"<<endl;
  myTree->Write();
  //Closing and writing the root file
  _HstFile->Write();
  _HstFile->Close();
  //Output to screen.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total events = "<<nEvtTotal<<endl;
  
 
 
}
void RHAna::Terminate()
{
  //   cout<<"Inside Terminate()"<<endl;
}



Bool_t RHAna::Process(Long64_t entry)
{

  // Decide to read in whole event or just some branches.
  // ReadLimited(0,..) reads in whole event
  // ReadLimited(1,..) reads in the branches defined in that function in header file.
  int readevent = ReadLimited(0,entry);
  if(readevent==0){ cout<<"Did not read in any branches.. quitting."<<endl; return kTRUE;}
  
  //Output processing information to screen based on verbosity level.
  if(_verbosity==0 && nEvtTotal%1000==0)cout<<"Processed "<<nEvtTotal<<" events..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" events..."<<endl;
 
  nEvtTotal++;


  // Some cleaning flags for running over different years.
  GoodEvt2018 = (_year==2018 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);

  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  //TTree *t = new TTree("EventTree","Modified_Tree");
  // If this event passes above filters, then we process it futher.
  if(GoodEvt){//&&nEvtRan<10001){
  

    nEvtRan++;

    

    /**************************************************************
     *                            MUONS                           *
     **************************************************************/

    //Clearing the arrays and initializing the integer variables.
    nMu_=0;
    muPt_.clear();
    muEta_.clear();
    

    
    //Defining the loop on the MUon collection
    for(unsigned int i=0; i< (nMuon); i++)
      {
	nMu_++;
	muPt_.push_back(Muon_pt[i]);
	muEta_.push_back(Muon_eta[i]);
	    
	//Filling the tree
	myTree->Fill();
    }
	  
  }

  
  return kTRUE;
}
