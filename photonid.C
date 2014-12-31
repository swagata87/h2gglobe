///// LATEST /////
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCuts.h"
//#endif

//////void photonid(char* outputFileName = "TMVA_photonid_General") {
void photonid(char* outputFileName = "TMVA_newTraining_onlyNewBDTsettings_noNewInputVar") {
  ///////void photonid(char* outputFileName = "TMVA_photonid_8TeVlegacy") {
  TMVA::Tools::Instance();
  
  char a[800];
  sprintf(a, "%s.root", outputFileName);
  TFile* outputFile = TFile::Open(a, "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory(outputFileName, outputFile, "!V:!Silent");
  
  ////// Variable ///// 
  /////  factory->AddVariable("sieie_cleaned_pho", 'F');
  factory->AddVariable("sieie_pho", 'F');
  
  //////factory->AddVariable("sieip_cleaned_pho", 'F');
  factory->AddVariable("sieip_pho", 'F');
  
  ///////factory->AddVariable("r9_cleaned_pho", 'F');
  factory->AddVariable("r9_pho", 'F');

  //////factory->AddVariable("E2byE5_cleaned_pho", 'F');
  factory->AddVariable("E2byE5_pho", 'F');

  factory->AddVariable("etawidth_pho", 'F');
  factory->AddVariable("phiwidth_pho", 'F');
  factory->AddVariable("pf_pho_iso", 'F');
  factory->AddVariable("pf_charged_iso_chosenvtx", 'F');
  // factory->AddVariable("pf_charged_iso_vtx1", 'F');
  // factory->AddVariable("pf_charged_iso_vtx2", 'F');
  factory->AddVariable("pf_charged_iso_badvtx", 'F');
  factory->AddVariable("rho", 'F');
  factory->AddVariable("sc_eta_pho", 'F');
  factory->AddVariable("scRaw_pho", 'F');

  ///// Spectator /////
  factory->AddSpectator("itype", 'F'); 
  factory->AddSpectator("genmatched_pho", 'I');
  factory->AddSpectator("nPho", 'I');
  factory->AddSpectator("nvtx",'I');
  factory->AddSpectator("pT_pho", 'F');
  /////  factory->AddSpectator("weight", 'F');
  TFile* input = TFile::Open("UCSDplotter/2DReweighted_opttree_8TeV_reproduceLegacy.root");
  
  TTree* inputTree = (TTree*)input->Get("opttree_2Dwt");
  outputFile->cd();
  
  factory->SetInputTrees(inputTree, "genmatched_pho==1 && (abs(sc_eta_pho))<1.479", "genmatched_pho==0 && (abs(sc_eta_pho)) < 1.479"); 
  factory->SetWeightExpression("weight");

  factory->PrepareTrainingAndTestTree("", "", "nTrain_Signal=50000:nTrain_Background=50000::nTest_Signal=50000:nTest_Background=50000:SplitMode=Random:!V" );

  ///////  TString vars = "sieie_cleaned_pho:sieip_cleaned_pho:r9_cleaned_pho:E2byE5_cleaned_pho:etawidth_pho:phiwidth_pho:pf_pho_iso:pf_charged_iso_chosenvtx:pf_charged_iso_vtx1:pf_charged_iso_vtx2:pf_charged_iso_badvtx:rho:sc_eta_pho:scRaw_pho";
  //TString vars = "sieie_pho:sieip_pho:r9_pho:E2byE5_pho:etawidth_pho:phiwidth_pho:pf_pho_iso:pf_charged_iso_chosenvtx:pf_charged_iso_vtx1:pf_charged_iso_vtx2:pf_charged_iso_badvtx:rho:sc_eta_pho:scRaw_pho";

  TString vars = "sieie_pho:sieip_pho:r9_pho:E2byE5_pho:etawidth_pho:phiwidth_pho:pf_pho_iso:pf_charged_iso_chosenvtx:pf_charged_iso_badvtx:rho:sc_eta_pho:scRaw_pho";   
  
  TMVA::MethodBase* fiCat2 = factory->BookMethod(TMVA::Types::kBDT, "Gradient", "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedGrad=F:nCuts=2000:MaxDepth=3:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=1000"); 
 
  //////TMVA::MethodBase* fiCat2 = factory->BookMethod(TMVA::Types::kBDT, "AdaBoost", "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=giniindex:nCuts=2000:PruneMethod=nopruning:MaxDepth=3");
 
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outputFile->Close();
    
  delete factory;
}
