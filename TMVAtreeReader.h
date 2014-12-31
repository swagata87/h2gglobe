//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 24 11:13:56 2014 by ROOT version 5.34/03
// from TTree TrainTree/TrainTree
// found on file: TMVA_photonid.root
//////////////////////////////////////////////////////////

#ifndef TMVAtreeReader_h
#define TMVAtreeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TMVAtreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           classID;
   Char_t          className[11];
   Float_t         r9_pho;
   Float_t         r9_cleaned_pho;
   Float_t         pT_pho; 
   Float_t         hoe_pho;
   Float_t         sieie_pho;
   Float_t         sieip_pho;
   Float_t         sieie_cleaned_pho;
   Float_t         E2byE5_pho;
   Float_t         E2byE5_cleaned_pho;
   Float_t         sieip_cleaned_pho;
   Float_t         scRaw_pho;
   Float_t         etawidth_pho;
   Float_t         phiwidth_pho;
   Float_t         sc_eta_pho;
   Float_t         itype;
   Float_t         genmatched_pho;
   Float_t         nPho;
   Float_t         nvtx;
   Float_t         Gradient_cat1;
   Float_t         Gradient_cat2;
   Float_t         weight;
   Float_t         Gradient;
   Float_t         trkiso03_pho;
   Float_t         sipip_pho;
   Float_t         pf_charged_iso_chosenvtx;
   Float_t         pf_charged_iso_vtx1;
   Float_t         pf_charged_iso_vtx2;
   Float_t         pf_charged_iso_badvtx;
   Float_t         pf_pho_iso;
   Float_t         rho;

   // List of branches
   TBranch        *b_classID;   //!
   TBranch        *b_className;   //!
   TBranch        *b_r9_cleaned_pho;   //!
   TBranch        *b_r9_pho;   //!
   TBranch        *b_pT_pho;   //!
   TBranch        *b_hoe_pho;   //!
   TBranch        *b_sieie_cleaned_pho;   //!
   TBranch        *b_sieip_cleaned_pho;   //!
   TBranch        *b_sieie_pho;   //!
   TBranch        *b_sieip_pho;   //!
   TBranch        *b_scRaw_pho;   //!
   TBranch        *b_etawidth_pho;   //!
   TBranch        *b_phiwidth_pho;   //!
   TBranch        *b_E2byE5_pho;
   TBranch        *b_E2byE5_cleaned_pho;
   TBranch        *b_sc_eta_pho;   //!
   TBranch        *b_itype;   //!
   TBranch        *b_genmatched_pho;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_Gradient_cat1;   //!
   TBranch        *b_Gradient_cat2;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_Gradient;   //!
   TBranch        *b_sipip_pho;   //!
   TBranch        *b_pf_charged_iso_chosenvtx;   //!
   TBranch        *b_pf_charged_iso_vtx1;
   TBranch        *b_pf_charged_iso_vtx2;
   TBranch        *b_pf_charged_iso_badvtx;   //!
   TBranch        *b_pf_pho_iso;   //!
   TBranch        *b_trkiso03_pho;   //!
   TBranch        *b_rho;   //!


   TMVAtreeReader(TTree *tree=0);
   virtual ~TMVAtreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TMVAtreeReader_cxx
TMVAtreeReader::TMVAtreeReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    /////    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TMVA_photonid_General.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TMVA_newTraining_onlyNewBDTsettings_noNewInputVar.root");
    if (!f || !f->IsOpen()) {
      //////  f = new TFile("TMVA_photonid_General.root");
      f = new TFile("TMVA_newTraining_onlyNewBDTsettings_noNewInputVar.root");
    }
    f->GetObject("TrainTree",tree);
    ////    f->GetObject("TestTree",tree);
    
  }
  Init(tree);
}

TMVAtreeReader::~TMVAtreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TMVAtreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TMVAtreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TMVAtreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("classID", &classID, &b_classID);
   fChain->SetBranchAddress("className", className, &b_className);
   fChain->SetBranchAddress("r9_cleaned_pho", &r9_cleaned_pho, &b_r9_cleaned_pho);
   fChain->SetBranchAddress("r9_pho", &r9_pho, &b_r9_pho);
   fChain->SetBranchAddress("hoe_pho", &hoe_pho, &b_hoe_pho);
   fChain->SetBranchAddress("pT_pho", &pT_pho, &b_pT_pho);
   fChain->SetBranchAddress("sieie_cleaned_pho", &sieie_cleaned_pho, &b_sieie_cleaned_pho);
   fChain->SetBranchAddress("sieip_cleaned_pho", &sieip_cleaned_pho, &b_sieip_cleaned_pho);
   fChain->SetBranchAddress("sieie_pho", &sieie_pho, &b_sieie_pho);
   fChain->SetBranchAddress("sieip_pho", &sieip_pho, &b_sieip_pho);
   fChain->SetBranchAddress("scRaw_pho", &scRaw_pho, &b_scRaw_pho);
   fChain->SetBranchAddress("etawidth_pho", &etawidth_pho, &b_etawidth_pho);
   fChain->SetBranchAddress("phiwidth_pho", &phiwidth_pho, &b_phiwidth_pho);
   fChain->SetBranchAddress("sipip_pho", &sipip_pho, &b_sipip_pho);
   fChain->SetBranchAddress("pf_charged_iso_chosenvtx", &pf_charged_iso_chosenvtx, &b_pf_charged_iso_chosenvtx);
   fChain->SetBranchAddress("pf_charged_iso_vtx1", &pf_charged_iso_vtx1, &b_pf_charged_iso_vtx1);
   fChain->SetBranchAddress("pf_charged_iso_vtx2", &pf_charged_iso_vtx2, &b_pf_charged_iso_vtx2);
   fChain->SetBranchAddress("pf_charged_iso_badvtx", &pf_charged_iso_badvtx, &b_pf_charged_iso_badvtx);
   fChain->SetBranchAddress("pf_pho_iso", &pf_pho_iso, &b_pf_pho_iso);
   fChain->SetBranchAddress("E2byE5_pho", &E2byE5_pho, &b_E2byE5_pho);
   fChain->SetBranchAddress("E2byE5_cleaned_pho", &E2byE5_cleaned_pho, &b_E2byE5_cleaned_pho);
   fChain->SetBranchAddress("trkiso03_pho", &trkiso03_pho, &b_trkiso03_pho);
   fChain->SetBranchAddress("sc_eta_pho", &sc_eta_pho, &b_sc_eta_pho);
   fChain->SetBranchAddress("itype", &itype, &b_itype);
   fChain->SetBranchAddress("genmatched_pho", &genmatched_pho, &b_genmatched_pho);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("Gradient_cat1", &Gradient_cat1, &b_Gradient_cat1);
   fChain->SetBranchAddress("Gradient_cat2", &Gradient_cat2, &b_Gradient_cat2);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("Gradient", &Gradient, &b_Gradient);
   fChain->SetBranchAddress("rho", &rho, &b_rho);

   Notify();
}

Bool_t TMVAtreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TMVAtreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TMVAtreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TMVAtreeReader_cxx
