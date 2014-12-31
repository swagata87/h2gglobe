//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 21 06:03:31 2014 by ROOT version 5.34/18
// from TTree opttree_2Dwt/opttree_2Dwt
// found on file: 2DReweighted_opttree_8TeV.root
//////////////////////////////////////////////////////////

#ifndef ROCLegacy_h
#define ROCLegacy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class ROCLegacy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumis;
   Int_t           event;
   Int_t           itype;
   Int_t           nvtx;
   Int_t           nPU;
   Float_t         weight;
   Float_t         rho;
   Int_t           nPho;
   Float_t         pT_pho;
   Float_t         p4_eta_pho;
   Float_t         sc_eta_pho;
   Float_t         Et_pho;
   Float_t         E2byE5_pho;
   Float_t         E2byE5_cleaned_pho;
   Float_t         r9_pho;
   Float_t         r9_cleaned_pho;
   Float_t         sieie_pho;
   Float_t         sieie_cleaned_pho;
   Float_t         hoe_pho;
   Float_t         isEB_pho;
   Float_t         trkiso03_pho;
   Float_t         ecaliso03_pho;
   Float_t         hcaliso03_pho;
   Float_t         pf_pho_iso;
   Float_t         pfPhoIso_cleaned_t15;
   Float_t         pfPhoIso_cleaned_t10;
   Int_t           genmatched_pho;
   Float_t         etawidth_pho;
   Float_t         phiwidth_pho;
   Float_t         sieip_pho;
   Float_t         sieip_cleaned_pho;
   Float_t         scRaw_pho;
   Int_t           hasPixelSeed_pho;
   Float_t         pf_charged_iso_chosenvtx;
   Float_t         pf_charged_iso_badvtx;
   Float_t         sipip_pho;
   Float_t         mva_legacy;
   
   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_event;   //!
   TBranch        *b_itype;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_new_wt;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_pT_pho;   //!
   TBranch        *b_p4_eta_pho;   //!
   TBranch        *b_sc_eta_pho;   //!
   TBranch        *b_Et_pho;   //!
   TBranch        *b_E2byE5_pho;   //!
   TBranch        *b_E2byE5_cleaned_pho;   //!
   TBranch        *b_r9_pho;   //!
   TBranch        *b_r9_cleaned_pho;   //!
   TBranch        *b_sieie_pho;   //!
   TBranch        *b_sieie_cleaned_pho;   //!
   TBranch        *b_hoe_pho;   //!
   TBranch        *b_isEB_pho;   //!
   TBranch        *b_trkiso03_pho;   //!
   TBranch        *b_ecaliso03_pho;   //!
   TBranch        *b_hcaliso03_pho;   //!
   TBranch        *b_pf_pho_iso;   //!
   TBranch        *b_pfPhoIso_cleaned_t15;   //!
   TBranch        *b_pfPhoIso_cleaned_t10;   //!
   TBranch        *b_genmatched_pho;   //!
   TBranch        *b_etawidth_pho;   //!
   TBranch        *b_phiwidth_pho;   //!
   TBranch        *b_sieip_pho;   //!
   TBranch        *b_sieip_cleaned_pho;   //!
   TBranch        *b_scRaw_pho;   //!
   TBranch        *b_hasPixelSeed_pho;   //!
   TBranch        *b_pf_charged_iso_chosenvtx;   //!
   TBranch        *b_pf_charged_iso_badvtx;   //!
   TBranch        *b_sipip_pho;   //!
   TBranch        *b_mva_legacy; //!

   ROCLegacy(TTree *tree=0);
   virtual ~ROCLegacy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ROCLegacy_cxx
ROCLegacy::ROCLegacy(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("UCSDplotter/2DReweighted_opttree_8TeV_reproduceLegacy.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("UCSDplotter/2DReweighted_opttree_8TeV_reproduceLegacy.root");
      }
      f->GetObject("opttree_2Dwt",tree);

   }
   Init(tree);
}

ROCLegacy::~ROCLegacy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ROCLegacy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ROCLegacy::LoadTree(Long64_t entry)
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

void ROCLegacy::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("itype", &itype, &b_itype);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("weight", &weight, &b_new_wt);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("pT_pho", &pT_pho, &b_pT_pho);
   fChain->SetBranchAddress("p4_eta_pho", &p4_eta_pho, &b_p4_eta_pho);
   fChain->SetBranchAddress("sc_eta_pho", &sc_eta_pho, &b_sc_eta_pho);
   fChain->SetBranchAddress("Et_pho", &Et_pho, &b_Et_pho);
   fChain->SetBranchAddress("E2byE5_pho", &E2byE5_pho, &b_E2byE5_pho);
   fChain->SetBranchAddress("E2byE5_cleaned_pho", &E2byE5_cleaned_pho, &b_E2byE5_cleaned_pho);
   fChain->SetBranchAddress("r9_pho", &r9_pho, &b_r9_pho);
   fChain->SetBranchAddress("r9_cleaned_pho", &r9_cleaned_pho, &b_r9_cleaned_pho);
   fChain->SetBranchAddress("sieie_pho", &sieie_pho, &b_sieie_pho);
   fChain->SetBranchAddress("sieie_cleaned_pho", &sieie_cleaned_pho, &b_sieie_cleaned_pho);
   fChain->SetBranchAddress("hoe_pho", &hoe_pho, &b_hoe_pho);
   fChain->SetBranchAddress("isEB_pho", &isEB_pho, &b_isEB_pho);
   fChain->SetBranchAddress("trkiso03_pho", &trkiso03_pho, &b_trkiso03_pho);
   fChain->SetBranchAddress("ecaliso03_pho", &ecaliso03_pho, &b_ecaliso03_pho);
   fChain->SetBranchAddress("hcaliso03_pho", &hcaliso03_pho, &b_hcaliso03_pho);
   fChain->SetBranchAddress("pf_pho_iso", &pf_pho_iso, &b_pf_pho_iso);
   fChain->SetBranchAddress("pfPhoIso_cleaned_t15", &pfPhoIso_cleaned_t15, &b_pfPhoIso_cleaned_t15);
   fChain->SetBranchAddress("pfPhoIso_cleaned_t10", &pfPhoIso_cleaned_t10, &b_pfPhoIso_cleaned_t10);
   fChain->SetBranchAddress("genmatched_pho", &genmatched_pho, &b_genmatched_pho);
   fChain->SetBranchAddress("etawidth_pho", &etawidth_pho, &b_etawidth_pho);
   fChain->SetBranchAddress("phiwidth_pho", &phiwidth_pho, &b_phiwidth_pho);
   fChain->SetBranchAddress("sieip_pho", &sieip_pho, &b_sieip_pho);
   fChain->SetBranchAddress("sieip_cleaned_pho", &sieip_cleaned_pho, &b_sieip_cleaned_pho);
   fChain->SetBranchAddress("scRaw_pho", &scRaw_pho, &b_scRaw_pho);
   fChain->SetBranchAddress("hasPixelSeed_pho", &hasPixelSeed_pho, &b_hasPixelSeed_pho);
   fChain->SetBranchAddress("pf_charged_iso_chosenvtx", &pf_charged_iso_chosenvtx, &b_pf_charged_iso_chosenvtx);
   fChain->SetBranchAddress("pf_charged_iso_badvtx", &pf_charged_iso_badvtx, &b_pf_charged_iso_badvtx);
   fChain->SetBranchAddress("sipip_pho", &sipip_pho, &b_sipip_pho);
   fChain->SetBranchAddress("mva_legacy", &mva_legacy, &b_mva_legacy);
   Notify();
}

Bool_t ROCLegacy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ROCLegacy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ROCLegacy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ROCLegacy_cxx
