//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 12 17:31:27 2015 by ROOT version 5.34/18
// from TTree event/Reduced tree
// found on file: /afs/cern.ch/work/s/sethzenz/GJ_v8/1_1.root
//////////////////////////////////////////////////////////

#ifndef RegressionTreeMaker2_h
#define RegressionTreeMaker2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class RegressionTreeMaker2 : public TSelector {
public :
   TTree *fChain;   //!pointer to the analyzed TTree or TChain
   TTree *outTree;
   TFile *outFile;

   int evt;
   float scRaw_pho;
   float r9_pho;
   float sc_eta_pho;
   float sc_phi_pho;
   float be5x5_raw;
   float hoe_pho;
   float sc_etawidth;
   float sc_phiwidth;

   float bc_sc_deta;
   float bc_sc_dphi;
   float seedOverRaw;
   float be3x3_eseed;
   float be5x5_eseed;
   float sieie_pho;
   float sieip_pho;
   float sipip_pho;

   float bemax_be;
   float be2nd_bemax;
   float betop_bemax;
   float bebottom_bemax;
   float beleft_bemax;
   float beright_bemax;
   float betopbottom_sym;
   float beleftright_sym;

   float bc2_sc_deta;
   float bc2_sc_dphi;
   float bc2OverRaw;
   float bc2e3x3_bc2e;
   float bc2e5x5_bc2e;
   float bc2_sieie;
   float bc2_sieip;
   float bc2_sipip;

   float bc2emax_bc2e;
   float bc2e2nd_bc2emax;
   float bc2etop_bc2emax;
   float bc2ebottom_bc2emax;
   float bc2eleft_bc2emax;
   float bc2eright_bc2emax;
   float bc2etopbottom_sym;
   float bc2eleftright_sym;

   float bclast_sc_deta;
   float bclast_sc_dphi;
   float bclastOverRaw;
   float bclaste3x3_bclaste;
   float bclaste5x5_bclaste;
   float bclast_sieie;
   float bclast_sieip;
   float bclast_sipip;

   float bclast2_sc_deta;
   float bclast2_sc_dphi;
   float bclast2OverRaw;
   float bclast2e3x3_bclast2e;
   float bclast2e5x5_bclast2e;
   float bclast2_sieie;
   float bclast2_sieip;
   float bclast2_sipip;

   int bieta_pho;
   int biphi_pho;
   float betacry_pho;
   float bphicry_pho;
   int bieta_var1;
   int biphi_var1;
   int bi_var2;
   int biphi_var2;

   int bc2ieta_pho;
   int bc2iphi_pho;
   float bc2etacry_pho;
   float bc2phicry_pho;
   int bc2ieta_var1;
   int bc2iphi_var1;
   int bc2i_var2;
   int bc2iphi_var2;

   int sc_nclu;
   float be3x3_be5x5;
   float bemax_be5x5;
   float be2nd_be5x5;
   float betop_be5x5;
   float bebottom_be5x5;
   float beleft_be5x5;
   float beright_be5x5;
   float be2x5max_be5x5;
   float be2x5top_be5x5;
   float be2x5bottom_be5x5;
   float be2x5left_be5x5;
   float be2x5right_be5x5;
 
   int nVtx;
   float rho;
   float genenergy_pho;
   int genmatched_pho;
   int isbarrel_pho;
   float tgtvar;

   #include "../branchdef/Limits.h"
   #include "../branchdef/branchdef.h"
   #include "../branchdef/treedef.h"
   
   Bool_t pho_genmatched[MAX_PHOTONS];
   Float_t pho_genenergy[MAX_PHOTONS];
   TBranch* b_pho_genmatched;
   TBranch* b_pho_genenergy;


 RegressionTreeMaker2(TTree * /*tree*/ =0) : fChain(0) { outTree = 0;}
   virtual ~RegressionTreeMaker2() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(RegressionTreeMaker2,0);
};

#endif

#ifdef RegressionTreeMaker2_cxx
void RegressionTreeMaker2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
  sc_xyz = new TClonesArray("TVector3");
  pho_p4 = new TClonesArray("TLorentzVector");
  bc_p4 = new TClonesArray("TLorentzVector");
  sc_p4 = new TClonesArray("TLorentzVector");
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("sc_xyz", &sc_xyz, &b_sc_xyz);
  fChain->SetBranchAddress("sc_n", &sc_n, &b_sc_n);
  fChain->SetBranchAddress("sc_sphi", sc_sphi, &b_sc_sphi);
  fChain->SetBranchAddress("sc_seta", sc_seta, &b_sc_seta);
  fChain->SetBranchAddress("sc_raw", sc_raw, &b_sc_raw);
  //fChain->SetBranchAddress("sc_pre", sc_pre, &b_sc_pre);
  fChain->SetBranchAddress("sc_bcseedind", sc_bcseedind, &b_sc_bcseedind);
  fChain->SetBranchAddress("sc_nclu", &sc_nclu);
  fChain->SetBranchAddress("bc_n", &bc_n, &b_bc_n);
  fChain->SetBranchAddress("bc_p4", &bc_p4, &b_bc_p4);
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  fChain->SetBranchAddress("pho_sieie", pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_sieip", pho_sieip, &b_pho_sieip);
  fChain->SetBranchAddress("pho_sipip", pho_sipip, &b_pho_sipip);
  //fChain->SetBranchAddress("pho_e1x5", pho_e1x5, &b_pho_e1x5);
  fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
  fChain->SetBranchAddress("pho_e5x5", pho_e5x5, &b_pho_e5x5);
  fChain->SetBranchAddress("pho_emaxxtal", pho_emaxxtal, &b_pho_emaxxtal);
  fChain->SetBranchAddress("pho_r9", pho_r9, &b_pho_r9);
  fChain->SetBranchAddress("pho_e2nd", pho_e2nd, &b_pho_e2nd);
  fChain->SetBranchAddress("pho_eright", pho_eright, &b_pho_eright);
  fChain->SetBranchAddress("pho_eleft", pho_eleft, &b_pho_eleft);
  fChain->SetBranchAddress("pho_etop", pho_etop, &b_pho_etop);
  fChain->SetBranchAddress("pho_ebottom", pho_ebottom, &b_pho_ebottom);
  fChain->SetBranchAddress("pho_etawidth", pho_etawidth, &b_pho_etawidth);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
  fChain->SetBranchAddress("sc_p4", &sc_p4, &b_sc_p4);
  fChain->SetBranchAddress("pho_scind", pho_scind, &b_pho_scind);
  fChain->SetBranchAddress("pho_hoe", pho_hoe_bc, &b_pho_hoe_bc);
  fChain->SetBranchAddress("vtx_std_n", &vtx_std_n, &b_vtx_std_n);
  fChain->SetBranchAddress("rho_algo1", &rho_algo1, &b_rho_algo1);
  fChain->SetBranchAddress("pho_bieta", pho_bieta, &b_pho_bieta);
  fChain->SetBranchAddress("pho_biphi", pho_biphi, &b_pho_biphi);
  fChain->SetBranchAddress("pho_betacry", pho_betacry, &b_pho_betacry);
  fChain->SetBranchAddress("pho_bphicry", pho_bphicry, &b_pho_bphicry);
  fChain->SetBranchAddress("pho_bc2ieta",   pho_bc2ieta,   &b_pho_bc2ieta);
  fChain->SetBranchAddress("pho_bc2iphi",   pho_bc2iphi,   &b_pho_bc2iphi);
  fChain->SetBranchAddress("pho_bc2etacry", pho_bc2etacry, &b_pho_bc2etacry);
  fChain->SetBranchAddress("pho_bc2phicry", pho_bc2phicry, &b_pho_bc2phicry);
  //fChain->SetBranchAddress("dipho_n", &dipho_n, &b_dipho_n);
  //fChain->SetBranchAddress("dipho_vtxind", &dipho_vtxind, &b_dipho_vtxind);
  fChain->SetBranchAddress("pho_genmatched", &pho_genmatched, &b_pho_genmatched);
  fChain->SetBranchAddress("pho_genenergy",  &pho_genenergy, &b_pho_genenergy);
  
  fChain->SetBranchAddress("pho_bc2e3x3",      pho_bc2e3x3,     &b_pho_bc2e3x3);
  fChain->SetBranchAddress("pho_bc2eta",       pho_bc2eta,      &b_pho_bc2eta);
  fChain->SetBranchAddress("pho_bc2phi",       pho_bc2phi,      &b_pho_bc2phi);
  fChain->SetBranchAddress("pho_bc2e",	       pho_bc2e,        &b_pho_bc2e);
  fChain->SetBranchAddress("pho_bc2e5x5",      pho_bc2e5x5,     &b_pho_bc2e5x5);
  fChain->SetBranchAddress("pho_bc2sieie",     pho_bc2sieie,    &b_pho_bc2sieie);
  fChain->SetBranchAddress("pho_bc2sieip",     pho_bc2sieip,    &b_pho_bc2sieip);
  fChain->SetBranchAddress("pho_bc2sipip",     pho_bc2sipip,    &b_pho_bc2sipip);
  fChain->SetBranchAddress("pho_bc2emax",      pho_bc2emax,     &b_pho_bc2emax);
  fChain->SetBranchAddress("pho_bc2e2nd",      pho_bc2e2nd,     &b_pho_bc2e2nd);
  fChain->SetBranchAddress("pho_bc2etop",      pho_bc2etop,     &b_pho_bc2etop);
  fChain->SetBranchAddress("pho_bc2ebottom",   pho_bc2ebottom,  &b_pho_bc2ebottom);
  fChain->SetBranchAddress("pho_bc2eleft",     pho_bc2eleft,    &b_pho_bc2eleft);
  fChain->SetBranchAddress("pho_bc2eright",    pho_bc2eright,   &b_pho_bc2eright);

  fChain->SetBranchAddress("pho_bclaste3x3",      pho_bclaste3x3,     &b_pho_bclaste3x3);
  fChain->SetBranchAddress("pho_bclasteta",       pho_bclasteta,      &b_pho_bclasteta);
  fChain->SetBranchAddress("pho_bclastphi",       pho_bclastphi,      &b_pho_bclastphi);
  fChain->SetBranchAddress("pho_bclaste",	  pho_bclaste,        &b_pho_bclaste);
  fChain->SetBranchAddress("pho_bclaste5x5",      pho_bclaste5x5,     &b_pho_bclaste5x5);
  fChain->SetBranchAddress("pho_bclastsieie",     pho_bclastsieie,    &b_pho_bclastsieie);
  fChain->SetBranchAddress("pho_bclastsieip",     pho_bclastsieip,    &b_pho_bclastsieip);
  fChain->SetBranchAddress("pho_bclastsipip",     pho_bclastsipip,    &b_pho_bclastsipip);

  fChain->SetBranchAddress("pho_bclast2e3x3",      pho_bclast2e3x3,     &b_pho_bclast2e3x3);
  fChain->SetBranchAddress("pho_bclast2eta",       pho_bclast2eta,      &b_pho_bclast2eta);
  fChain->SetBranchAddress("pho_bclast2phi",       pho_bclast2phi,      &b_pho_bclast2phi);
  fChain->SetBranchAddress("pho_bclast2e",	   pho_bclast2e,        &b_pho_bclast2e);
  fChain->SetBranchAddress("pho_bclast2e5x5",      pho_bclast2e5x5,     &b_pho_bclast2e5x5);
  fChain->SetBranchAddress("pho_bclast2sieie",     pho_bclast2sieie,    &b_pho_bclast2sieie);
  fChain->SetBranchAddress("pho_bclast2sieip",     pho_bclast2sieip,    &b_pho_bclast2sieip);
  fChain->SetBranchAddress("pho_bclast2sipip",     pho_bclast2sipip,    &b_pho_bclast2sipip);

  fChain->SetBranchAddress("pho_e2x5max",   &pho_e2x5max);
  fChain->SetBranchAddress("pho_e2x5top",   &pho_e2x5top);
  fChain->SetBranchAddress("pho_e2x5bottom",&pho_e2x5bottom);
  fChain->SetBranchAddress("pho_e2x5left",  &pho_e2x5left);
  fChain->SetBranchAddress("pho_e2x5right", &pho_e2x5right);
  
}

Bool_t RegressionTreeMaker2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef RegressionTreeMaker2_cxx
