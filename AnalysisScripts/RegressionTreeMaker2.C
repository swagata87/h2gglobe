#define RegressionTreeMaker2_cxx
// The class definition in RegressionTreeMaker2.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("RegressionTreeMaker2.C")
// Root > T->Process("RegressionTreeMaker2.C","some options")
// Root > T->Process("RegressionTreeMaker2.C+")
//

#include "RegressionTreeMaker2.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TMath.h>

float deltaPhi(float p1, float p2) {
  float dp = p1 - p2;
  if (dp > TMath::Pi()) 
    dp -= (2*TMath::Pi());

  return dp;
}

void RegressionTreeMaker2::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   std::cout << "Producing ntuple for: " << option.Data() << std::endl;
   TString dirname("/tmp/sani/regressionTree_");
   outFile = new TFile(dirname+option.Data()+TString(".root"), "RECREATE");

}

void RegressionTreeMaker2::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   outTree = new TTree("regressionTree","regressionTree");

   outTree->Branch("evt", &evt, "evt/I");
   outTree->Branch("scRaw_pho",&scRaw_pho,"scRaw_pho/F");
   outTree->Branch("r9_pho",&r9_pho,"r9_pho/F");
   outTree->Branch("sc_eta_pho",&sc_eta_pho,"sc_eta_pho/F");
   outTree->Branch("sc_phi_pho",&sc_phi_pho,"sc_phi_pho/F");
   outTree->Branch("be5x5_raw",&be5x5_raw,"be5x5_raw/F");
   outTree->Branch("hoe_pho",&hoe_pho,"_hoe_pho/F");
   outTree->Branch("sc_etawidth",&sc_etawidth,"sc_etawidth/F");
   outTree->Branch("sc_phiwidth",&sc_phiwidth,"sc_phiwidth/F");

   outTree->Branch("bc_sc_deta",&bc_sc_deta,"bc_sc_deta/F");
   outTree->Branch("bc_sc_dphi",&bc_sc_dphi,"bc_sc_dphi/F");
   outTree->Branch("eseed_raw",&seedOverRaw,"eseed_raw/F");
   outTree->Branch("be3x3_eseed",&be3x3_eseed,"be3x3_eseed/F");
   outTree->Branch("be5x5_eseed",&be5x5_eseed,"be5x5_eseed/F");
   outTree->Branch("sieie_pho",&sieie_pho,"sieie_pho/F");
   outTree->Branch("sieip_pho",&sieip_pho,"sieip_pho/F");
   outTree->Branch("sipip_pho",&sipip_pho,"sipip_pho/F");
   
   outTree->Branch("bemax_be",&bemax_be,"bemax_be/F");
   outTree->Branch("be2nd_bemax",&be2nd_bemax,"be2nd_bemax/F");
   outTree->Branch("betop_bemax",&betop_bemax,"betop_bemax/F");
   outTree->Branch("bebottom_bemax",&bebottom_bemax,"bebottom_bemax/F");
   outTree->Branch("beleft_bemax",&beleft_bemax,"beleft_bemax/F");
   outTree->Branch("beright_bemax",&beright_bemax,"beright_bemax/F");
   outTree->Branch("betopbottom_sym",&betopbottom_sym,"betopbottom_sym/F");
   outTree->Branch("beleftright_sym",&beleftright_sym,"beleftright_sym/F");

   outTree->Branch("bc2_sc_deta"  ,&bc2_sc_deta, "bc2_sc_deta/F");
   outTree->Branch("bc2_sc_dphi"  ,&bc2_sc_dphi, "bc2_sc_dphi/F");
   outTree->Branch("bc2e_raw"   ,&bc2OverRaw,  "bc2e_raw/F");
   outTree->Branch("bc2e3x3_bc2e" ,&bc2e3x3_bc2e,"bc2e3x3_bc2e/F");
   outTree->Branch("bc2e5x5_bc2e" ,&bc2e5x5_bc2e,"bc2e5x5_bc2e/F");
   outTree->Branch("bc2_sieie"    ,&bc2_sieie,   "bc2_sieie/F");
   outTree->Branch("bc2_sieip"    ,&bc2_sieip,   "bc2_sieip/F");
   outTree->Branch("bc2_sipip"    ,&bc2_sipip,   "bc2_sipip/F");
   
   outTree->Branch("bc2emax_bc2e",        &bc2emax_bc2e,        "bc2emax_bc2e/F");
   outTree->Branch("bc2e2nd_bc2emax",     &bc2e2nd_bc2emax,     "bc2e2nd_bc2emax/F");
   outTree->Branch("bc2etop_bc2emax",     &bc2etop_bc2emax,     "bc2etop_bc2emax/F");
   outTree->Branch("bc2ebottom_bc2emax",  &bc2ebottom_bc2emax,  "bc2ebottom_bc2emax/F");
   outTree->Branch("bc2eleft_bc2emax",    &bc2eleft_bc2emax,    "bc2eleft_bc2emax/F");
   outTree->Branch("bc2eright_bc2emax",   &bc2eright_bc2emax,   "bc2eright_bc2emax/F");
   outTree->Branch("bc2etopbottom_sym", &bc2etopbottom_sym, "bc2etopbottom_sym/F");
   outTree->Branch("bc2eleftright_sym", &bc2eleftright_sym, "bc2eleftright_sym/F");
   
   outTree->Branch("bclast_sc_deta"  ,&bclast_sc_deta, "bclast_sc_deta/F");
   outTree->Branch("bclast_sc_dphi"  ,&bclast_sc_dphi, "bclast_sc_dphi/F");
   outTree->Branch("bclaste_raw"   ,&bclastOverRaw,  "bclaste_raw/F");
   outTree->Branch("bclaste3x3_bclaste" ,&bclaste3x3_bclaste,"bclaste3x3_bclaste/F");
   outTree->Branch("bclaste5x5_bclaste" ,&bclaste5x5_bclaste,"bclaste5x5_bclaste/F");
   outTree->Branch("bclast_sieie"    ,&bclast_sieie,   "bclast_sieie/F");
   outTree->Branch("bclast_sieip"    ,&bclast_sieip,   "bclast_sieip/F");
   outTree->Branch("bclast_sipip"    ,&bclast_sipip,   "bclast_sipip/F");
   
   outTree->Branch("bclast2_sc_deta"  ,&bclast2_sc_deta, "bclast2_sc_deta/F");
   outTree->Branch("bclast2_sc_dphi"  ,&bclast2_sc_dphi, "bclast2_sc_dphi/F");
   outTree->Branch("bclast2e_raw"   ,&bclast2OverRaw,  "bclast2e_raw/F");
   outTree->Branch("bclast2e3x3_bclaste" ,&bclast2e3x3_bclast2e,"bclast2e3x3_bclast2e/F");
   outTree->Branch("bclast2e5x5_bclaste" ,&bclast2e5x5_bclast2e,"bclast2e5x5_bclast2e/F");
   outTree->Branch("bclast2_sieie"    ,&bclast2_sieie,   "bclast2_sieie/F");
   outTree->Branch("bclast2_sieip"    ,&bclast2_sieip,   "bclast2_sieip/F");
   outTree->Branch("bclast2_sipip"    ,&bclast2_sipip,   "bclast2_sipip/F");
   
   outTree->Branch("bieta_pho", &bieta_pho, "bieta_pho/I");
   outTree->Branch("biphi_pho", &biphi_pho, "biphi_pho/I");
   outTree->Branch("betacry_pho", &betacry_pho, "betacry_pho/F");
   outTree->Branch("bieta_var1", &bieta_var1, "bieta_var1/I");
   outTree->Branch("biphi_var1", &biphi_var1, "biphi_var1/I");
   outTree->Branch("bi_var2", &bi_var2, "bi_var2/I");
   outTree->Branch("biphi_var2", &biphi_var2, "biphi_var2/I");
   outTree->Branch("bphicry_pho", &bphicry_pho, "&bphicry_pho/F");

   outTree->Branch("bc2ieta_pho",   &bc2ieta_pho,   "bc2ieta_pho/I");
   outTree->Branch("bc2iphi_pho",   &bc2iphi_pho,   "bc2iphi_pho/I");
   outTree->Branch("bc2etacry_pho", &bc2etacry_pho, "bc2etacry_pho/F");
   outTree->Branch("bc2ieta_var1",  &bc2ieta_var1,  "bc2ieta_var1/I");
   outTree->Branch("bc2iphi_var1",  &bc2iphi_var1,  "bc2iphi_var1/I");
   outTree->Branch("bc2i_var2",     &bc2i_var2,     "bc2i_var2/I");
   outTree->Branch("bc2iphi_var2",  &bc2iphi_var2,  "bc2iphi_var2/I");
   outTree->Branch("bc2phicry_pho", &bc2phicry_pho, "bc2phicry_pho/F");

   outTree->Branch("sc_nclu", &sc_nclu, "sc_nclu/I");
   outTree->Branch("be3x3_be5x5", &be3x3_be5x5, "be3x3_be5x5/F");
   outTree->Branch("bemax_be5x5"	    , &bemax_be5x5	 , "bemax_be5x5/F");
   outTree->Branch("be2nd_be5x5"	    , &be2nd_be5x5	 , "be2nd_be5x5/F");
   outTree->Branch("betop_be5x5"	    , &betop_be5x5	 , "betop_be5x5/F");
   outTree->Branch("bebottom_be5x5"    , &bebottom_be5x5    , "bebottom_be5x5/F");
   outTree->Branch("beleft_be5x5"	    , &beleft_be5x5	 , "beleft_be5x5/F");
   outTree->Branch("beright_be5x5"	    , &beright_be5x5     , "beright_be5x5/F");
   outTree->Branch("be2x5max_be5x5"    , &be2x5max_be5x5    , "be2x5max_be5x5/F");
   outTree->Branch("be2x5top_be5x5"    , &be2x5top_be5x5    , "be2x5top_be5x5/F");
   outTree->Branch("be2x5bottom_be5x5" , &be2x5bottom_be5x5 , "be2x5bottom_be5x5/F");
   outTree->Branch("be2x5left_be5x5"   , &be2x5left_be5x5   , "be2x5left_be5x5/F");
   outTree->Branch("be2x5right_be5x5"  , &be2x5right_be5x5  , "be2x5right_be5x5/F");
   
   outTree->Branch("rho", &rho, "rho/F");
   outTree->Branch("vtx_n",&nVtx,"vtx_n/I");
   outTree->Branch("genenergy_pho",&genenergy_pho,"genenergy_pho/F");
   outTree->Branch("genmatched_pho",&genmatched_pho,"genmatched_pho/I");
   outTree->Branch("isbarrel_pho",&isbarrel_pho,"isbarrel_pho/I");

   outTree->Branch("tgtvar",&tgtvar,"tgtvar/F");
   //outTree->Branch("pre",&pre,"pre/F"); // endcp only 
 
   fOutput->Add(outTree);
}

Bool_t RegressionTreeMaker2::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either RegressionTreeMaker2::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  std::map<int,std::string> bad;
  //  bad[6614] = std::string("2_2.root");
  //  bad[6954] = std::string("2_6.root");

  if (bad.count(entry)) {
    const TFile *current_file = ((TChain*)(fChain))->GetCurrentFile();
    if (current_file) {
      const char *current_file_name = current_file->GetName();
      if (current_file_name) {
        if (strstr(current_file_name,bad[entry].c_str())) {
          return kFALSE;
        } else {
	  std::cout << " Bad entry number " << entry << " but not the right file (current_file_name=" << current_file_name <<"), continuing..." << std::endl;
        }
      } else {
	std::cout << " Could not get current_file_name" << std::endl;
      }
    } else {
      std::cout << " Could not get current_file" << std::endl;
    }
  }

  GetEntry(entry);
  for (int i = 0 ; i < pho_n ; i++) {
    TLorentzVector *phop4 = (TLorentzVector*)(pho_p4->At(i));
    if (!phop4) {
      //std::cout << "Photon " << i << "/" << pho_n << " seems to be null?" << std::endl;
    } else {

      if (pho_genmatched[i] && pho_genenergy[i] > 20. 
          && fabs(phop4->Eta()) < 3.0 && !(fabs(phop4->Eta()) > 1.4442 && fabs(phop4->Eta()) < 1.556)) {

	TLorentzVector *bcpos =(TLorentzVector*)bc_p4->At(sc_bcseedind[pho_scind[i]]);
	TVector3 *scpos = (TVector3*)sc_xyz->At(pho_scind[i]);

	evt = event;
        scRaw_pho = sc_raw[pho_scind[i]];
        r9_pho = pho_r9[i];
	sc_eta_pho = scpos->Eta();
        sc_phi_pho = scpos->Phi();
	be5x5_raw = pho_e5x5[i]/scRaw_pho;
	hoe_pho = pho_hoe_bc[i];
	sc_etawidth = sc_seta[pho_scind[i]];
	sc_phiwidth = sc_sphi[pho_scind[i]];

	bc_sc_deta = bcpos->Eta() - scpos->Eta();
	bc_sc_dphi = bcpos->Vect().DeltaPhi(*scpos);
	seedOverRaw = bcpos->E()/sc_raw[pho_scind[i]];
	be3x3_eseed = pho_e3x3[i]/bcpos->E();
	be5x5_eseed = pho_e5x5[i]/bcpos->E();
        sieie_pho = pho_sieie[i];
        sieip_pho = pho_sieip[i];
        sipip_pho = pho_sipip[i];

	bemax_be        = pho_emaxxtal[i]/bcpos->E();
	be2nd_bemax     = log(pho_e2nd[i]/pho_emaxxtal[i]);
	betop_bemax     = log(pho_etop[i]/pho_emaxxtal[i]);
	bebottom_bemax  = log(pho_ebottom[i]/pho_emaxxtal[i]);
	beleft_bemax    = log(pho_eleft[i]/pho_emaxxtal[i]);
	beright_bemax   = log(pho_eright[i]/pho_emaxxtal[i]);
	betopbottom_sym = (pho_etop[i]-pho_ebottom[i])/(pho_etop[i]+pho_ebottom[i]);
	beleftright_sym = (pho_eleft[i]-pho_eright[i])/(pho_eleft[i]+pho_eright[i]);	

	sc_nclu          = sc_nbc[pho_scind[i]]; 
	bemax_be5x5      = pho_emaxxtal[i]/pho_e5x5[i];
	be2nd_be5x5      = pho_e2nd[i]/pho_e5x5[i];
	betop_be5x5      = pho_etop[i]/pho_e5x5[i];
	bebottom_be5x5    = pho_ebottom[i]/pho_e5x5[i];
	beleft_be5x5      = pho_eleft[i]/pho_e5x5[i];
	beright_be5x5     = pho_eright[i]/pho_e5x5[i];
	be2x5max_be5x5    = pho_e2x5max[i]/pho_e5x5[i];
	be2x5top_be5x5    = pho_e2x5top[i]/pho_e5x5[i];
	be2x5bottom_be5x5 = pho_e2x5bottom[i]/pho_e5x5[i];
	be2x5left_be5x5   = pho_e2x5left[i]/pho_e5x5[i];
	be2x5right_be5x5  = pho_e2x5right[i]/pho_e5x5[i];
	
	bc2_sc_deta  = pho_bc2e3x3[i] ? (pho_bc2eta[i] - scpos->Eta()) : 0.;
	bc2_sc_dphi  = pho_bc2e3x3[i] ? deltaPhi(pho_bc2phi[i], scpos->Phi()) : 0.;
	bc2OverRaw   = pho_bc2e3x3[i] ? pho_bc2e[i]/sc_raw[pho_scind[i]] : 0;
	bc2e3x3_bc2e = pho_bc2e3x3[i] ? pho_bc2e3x3[i]/pho_bc2e[i] : 0.;
	bc2e5x5_bc2e = pho_bc2e3x3[i] ? pho_bc2e5x5[i]/pho_bc2e[i] : 0.;
	bc2_sieie    = pho_bc2e3x3[i] ? pho_bc2sieie[i] : 0.;
	bc2_sieip    = pho_bc2e3x3[i] ? pho_bc2sieip[i] : 0.;
	bc2_sipip    = pho_bc2e3x3[i] ? pho_bc2sipip[i] : 0.;

	bc2emax_bc2e       = pho_bc2e3x3[i] ? pho_bc2emax[i]/pho_bc2e[i] : 0.;
	bc2e2nd_bc2emax    = pho_bc2e3x3[i] ? log(pho_bc2e2nd[i]/pho_bc2emax[i]) : 0.;
	bc2etop_bc2emax    = pho_bc2e3x3[i] ? log(pho_bc2etop[i]/pho_bc2emax[i]) : 0.;
	bc2ebottom_bc2emax = pho_bc2e3x3[i] ? log(pho_bc2ebottom[i]/pho_bc2emax[i]) : 0.;
	bc2eleft_bc2emax   = pho_bc2e3x3[i] ? log(pho_bc2eleft[i]/pho_bc2emax[i]) : 0.;
	bc2eright_bc2emax  = pho_bc2e3x3[i] ? log(pho_bc2eright[i]/pho_bc2emax[i]) : 0.;
	bc2etopbottom_sym  = pho_bc2e3x3[i] ? (pho_bc2etop[i]-pho_bc2ebottom[i])/(pho_bc2etop[i]+pho_bc2ebottom[i]) : 0.;
	bc2eleftright_sym  = pho_bc2e3x3[i] ? (pho_bc2eleft[i]-pho_bc2eright[i])/(pho_bc2eleft[i]+pho_bc2eright[i]) : 0.;

	bclast_sc_deta  = pho_bclaste3x3[i] ? (pho_bclasteta[i] - scpos->Eta()) : 0.; 
	bclast_sc_dphi  = pho_bclaste3x3[i] ? deltaPhi(pho_bclastphi[i], scpos->Phi()) : 0.;
	bclastOverRaw   = pho_bclaste3x3[i] ? pho_bclaste[i]/sc_raw[pho_scind[i]] : 0;
	bclaste3x3_bclaste = pho_bclaste3x3[i] ? pho_bclaste3x3[i]/pho_bclaste[i] : 0.;
	bclaste5x5_bclaste = pho_bclaste3x3[i] ? pho_bclaste5x5[i]/pho_bclaste[i] : 0.;
	bclast_sieie    = pho_bclaste3x3[i] ? pho_bclastsieie[i] : 0.;
	bclast_sieip    = pho_bclaste3x3[i] ? pho_bclastsieip[i] : 0.;
	bclast_sipip    = pho_bclaste3x3[i] ? pho_bclastsipip[i] : 0.;

	bclast2_sc_deta  = pho_bclast2e3x3[i] ? (pho_bclast2eta[i] - scpos->Eta()) : 0.; 
	bclast2_sc_dphi  = pho_bclast2e3x3[i] ? deltaPhi(pho_bclast2phi[i], scpos->Phi()) : 0.;
	bclast2OverRaw   = pho_bclast2e3x3[i] ? pho_bclast2e[i]/sc_raw[pho_scind[i]] : 0;
	bclast2e3x3_bclast2e = pho_bclast2e3x3[i] ? pho_bclast2e3x3[i]/pho_bclast2e[i] : 0.;
	bclast2e5x5_bclast2e = pho_bclast2e3x3[i] ? pho_bclast2e5x5[i]/pho_bclast2e[i] : 0.;
	bclast2_sieie    = pho_bclast2e3x3[i] ? pho_bclast2sieie[i] : 0.;
	bclast2_sieip    = pho_bclast2e3x3[i] ? pho_bclast2sieip[i] : 0.;
	bclast2_sipip    = pho_bclast2e3x3[i] ? pho_bclast2sipip[i] : 0.;

	bieta_pho = pho_bieta[i];
	biphi_pho = pho_biphi[i];
	betacry_pho = pho_betacry[i];
	bphicry_pho = pho_bphicry[i];
	bieta_var1 = (pho_bieta[i]-1*abs(pho_bieta[i])/pho_bieta[i])%5;
	biphi_var1 = (pho_biphi[i]-1)%2;
	bi_var2 = (abs(pho_bieta[i])<=25)*((pho_bieta[i]-1*abs(pho_bieta[i])/pho_bieta[i])%25) + (abs(pho_bieta[i])>25)*((pho_bieta[i]-26*abs(pho_bieta[i])/pho_bieta[i])%20);
	biphi_var2 = (pho_biphi[i]-1)%20;

	bc2ieta_pho   = pho_bc2ieta[i];
	bc2iphi_pho   = pho_bc2iphi[i];
	bc2etacry_pho = pho_bc2etacry[i];
	bc2phicry_pho = pho_bc2phicry[i];

	bc2ieta_var1  = pho_bc2ieta[i] ? (pho_bc2ieta[i]-1*abs(pho_bc2ieta[i])/pho_bc2ieta[i])%5 : 0;
	bc2iphi_var1  = (pho_bc2iphi[i]-1)%2;
	bc2i_var2     = pho_bc2ieta[i] ? (abs(pho_bc2ieta[i])<=25)*((pho_bc2ieta[i]-1*abs(pho_bc2ieta[i])/pho_bc2ieta[i])%25) + (abs(pho_bc2ieta[i])>25)*((pho_bc2ieta[i]-26*abs(pho_bc2ieta[i])/pho_bc2ieta[i])%20) : 0;
	bc2iphi_var2  = (pho_bc2iphi[i]-1)%20;

	rho = rho_algo1;
	nVtx = vtx_std_n;
        genenergy_pho = pho_genenergy[i];
	genmatched_pho = int(pho_genmatched[i]);
	isbarrel_pho = (fabs(phop4->Eta()) < 1.4442);
        tgtvar = genenergy_pho/scRaw_pho;

        outTree->Fill();
      }
    }
  }

   return kTRUE;
}

void RegressionTreeMaker2::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void RegressionTreeMaker2::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  outTree = dynamic_cast<TTree*>(fOutput->FindObject("regressionTree"));
  //if (outTree) {
  outFile->cd();
  outTree->Write();
  outFile->Close();
  //}
    
}
