#define TMVAtreeReader_cxx
#include "TMVAtreeReader.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TGraph.h"

void TMVAtreeReader::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  //////  TFile* outputFile = new TFile("Overtraining_general.root","RECREATE");
  //////  TFile* outputFile = new TFile("Hgg_StudyTMVA_timeCl_2thEvt_1000trees.root","RECREATE");
  TFile* outputFile = new TFile("Hgg_General_8TeV_NewTrainingNewSettings_noNewVar.root","RECREATE");

  TH1F* h_weight_sig = new TH1F("weight_Sig", "h_weight_Sig", 200, 0.0, 20.0);
  TH1F* h_weight_bkg = new TH1F("weight_Bkg", "h_weight_Bkg", 200, 0.0, 20.0);

  TH1F* h_rho_sig_barr = new TH1F("rho_Sig_b", "h_rho_Sig_b", 200, 0.0, 200.0);
  TH1F* h_rho_bkg_barr = new TH1F("rho_Bkg_b", "h_rho_Bkg_b", 200, 0.0, 200.0);
  TH1F* h_rho_sig_endC = new TH1F("rho_Sig_e", "h_rho_Sig_e", 200, 0.0, 200.0);
  TH1F* h_rho_bkg_endC = new TH1F("rho_Bkg_e", "h_rho_Bkg_e", 200, 0.0, 200.0);

  TH1F* h_pt_sig_barr = new TH1F("pt_Sig_b", "h_pt_Sig_b", 400, 0.0, 400.0);
  TH1F* h_pt_bkg_barr = new TH1F("pt_Bkg_b", "h_pt_Bkg_b", 400, 0.0, 400.0);
  TH1F* h_pt_sig_endC = new TH1F("pt_Sig_e", "h_pt_Sig_e", 400, 0.0, 400.0);
  TH1F* h_pt_bkg_endC = new TH1F("pt_Bkg_e", "h_pt_Bkg_e", 400, 0.0, 400.0);

  TH1F* h_sc_eta_pho_sig_barr = new TH1F("sc_eta_pho_Sig_b", "h_sc_eta_pho_Sig_b", 150, 0.0, 2.5);
  TH1F* h_sc_eta_pho_bkg_barr = new TH1F("sc_eta_pho_Bkg_b", "h_sc_eta_pho_Bkg_b", 150, 0.0, 2.5);
  TH1F* h_sc_eta_pho_sig_endC = new TH1F("sc_eta_pho_Sig_e", "h_sc_eta_pho_Sig_e", 150, 0.0, 2.5);
  TH1F* h_sc_eta_pho_bkg_endC = new TH1F("sc_eta_pho_Bkg_e", "h_sc_eta_pho_Bkg_e", 150, 0.0, 2.5);

  TH1F* h_E2byE5_sig_barr = new TH1F("E2byE5_Sig_b", "h_E2byE5_Sig_b", 200, 0.0, 2.0);
  TH1F* h_E2byE5_bkg_barr = new TH1F("E2byE5_Bkg_b", "h_E2byE5_Bkg_b", 200, 0.0, 2.0);
  TH1F* h_E2byE5_sig_endC = new TH1F("E2byE5_Sig_e", "h_E2byE5_Sig_e", 200, 0.0, 2.0);
  TH1F* h_E2byE5_bkg_endC = new TH1F("E2byE5_Bkg_e", "h_E2byE5_Bkg_e", 200, 0.0, 2.0);

  TH1F* h_etawidth_sig_barr = new TH1F("etawidth_Sig_b", "h_etawidth_Sig_b", 300, 0.0, 0.03);
  TH1F* h_etawidth_bkg_barr = new TH1F("etawidth_Bkg_b", "h_etawidth_Bkg_b", 300, 0.0, 0.03);
  TH1F* h_etawidth_sig_endC = new TH1F("etawidth_Sig_e", "h_etawidth_Sig_e", 800, 0.0, 0.08);
  TH1F* h_etawidth_bkg_endC = new TH1F("etawidth_Bkg_e", "h_etawidth_Bkg_e", 800, 0.0, 0.08);

  TH1F* h_phiwidth_sig_barr = new TH1F("phiwidth_Sig_b", "h_phiwidth_Sig_b", 150, 0.0, 0.15);
  TH1F* h_phiwidth_bkg_barr = new TH1F("phiwidth_Bkg_b", "h_phiwidth_Bkg_b", 150, 0.0, 0.15);
  TH1F* h_phiwidth_sig_endC = new TH1F("phiwidth_Sig_e", "h_phiwidth_Sig_e", 150, 0.0, 0.15);
  TH1F* h_phiwidth_bkg_endC = new TH1F("phiwidth_Bkg_e", "h_phiwidth_Bkg_e", 150, 0.0, 0.15);
  
  TH1F* h_scRaw_sig_barr = new TH1F("scRaw_Sig_b", "h_scRaw_Sig_b", 500, 0.0, 500.0);
  TH1F* h_scRaw_bkg_barr = new TH1F("scRaw_Bkg_b", "h_scRaw_Bkg_b", 500, 0.0, 500.0);
  TH1F* h_scRaw_sig_endC = new TH1F("scRaw_Sig_e", "h_scRaw_Sig_e", 1500, 0.0, 1500.0);
  TH1F* h_scRaw_bkg_endC = new TH1F("scRaw_Bkg_e", "h_scRaw_Bkg_e", 1500, 0.0, 1500.0);

  TH1F* h_sieip_sig_barr = new TH1F("sieip_pho_Sig_b", "h_sieip_pho_Sig_b", 200, 0.0, 0.001);
  TH1F* h_sieip_bkg_barr = new TH1F("sieip_pho_Bkg_b", "h_sieip_pho_Bkg_b", 200, 0.0, 0.001);
  TH1F* h_sieip_sig_endC = new TH1F("sieip_pho_Sig_e", "h_sieip_pho_Sig_e", 200, 0.0, 0.004);
  TH1F* h_sieip_bkg_endC = new TH1F("sieip_pho_Bkg_e", "h_sieip_pho_Bkg_e", 200, 0.0, 0.004);

  TH1F* h_pf_pho_iso_sig_barr = new TH1F("pf_pho_iso_Sig_b", "h_pf_pho_iso_Sig_b", 400, 0.0, 400.0);
  TH1F* h_pf_pho_iso_bkg_barr = new TH1F("pf_pho_iso_Bkg_b", "h_pf_pho_iso_Bkg_b", 400, 0.0, 400.0);
  TH1F* h_pf_pho_iso_sig_endC = new TH1F("pf_pho_iso_Sig_e", "h_pf_pho_iso_Sig_e", 400, 0.0, 400.0);
  TH1F* h_pf_pho_iso_bkg_endC = new TH1F("pf_pho_iso_Bkg_e", "h_pf_pho_iso_Bkg_e", 400, 0.0, 400.0);

  TH1F* h_pf_pho_iso_by_pt_sig_barr = new TH1F("pf_pho_iso_Pt_Sig_b", "h_pf_pho_iso_Pt_Sig_b", 400, 0.0, 10.0);
  TH1F* h_pf_pho_iso_by_pt_bkg_barr = new TH1F("pf_pho_iso_Pt_Bkg_b", "h_pf_pho_iso_Pt_Bkg_b", 400, 0.0, 10.0);
  TH1F* h_pf_pho_iso_by_pt_sig_endC = new TH1F("pf_pho_iso_Pt_Sig_e", "h_pf_pho_iso_Pt_Sig_e", 400, 0.0, 10.0);
  TH1F* h_pf_pho_iso_by_pt_bkg_endC = new TH1F("pf_pho_iso_Pt_Bkg_e", "h_pf_pho_iso_Pt_Bkg_e", 400, 0.0, 10.0);

  TH1F* h_r9_sig_barr = new TH1F("r9_Sig_b", "h_r9_Sig_b", 200, 0.0, 2.0);
  TH1F* h_r9_bkg_barr = new TH1F("r9_Bkg_b", "h_r9_Bkg_b", 200, 0.0, 2.0);
  TH1F* h_r9_sig_endC = new TH1F("r9_Sig_e", "h_r9_Sig_e", 200, 0.0, 2.0);
  TH1F* h_r9_bkg_endC = new TH1F("r9_Bkg_e", "h_r9_Bkg_e", 200, 0.0, 2.0);

  TH1F* h_sieie_sig_barr = new TH1F("sieie_Sig_b", "h_sieie_Sig_b", 500, 0.0, 0.05);
  TH1F* h_sieie_bkg_barr = new TH1F("sieie_Bkg_b", "h_sieie_Bkg_b", 500, 0.0, 0.05);
  TH1F* h_sieie_sig_endC = new TH1F("sieie_Sig_e", "h_sieie_Sig_e", 500, 0.0, 0.1);
  TH1F* h_sieie_bkg_endC = new TH1F("sieie_Bkg_e", "h_sieie_Bkg_e", 500, 0.0, 0.1);

  TH1F* h_chiso_chosenvtx_sig_barr = new TH1F("pf_charged_iso_chosenvtx_Sig_b", "h_pf_charged_iso_chosenvtx_Sig_b", 500, 0.0, 20.0);
  TH1F* h_chiso_chosenvtx_bkg_barr = new TH1F("pf_charged_iso_chosenvtx_Bkg_b", "h_pf_charged_iso_chosenvtx_Bkg_b", 500, 0.0, 20.0);
  TH1F* h_chiso_chosenvtx_sig_endC = new TH1F("pf_charged_iso_chosenvtx_Sig_e", "h_pf_charged_iso_chosenvtx_Sig_e", 500, 0.0, 20.0);
  TH1F* h_chiso_chosenvtx_bkg_endC = new TH1F("pf_charged_iso_chosenvtx_Bkg_e", "h_pf_charged_iso_chosenvtx_Bkg_e", 500, 0.0, 20.0);

  TH1F* h_chiso_vtx1_sig_barr = new TH1F("pf_charged_iso_vtx1_Sig_b", "h_pf_charged_iso_vtx1_Sig_b", 500, 0.0, 20.0);
  TH1F* h_chiso_vtx1_bkg_barr = new TH1F("pf_charged_iso_vtx1_Bkg_b", "h_pf_charged_iso_vtx1_Bkg_b", 500, 0.0, 20.0);

  TH1F* h_chiso_vtx2_sig_barr = new TH1F("pf_charged_iso_vtx2_Sig_b", "h_pf_charged_iso_vtx2_Sig_b", 500, 0.0, 20.0);
  TH1F* h_chiso_vtx2_bkg_barr = new TH1F("pf_charged_iso_vtx2_Bkg_b", "h_pf_charged_iso_vtx2_Bkg_b", 500, 0.0, 20.0);

  TH1F* h_pf_charged_iso_badvtx_sig_barr = new TH1F("pf_charged_iso_badvtx_Sig_b", "h_pf_charged_iso_badvtx_Sig_b", 500, 0.0, 50.0);
  TH1F* h_pf_charged_iso_badvtx_bkg_barr = new TH1F("pf_charged_iso_badvtx_Bkg_b", "h_pf_charged_iso_badvtx_Bkg_b", 500, 0.0, 50.0);
  TH1F* h_pf_charged_iso_badvtx_sig_endC = new TH1F("pf_charged_iso_badvtx_Sig_e", "h_pf_charged_iso_badvtx_Sig_e", 500, 0.0, 50.0);
  TH1F* h_pf_charged_iso_badvtx_bkg_endC = new TH1F("pf_charged_iso_badvtx_Bkg_e", "h_pf_charged_iso_badvtx_Bkg_e", 500, 0.0, 50.0);
 
  ////// ------ ///////
  TH1F* h_BDToutput_sig_barr = new TH1F("BDToutput_Sig_b", "h_BDToutput_Sig_b", 300, -1.5, 1.5);
  TH1F* h_BDToutput_bkg_barr = new TH1F("BDToutput_Bkg_b", "h_BDToutput_Bkg_b", 300, -1.5, 1.5);
  TH1F* h_BDToutput_sig_endC = new TH1F("BDToutput_Sig_e", "h_BDToutput_Sig_e", 300, -1.5, 1.5);
  TH1F* h_BDToutput_bkg_endC = new TH1F("BDToutput_Bkg_e", "h_BDToutput_Bkg_e", 300, -1.5, 1.5);

  TH1F*  SigEff_vs_eta_deno_barr = new TH1F("Sig_eff_deno_vs_eta_barr", "h_Sig_deno_barr", 40, -2.0, 2.0);
  TH1F*  SigEff_vs_eta_num_barr  = new TH1F("Sig_eff_num_vs_eta_barr",  "h_Sig_num_barr",  40, -2.0, 2.0);
  TH1F*  BkgEff_vs_eta_deno_barr = new TH1F("Bkg_eff_deno_vs_eta_barr", "h_Bkg_deno_barr", 40, -2.0, 2.0);
  TH1F*  BkgEff_vs_eta_num_barr  = new TH1F("Bkg_eff_num_vs_eta_barr",  "h_Bkg_num_barr",  40, -2.0, 2.0);

  TH1F*  SigEff_vs_eta_deno_endC = new TH1F("Sig_eff_deno_vs_eta_endC", "h_Sig_deno_endC", 60, -3.0, 3.0);
  TH1F*  SigEff_vs_eta_num_endC  = new TH1F("Sig_eff_num_vs_eta_endC",  "h_Sig_num_endC",  60, -3.0, 3.0);
  TH1F*  BkgEff_vs_eta_deno_endC = new TH1F("Bkg_eff_deno_vs_eta_endC", "h_Bkg_deno_endC", 60, -3.0, 3.0);
  TH1F*  BkgEff_vs_eta_num_endC  = new TH1F("Bkg_eff_num_vs_eta_endC",  "h_Bkg_num_endC",  60, -3.0, 3.0);

  TH1F*  SigEff_vs_pt_deno_barr = new TH1F("Sig_eff_deno_vs_pt_barr", "h1_Sig_deno_barr", 100, 0.0, 400.0);
  TH1F*  SigEff_vs_pt_num_barr  = new TH1F("Sig_eff_num_vs_pt_barr",  "h1_Sig_num_barr",  100, 0.0, 400.0);
  TH1F*  BkgEff_vs_pt_deno_barr = new TH1F("Bkg_eff_deno_vs_pt_barr", "h1_Bkg_deno_barr", 100, 0.0, 400.0);
  TH1F*  BkgEff_vs_pt_num_barr  = new TH1F("Bkg_eff_num_vs_pt_barr",  "h1_Bkg_num_barr",  100, 0.0, 400.0);
 
  TH1F*  SigEff_vs_pt_deno_endC = new TH1F("Sig_eff_deno_vs_pt_endC", "h1_Sig_deno_endC", 100, 0.0, 400.0);
  TH1F*  SigEff_vs_pt_num_endC  = new TH1F("Sig_eff_num_vs_pt_endC",  "h1_Sig_num_endC",  100, 0.0, 400.0);
  TH1F*  BkgEff_vs_pt_deno_endC = new TH1F("Bkg_eff_deno_vs_pt_endC", "h1_Bkg_deno_endC", 100, 0.0, 400.0);
  TH1F*  BkgEff_vs_pt_num_endC  = new TH1F("Bkg_eff_num_vs_pt_endC",  "h1_Bkg_num_endC",  100, 0.0, 400.0);

  ///// 2D /////
  TH2D*  SigEff_etaPt_deno_barr = new TH2D("Sig_eff_deno_vs_etaPt_barr", "h2D_Sig_deno_barr", 100, -2.0, 3.0, 100, 0.0, 400.0);
  TH2D*  SigEff_etaPt_num_barr  = new TH2D("Sig_eff_num_vs_etaPt_barr",  "h2D_Sig_num_barr",  100, -2.0, 3.0, 100, 0.0, 400.0);
  TH2D*  BkgEff_etaPt_deno_barr = new TH2D("Bkg_eff_deno_vs_etaPt_barr", "h2D_Bkg_deno_barr", 100, -2.0, 3.0, 100, 0.0, 400.0);
  TH2D*  BkgEff_etaPt_num_barr  = new TH2D("Bkg_eff_num_vs_etaPt_barr",  "h2D_Bkg_num_barr",  100, -2.0, 3.0, 100, 0.0, 400.0);
  /////// ** ///////

  TH1F* R9_siglikeBkg_barr       = new TH1F("R9_SigLikeBkg_barrel",       "SigLikeBkg_R9_barrel",       200, 0.0, 2.0);
  TH1F* BDT_siglikeBkg_barr      = new TH1F("BDT_SigLikeBkg_barrel",      "SigLikeBkg_BDT_barrel",      300, 0.8, 1.1);
  TH1F* scRaw_siglikeBkg_barr    = new TH1F("scRaw_SigLikeBkg_barrel",    "SigLikeBkg_scRaw_barrel",    500, 0.0, 500.0);
  TH1F* scEta_siglikeBkg_barr    = new TH1F("scEta_SigLikeBkg_barrel",    "SigLikeBkg_scEta_barrel",    250, 0.0, 2.5);
  TH1F* sieie_siglikeBkg_barr    = new TH1F("sieie_SigLikeBkg_barrel",    "SigLikeBkg_sieie_barrel",    500, 0.0, 0.05);
  TH1F* sieip_siglikeBkg_barr    = new TH1F("sieip_SigLikeBkg_barrel",    "SigLikeBkg_sieip_barrel",    200, 0.0, 0.001);
  TH1F* etawidth_siglikeBkg_barr = new TH1F("etawidth_SigLikeBkg_barrel", "SigLikeBkg_etawidth_barrel", 300, 0.0, 0.03);
  TH1F* phiwidth_siglikeBkg_barr = new TH1F("phiwidth_SigLikeBkg_barrel", "SigLikeBkg_phiwidth_barrel", 150, 0.0, 0.15);
  TH1F* ChIsoChosen_siglikeBkg_barr = new TH1F("chisochosen_SigLikeBkg_barrel", "SigLikeBkg_chisochosen_barrel", 500, 0.0, 20.0);

  
  SigEff_vs_eta_deno_barr->Sumw2();
  SigEff_vs_pt_deno_barr->Sumw2();
  BkgEff_vs_eta_deno_barr->Sumw2();
  SigEff_vs_eta_num_barr->Sumw2();
  SigEff_vs_pt_num_barr->Sumw2();
  BkgEff_vs_eta_num_barr->Sumw2();
  BkgEff_vs_pt_deno_barr->Sumw2();
  BkgEff_vs_pt_num_barr->Sumw2();
  SigEff_vs_eta_deno_endC->Sumw2();
  SigEff_vs_eta_num_endC->Sumw2();
  SigEff_vs_pt_deno_endC->Sumw2();
  SigEff_vs_pt_num_endC->Sumw2();   
  BkgEff_vs_eta_deno_endC->Sumw2();
  BkgEff_vs_eta_num_endC->Sumw2();
  BkgEff_vs_pt_deno_endC->Sumw2();
  BkgEff_vs_pt_num_endC->Sumw2();
     
  int counter_barr_sig_n[41]; //=0;
  int counter_barr_bkg_n[41]; //=0;
  int counter_barr_sig_d[41]; //=0;
  int counter_barr_bkg_d[41]; //=0;

  int counter_endc_sig_n[41]; //=0;
  int counter_endc_bkg_n[41]; //=0;
  int counter_endc_sig_d[41]; //=0;
  int counter_endc_bkg_d[41]; //=0;

  for (int j=0; j<41; j++) {
    counter_barr_sig_n[j]=0;
    counter_barr_bkg_n[j]=0;
    counter_barr_sig_d[j]=0;
    counter_barr_bkg_d[j]=0;

    counter_endc_sig_n[j]=0;
    counter_endc_bkg_n[j]=0;
    counter_endc_sig_d[j]=0;
    counter_endc_bkg_d[j]=0;
  }

  int count_siglike_bkg = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //  for (Long64_t jentry=0; jentry<15; jentry++) {
    // std::cout << "\n Event " << jentry << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    ////////---------------------------------------///////
    ////////-- Investigating about signal-like bkg --////// 
    ///////----------------------------------------///////

    if ( (fabs(sc_eta_pho)) < 1.479) {                   /// BARREL ///
      if ( (genmatched_pho==0) && (Gradient>0.9)) {     /// BKG with high BDT score ///  
	//	std::cout << "Evt no. " << jentry << "  Gradient=" << Gradient << "   r9_cleaned=" << r9_cleaned_pho << std::endl;
	/*
	//// call Sumw2() ////
	R9_siglikeBkg_barr->Sumw2();
        BDT_siglikeBkg_barr->Sumw2();
	scRaw_siglikeBkg_barr->Sumw2();
	scEta_siglikeBkg_barr->Sumw2();
	sieie_siglikeBkg_barr->Sumw2();
	sieip_siglikeBkg_barr->Sumw2();
	etawidth_siglikeBkg_barr->Sumw2();
	phiwidth_siglikeBkg_barr->Sumw2();
	ChIsoChosen_siglikeBkg_barr->Sumw2();
	*/

	R9_siglikeBkg_barr->Fill(r9_pho,weight);
	////R9_siglikeBkg_barr->Fill(r9_cleaned_pho,weight);
	BDT_siglikeBkg_barr->Fill(Gradient,weight);
	scRaw_siglikeBkg_barr->Fill(scRaw_pho,weight);
	scEta_siglikeBkg_barr->Fill((fabs(sc_eta_pho)),weight);
	sieie_siglikeBkg_barr->Fill(sieie_pho,weight);
        sieip_siglikeBkg_barr->Fill(sieip_pho,weight);
	////sieie_siglikeBkg_barr->Fill(sieie_cleaned_pho,weight);
	////sieip_siglikeBkg_barr->Fill(sieip_cleaned_pho,weight);
	etawidth_siglikeBkg_barr->Fill(etawidth_pho,weight);
	phiwidth_siglikeBkg_barr->Fill(phiwidth_pho,weight);
	ChIsoChosen_siglikeBkg_barr->Fill(pf_charged_iso_chosenvtx,weight);
	count_siglike_bkg++;
      }
    }
    /// ROC ///
    if ( (fabs(sc_eta_pho)) < 1.479) { /// BARREL ///
      float x=-1.00;
      for (int j=0; j<41; j++) {
	// std::cout << "\n cut value = " << x << " Gradient = " << Gradient << std::endl;
	if ((Gradient>x || Gradient==x) && genmatched_pho==1.0) counter_barr_sig_n[j]++;
	if ((Gradient>x || Gradient==x) && genmatched_pho==0.0) counter_barr_bkg_n[j]++;
	if (genmatched_pho==1.0) counter_barr_sig_d[j]++;
	if (genmatched_pho==0.0) counter_barr_bkg_d[j]++;
	//    std::cout << counter_barr_sig_n[j] << " " << counter_barr_bkg_n[j] << std::endl;
	x=x+0.05;
      }
    }
    /*
    if ( (fabs(sc_eta_pho)) > 1.479) { /// ENDCAP ///
      float x=-1.00;
      for (int j=0; j<41; j++) {
        // std::cout << "\n cut value = " << x << " Gradient = " << Gradient << std::endl;
        if ((Gradient>x || Gradient==x) && genmatched_pho==1.0) counter_endc_sig_n[j]++;
        if ((Gradient>x || Gradient==x) && genmatched_pho==0.0) counter_endc_bkg_n[j]++;
        if (genmatched_pho==1.0) counter_endc_sig_d[j]++;
        if (genmatched_pho==0.0) counter_endc_bkg_d[j]++;
        //    std::cout << counter_endc_sig_n[j] << " " << counter_endc_bkg_n[j] << std::endl;
        x=x+0.05;
      }
    }
    */

    //// Eff ////
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// signal barrel ///
      SigEff_vs_eta_deno_barr->Fill(sc_eta_pho,weight);
      SigEff_etaPt_deno_barr->Fill(sc_eta_pho,pT_pho,weight);
      if (Gradient>0.0)  {
	SigEff_vs_eta_num_barr->Fill(sc_eta_pho,weight);
	SigEff_etaPt_num_barr->Fill(sc_eta_pho,pT_pho,weight);
      //  h_weight_sig->Fill(weight);
      }
    }
    
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// signal in barrel ///
      SigEff_vs_pt_deno_barr->Fill(pT_pho,weight);
      if (Gradient>0.0)  SigEff_vs_pt_num_barr->Fill(pT_pho,weight);
    }
    
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))<1.479)  ) {   /// background barrel ///
      BkgEff_vs_eta_deno_barr->Fill(sc_eta_pho,weight);
      BkgEff_etaPt_deno_barr->Fill(sc_eta_pho,pT_pho,weight);
      if (Gradient>0.0) {
	BkgEff_vs_eta_num_barr->Fill(sc_eta_pho,weight);
	BkgEff_etaPt_num_barr->Fill(sc_eta_pho,pT_pho,weight);
      //   h_weight_bkg->Fill(weight);
  }
    }
    
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// background in barrel ///
      BkgEff_vs_pt_deno_barr->Fill(pT_pho,weight);
      if (Gradient>0.0)  BkgEff_vs_pt_num_barr->Fill(pT_pho,weight);
    }
    
    ///////// endcap ///////
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))>1.479) ) {   /// signal endcap ///
      SigEff_vs_eta_deno_endC->Fill(sc_eta_pho,weight);
      if (Gradient>0.0)  SigEff_vs_eta_num_endC->Fill(sc_eta_pho,weight);
      //  h_weight_sig->Fill(weight);
    }
    
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))>1.479) ) {   /// signal in endcap ///
      SigEff_vs_pt_deno_endC->Fill(pT_pho,weight);
      if (Gradient>0.0)  SigEff_vs_pt_num_endC->Fill(pT_pho,weight);
    }
    
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))>1.479)  ) {   /// background endcap ///
      BkgEff_vs_eta_deno_endC->Fill(sc_eta_pho,weight);
      if (Gradient>0.0)  BkgEff_vs_eta_num_endC->Fill(sc_eta_pho,weight);
      //   h_weight_bkg->Fill(weight);
    }
    
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))>1.479) ) {   /// background in endcap ///
      BkgEff_vs_pt_deno_endC->Fill(pT_pho,weight);
      if (Gradient>0.0)  BkgEff_vs_pt_num_endC->Fill(pT_pho,weight);
    }
  
  
    float pf_pho_iso_Pt = (pf_pho_iso/pT_pho);

    if ((fabs(sc_eta_pho))<1.479) {
      if(genmatched_pho==1) {
	h_r9_sig_barr->Fill(r9_pho);  
	////h_r9_sig_barr->Fill(r9_cleaned_pho,weight);
	h_rho_sig_barr->Fill(rho,weight);
        h_pt_sig_barr->Fill(pT_pho,weight);
	h_sieip_sig_barr->Fill(sieip_pho);
	///////h_sieip_sig_barr->Fill(sieip_cleaned_pho,weight);
	h_pf_pho_iso_sig_barr->Fill(pf_pho_iso,weight);
        h_pf_pho_iso_by_pt_sig_barr->Fill(pf_pho_iso_Pt,weight);
	h_E2byE5_sig_barr->Fill(E2byE5_pho,weight);
	/////h_E2byE5_sig_barr->Fill(E2byE5_cleaned_pho,weight);
	h_etawidth_sig_barr->Fill(etawidth_pho,weight);
	h_phiwidth_sig_barr->Fill(phiwidth_pho,weight);
	h_scRaw_sig_barr->Fill(scRaw_pho,weight);
	h_sieie_sig_barr->Fill(sieie_pho);
	///////h_sieie_sig_barr->Fill(sieie_cleaned_pho,weight);
     	h_chiso_chosenvtx_sig_barr->Fill(pf_charged_iso_chosenvtx,weight);
        h_chiso_vtx1_sig_barr->Fill(pf_charged_iso_vtx1,weight);
        h_chiso_vtx2_sig_barr->Fill(pf_charged_iso_vtx2,weight);
	h_pf_charged_iso_badvtx_sig_barr->Fill(pf_charged_iso_badvtx,weight);
	h_sc_eta_pho_sig_barr->Fill((fabs(sc_eta_pho)),weight);
	h_BDToutput_sig_barr->Fill(Gradient,weight);
      }
      
      if(genmatched_pho==0) {
	h_r9_bkg_barr->Fill(r9_pho);
        //////h_r9_bkg_barr->Fill(r9_cleaned_pho,weight);
        h_rho_bkg_barr->Fill(rho,weight);
        h_pt_bkg_barr->Fill(pT_pho,weight);
	h_sieip_bkg_barr->Fill(sieip_pho);
        ///////h_sieip_bkg_barr->Fill(sieip_cleaned_pho,weight);
	h_pf_pho_iso_bkg_barr->Fill(pf_pho_iso,weight);
        h_pf_pho_iso_by_pt_bkg_barr->Fill(pf_pho_iso_Pt,weight);
	h_E2byE5_bkg_barr->Fill(E2byE5_pho,weight);
        //////h_E2byE5_bkg_barr->Fill(E2byE5_cleaned_pho,weight);
        h_etawidth_bkg_barr->Fill(etawidth_pho,weight);
        h_phiwidth_bkg_barr->Fill(phiwidth_pho,weight);
        h_scRaw_bkg_barr->Fill(scRaw_pho,weight);
	h_sieie_bkg_barr->Fill(sieie_pho);
        /////h_sieie_bkg_barr->Fill(sieie_cleaned_pho,weight);
        h_chiso_chosenvtx_bkg_barr->Fill(pf_charged_iso_chosenvtx,weight);
        h_chiso_vtx1_bkg_barr->Fill(pf_charged_iso_vtx1,weight);
        h_chiso_vtx2_bkg_barr->Fill(pf_charged_iso_vtx2,weight);
        h_pf_charged_iso_badvtx_bkg_barr->Fill(pf_charged_iso_badvtx,weight);
        h_sc_eta_pho_bkg_barr->Fill((fabs(sc_eta_pho)),weight);
        h_BDToutput_bkg_barr->Fill(Gradient,weight);
      }
    }

    if ((fabs(sc_eta_pho))>1.479) {
      if(genmatched_pho==1) {
	h_r9_sig_endC->Fill(r9_pho);
	////h_r9_sig_endC->Fill(r9_cleaned_pho,weight);
        h_rho_sig_endC->Fill(rho,weight);
        h_pt_sig_endC->Fill(pT_pho,weight);
	h_sieip_sig_endC->Fill(sieip_pho);
        /////h_sieip_sig_endC->Fill(sieip_cleaned_pho,weight);
        h_pf_pho_iso_sig_endC->Fill(pf_pho_iso,weight);
        h_pf_pho_iso_by_pt_sig_endC->Fill(pf_pho_iso_Pt,weight);
	h_E2byE5_sig_endC->Fill(E2byE5_pho,weight);
	///////h_E2byE5_sig_endC->Fill(E2byE5_cleaned_pho,weight);
        h_etawidth_sig_endC->Fill(etawidth_pho,weight);
        h_phiwidth_sig_endC->Fill(phiwidth_pho,weight);
        h_scRaw_sig_endC->Fill(scRaw_pho,weight);
	h_sieie_sig_endC->Fill(sieie_pho);
	//////h_sieie_sig_endC->Fill(sieie_cleaned_pho,weight);
        h_chiso_chosenvtx_sig_endC->Fill(pf_charged_iso_chosenvtx,weight);
        h_pf_charged_iso_badvtx_sig_endC->Fill(pf_charged_iso_badvtx,weight);
        h_sc_eta_pho_sig_endC->Fill((fabs(sc_eta_pho)),weight);
        h_BDToutput_sig_endC->Fill(Gradient,weight);
      }

      if(genmatched_pho==0) {
	h_r9_bkg_endC->Fill(r9_pho);
	///////h_r9_bkg_endC->Fill(r9_cleaned_pho,weight);
        h_rho_bkg_endC->Fill(rho,weight);
        h_pt_bkg_endC->Fill(pT_pho,weight);
	h_sieip_bkg_endC->Fill(sieip_pho);
	///////h_sieip_bkg_endC->Fill(sieip_cleaned_pho,weight);
        h_pf_pho_iso_bkg_endC->Fill(pf_pho_iso,weight);
        h_pf_pho_iso_by_pt_bkg_endC->Fill(pf_pho_iso_Pt,weight);
	h_E2byE5_bkg_endC->Fill(E2byE5_pho,weight);
        //////h_E2byE5_bkg_endC->Fill(E2byE5_cleaned_pho,weight);
        h_etawidth_bkg_endC->Fill(etawidth_pho,weight);
        h_phiwidth_bkg_endC->Fill(phiwidth_pho,weight);
	h_scRaw_bkg_endC->Fill(scRaw_pho,weight);
	h_sieie_bkg_endC->Fill(sieie_pho);
	/////h_sieie_bkg_endC->Fill(sieie_cleaned_pho,weight);
        h_chiso_chosenvtx_bkg_endC->Fill(pf_charged_iso_chosenvtx,weight);
        h_pf_charged_iso_badvtx_bkg_endC->Fill(pf_charged_iso_badvtx,weight);
        h_sc_eta_pho_bkg_endC->Fill((fabs(sc_eta_pho)),weight);
        h_BDToutput_bkg_endC->Fill(Gradient,weight);
      }
    }
  }
  
  std::cout << "count_siglike_bkg = " <<  count_siglike_bkg << std::endl ; 

  std::vector<float> roc_barr_sig_eff;
  std::vector<float> roc_barr_bkg_eff;
  std::vector<float> roc_endc_sig_eff;
  std::vector<float> roc_endc_bkg_eff;
  TH1F* hist_roc = new TH1F("hist_roc", "hist_roc", 41, 0.0, 1.0);
  for(int j=0; j<41; j++) {
    std::cout << "\n\n j=" << j << std::endl;
    // std::cout << "counter_barr_sig_n[j]=" << counter_barr_sig_n[j] << "    counter_barr_sig_d[j]=" << counter_barr_sig_d[j] << std::endl;
    //std::cout << "counter_barr_bkg_n[j]=" << counter_barr_bkg_n[j] << "    counter_barr_bkg_d[j]=" << counter_barr_bkg_d[j] << std::endl;
    //std::cout << "sig_eff_barr=" << ((counter_barr_sig_n[j]*1.0)/(counter_barr_sig_d[j]*1.0)) << std::endl;
    //std::cout << "bkg_eff_barr=" << ((counter_barr_bkg_n[j]*1.0)/(counter_barr_bkg_d[j]*1.0)) << std::endl;
    float xval = ((counter_barr_sig_n[j]*1.0)/(counter_barr_sig_d[j]*1.0));
    int ibin = hist_roc->GetBin(xval);
    std::cout << "ibin=" << ibin << "   xval=" << xval <<std::endl;
    hist_roc->SetBinContent(ibin,xval);
    roc_barr_sig_eff.push_back(((counter_barr_sig_n[j]*1.0)/(counter_barr_sig_d[j]*1.0)));
    roc_barr_bkg_eff.push_back(((counter_barr_bkg_n[j]*1.0)/(counter_barr_bkg_d[j]*1.0)));
  }

  /*
  for(int j=0; j<41; j++) {
    //    std::cout << "\n\n j=" << j << std::endl;
    // std::cout << "counter_endc_sig_n[j]=" << counter_endc_sig_n[j] << "    counter_endc_sig_d[j]=" << counter_endc_sig_d[j] << std::endl;
    //std::cout << "counter_endc_bkg_n[j]=" << counter_endc_bkg_n[j] << "    counter_endc_bkg_d[j]=" << counter_endc_bkg_d[j] << std::endl;
    //std::cout << "sig_eff_endc=" << ((counter_endc_sig_n[j]*1.0)/(counter_endc_sig_d[j]*1.0)) << std::endl;
    //std::cout << "bkg_eff_endc=" << ((counter_endc_bkg_n[j]*1.0)/(counter_endc_bkg_d[j]*1.0)) << std::endl;
    roc_endc_sig_eff.push_back(((counter_endc_sig_n[j]*1.0)/(counter_endc_sig_d[j]*1.0)));
    roc_endc_bkg_eff.push_back(((counter_endc_bkg_n[j]*1.0)/(counter_endc_bkg_d[j]*1.0)));
    }*/

  TGraph* roc = new TGraph (roc_barr_sig_eff.size(), &roc_barr_sig_eff[0], &roc_barr_bkg_eff[0]);
  roc->SetMarkerStyle(21);
  roc->SetMarkerColor(kRed);
  roc->SetLineColor(kRed);
  roc->SetLineWidth(2);
  roc->GetHistogram()->SetMaximum(1.3);
  roc->GetHistogram()->SetMinimum(0.0);
  roc->SetTitle("ROC Curve (barrel)");
  roc->GetXaxis()->SetTitle("Signal Efficiency");
  roc->GetYaxis()->SetTitle("Background Efficiency");
 
  TCanvas* a1 = new TCanvas("eff", "eff_s_b");
  a1->cd();
  gPad->SetGrid();
  roc->Draw("AP");
  /*
  TGraph* roc1 = new TGraph (roc_endc_sig_eff.size(), &roc_endc_sig_eff[0], &roc_endc_bkg_eff[0]);
  roc1->SetMarkerStyle(21);
  roc1->SetMarkerColor(kRed);
  roc1->SetLineColor(kRed);
  roc1->SetLineWidth(2);
  roc1->GetHistogram()->SetMaximum(1.3);
  roc1->GetHistogram()->SetMinimum(0.0);
  roc1->SetTitle("ROC Curve (endcap)");
  roc1->GetXaxis()->SetTitle("Signal Efficiency");
  roc1->GetYaxis()->SetTitle("Background Efficiency");

  TCanvas* a2 = new TCanvas("eff1", "eff_s_e");
  a2->cd();
  gPad->SetGrid();
  roc1->Draw("AP");
  */
  outputFile->cd();
  a1->Write();
  roc->Write();
  /// a2->Write();
  ///////  roc1->Write();
  hist_roc->Write();

  TH2D SigEff2D_barr = ((*SigEff_etaPt_num_barr)/(*SigEff_etaPt_deno_barr)) ;
  TH2D BkgEff2D_barr = ((*BkgEff_etaPt_num_barr)/(*BkgEff_etaPt_deno_barr)) ;
  SigEff2D_barr.SetNameTitle("SigEff2D_divideHist", "SigEff2D_divideHist");
  BkgEff2D_barr.SetNameTitle("BkgEff2D_divideHist", "BkgEff2D_divideHist");

  SigEff2D_barr.SetOption("colz");
  BkgEff2D_barr.SetOption("colz");
 
  R9_siglikeBkg_barr->Write();
  BDT_siglikeBkg_barr->Write();
  scRaw_siglikeBkg_barr->Write();
  scEta_siglikeBkg_barr->Write();
  sieie_siglikeBkg_barr->Write();
  sieip_siglikeBkg_barr->Write();
  etawidth_siglikeBkg_barr->Write();
  phiwidth_siglikeBkg_barr->Write();
  ChIsoChosen_siglikeBkg_barr->Write();

  h_r9_sig_barr->Write();
  h_rho_sig_barr->Write();
  h_pt_sig_barr->Write();
  h_sieip_sig_barr->Write();
  h_pf_pho_iso_sig_barr->Write();
  h_pf_pho_iso_by_pt_sig_barr->Write();
  h_E2byE5_sig_barr->Write();
  h_etawidth_sig_barr->Write();
  h_phiwidth_sig_barr->Write();
  h_scRaw_sig_barr->Write();
  h_sieie_sig_barr->Write();
  h_chiso_chosenvtx_sig_barr->Write();
  h_chiso_vtx1_sig_barr->Write();
  h_chiso_vtx2_sig_barr->Write();
  h_pf_charged_iso_badvtx_sig_barr->Write();
  h_sc_eta_pho_sig_barr->Write();

  h_r9_sig_endC->Write();
  h_rho_sig_endC->Write();
  h_pt_sig_endC->Write();
  h_sieip_sig_endC->Write();
  h_pf_pho_iso_sig_endC->Write();
  h_pf_pho_iso_by_pt_sig_endC->Write();
  h_E2byE5_sig_endC->Write();
  h_etawidth_sig_endC->Write();
  h_phiwidth_sig_endC->Write();
  h_scRaw_sig_endC->Write();
  h_sieie_sig_endC->Write();
  h_chiso_chosenvtx_sig_endC->Write();
  h_pf_charged_iso_badvtx_sig_endC->Write();
  h_sc_eta_pho_sig_endC->Write();

  h_r9_bkg_barr->Write();
  h_rho_bkg_barr->Write();
  h_pt_bkg_barr->Write();
  h_sieip_bkg_barr->Write();
  h_pf_pho_iso_bkg_barr->Write();
  h_pf_pho_iso_by_pt_bkg_barr->Write();
  h_E2byE5_bkg_barr->Write();
  h_etawidth_bkg_barr->Write();
  h_phiwidth_bkg_barr->Write();
  h_scRaw_bkg_barr->Write();
  h_sieie_bkg_barr->Write();
  h_chiso_chosenvtx_bkg_barr->Write();
  h_chiso_vtx1_bkg_barr->Write();
  h_chiso_vtx2_bkg_barr->Write();
  h_pf_charged_iso_badvtx_bkg_barr->Write();
  h_sc_eta_pho_bkg_barr->Write();

  h_r9_bkg_endC->Write();
  h_rho_bkg_endC->Write();
  h_pt_bkg_endC->Write();
  h_sieip_bkg_endC->Write();
  h_pf_pho_iso_bkg_endC->Write();
  h_pf_pho_iso_by_pt_bkg_endC->Write();
  h_E2byE5_bkg_endC->Write();
  h_etawidth_bkg_endC->Write();
  h_phiwidth_bkg_endC->Write();
  h_scRaw_bkg_endC->Write();
  h_sieie_bkg_endC->Write();
  h_chiso_chosenvtx_bkg_endC->Write();
  h_pf_charged_iso_badvtx_bkg_endC->Write();
  h_sc_eta_pho_bkg_endC->Write();

  h_BDToutput_sig_barr->Write();
  h_BDToutput_bkg_barr->Write();

  h_BDToutput_sig_endC->Write();
  h_BDToutput_bkg_endC->Write();

  SigEff_vs_eta_deno_barr->Write();
  SigEff_vs_eta_num_barr->Write();
  BkgEff_vs_eta_deno_barr->Write();
  BkgEff_vs_eta_num_barr->Write();
  SigEff_vs_pt_deno_barr->Write();
  SigEff_vs_pt_num_barr->Write();
  BkgEff_vs_pt_deno_barr->Write();
  BkgEff_vs_pt_num_barr->Write();

  SigEff_vs_eta_deno_endC->Write();
  SigEff_vs_eta_num_endC->Write();
  BkgEff_vs_eta_deno_endC->Write();
  BkgEff_vs_eta_num_endC->Write();
  SigEff_vs_pt_deno_endC->Write();
  SigEff_vs_pt_num_endC->Write();
  BkgEff_vs_pt_deno_endC->Write();
  BkgEff_vs_pt_num_endC->Write();

  SigEff_etaPt_deno_barr->Write();
  SigEff_etaPt_num_barr->Write(); 
  BkgEff_etaPt_deno_barr->Write();
  BkgEff_etaPt_num_barr->Write(); 

  SigEff2D_barr.Write("SigEff2D_barr");
  BkgEff2D_barr.Write("BkgEff2D_barr");

  h_weight_sig->Write();
  h_weight_bkg->Write();

  outputFile->Close();
}


						 
