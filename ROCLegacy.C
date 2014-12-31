#define ROCLegacy_cxx
#include "ROCLegacy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TGraph.h"

void ROCLegacy::Loop()
{

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  TFile* outputFile = new TFile("Hgg_StudyTMVA_8TeV_Legacy.root","RECREATE");

  TH1F*  SigEff_vs_eta_deno_barr = new TH1F("Sig_eff_deno_vs_eta_barr", "h_Sig_deno_barr", 40, -2.0, 2.0);
  TH1F*  SigEff_vs_eta_num_barr  = new TH1F("Sig_eff_num_vs_eta_barr",  "h_Sig_num_barr",  40, -2.0, 2.0);
  TH1F*  BkgEff_vs_eta_deno_barr = new TH1F("Bkg_eff_deno_vs_eta_barr", "h_Bkg_deno_barr", 40, -2.0, 2.0);
  TH1F*  BkgEff_vs_eta_num_barr  = new TH1F("Bkg_eff_num_vs_eta_barr",  "h_Bkg_num_barr",  40, -2.0, 2.0);

  TH1F*  SigEff_vs_pt_deno_barr = new TH1F("Sig_eff_deno_vs_pt_barr", "h1_Sig_deno_barr", 100, 0.0, 200.0);
  TH1F*  SigEff_vs_pt_num_barr  = new TH1F("Sig_eff_num_vs_pt_barr",  "h1_Sig_num_barr",  100, 0.0, 200.0);
  TH1F*  BkgEff_vs_pt_deno_barr = new TH1F("Bkg_eff_deno_vs_pt_barr", "h1_Bkg_deno_barr", 100, 0.0, 200.0);
  TH1F*  BkgEff_vs_pt_num_barr  = new TH1F("Bkg_eff_num_vs_pt_barr",  "h1_Bkg_num_barr",  100, 0.0, 200.0);

  TH1F* h_BDToutput_sig_barr = new TH1F("BDToutput_Sig_b", "h_BDToutput_Sig_b", 300, -1.5, 1.5);
  TH1F* h_BDToutput_bkg_barr = new TH1F("BDToutput_Bkg_b", "h_BDToutput_Bkg_b", 300, -1.5, 1.5);

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

    /*    
    counter_endc_sig_n[j]=0;
    counter_endc_bkg_n[j]=0;
    counter_endc_sig_d[j]=0;
    counter_endc_bkg_d[j]=0;
    */
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    /// ROC ///                                                                                                                         
    if ( (fabs(sc_eta_pho)) < 1.479) { /// BARREL ///                                                                                   
      float x=-1.00;
      for (int j=0; j<41; j++) {
	if ((mva_legacy>x || mva_legacy==x) && genmatched_pho==1.0) counter_barr_sig_n[j]++;
	if ((mva_legacy>x || mva_legacy==x) && genmatched_pho==0.0) counter_barr_bkg_n[j]++;
	if (genmatched_pho==1.0) counter_barr_sig_d[j]++;
	if (genmatched_pho==0.0) counter_barr_bkg_d[j]++;
	x=x+0.05;
      }
    }

    /*
    if ( (fabs(sc_eta_pho)) > 1.479) { /// ENDCAP ///                                                                                   
      float x=-1.00;
      for (int j=0; j<41; j++) {
	if ((mva_legacy>x || mva_legacy==x) && genmatched_pho==1.0) counter_endc_sig_n[j]++;
	if ((mva_legacy>x || mva_legacy==x) && genmatched_pho==0.0) counter_endc_bkg_n[j]++;
	if (genmatched_pho==1.0) counter_endc_sig_d[j]++;
	if (genmatched_pho==0.0) counter_endc_bkg_d[j]++;
	x=x+0.05;
      }
    }
    */

    //// Eff ////                                                                                                                                         
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// signal barrel ///                                                                   
      SigEff_vs_eta_deno_barr->Fill(sc_eta_pho,weight);
      if (mva_legacy>0.0)  {
        SigEff_vs_eta_num_barr->Fill(sc_eta_pho,weight);
   	//  h_weight_sig->Fill(weight);                                                                                                                   
      }
    }
    if ((genmatched_pho==1.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// signal in barrel ///                                                                
      SigEff_vs_pt_deno_barr->Fill(pT_pho,weight);
      if (mva_legacy>0.0)  SigEff_vs_pt_num_barr->Fill(pT_pho,weight);
    }
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))<1.479)  ) {   /// background barrel ///                                                              
      BkgEff_vs_eta_deno_barr->Fill(sc_eta_pho,weight);
      if (mva_legacy>0.0) {
        BkgEff_vs_eta_num_barr->Fill(sc_eta_pho,weight);
      	//   h_weight_bkg->Fill(weight);                                                                                                                  
      }
    }
    
    if ((genmatched_pho==0.0) && ((fabs(sc_eta_pho))<1.479) ) {   /// background in barrel ///                                                            
      BkgEff_vs_pt_deno_barr->Fill(pT_pho,weight);
      if (mva_legacy>0.0)  BkgEff_vs_pt_num_barr->Fill(pT_pho,weight);
    }


    if ((fabs(sc_eta_pho))<1.479) {
      if(genmatched_pho==1) {
        h_BDToutput_sig_barr->Fill(mva_legacy,weight);
      }
      if(genmatched_pho==0) {
        h_BDToutput_bkg_barr->Fill(mva_legacy,weight);
      }
    }
  }
  
  std::vector<float> roc_barr_sig_eff;
  std::vector<float> roc_barr_bkg_eff;
  //  std::vector<float> roc_endc_sig_eff;
  //  std::vector<float> roc_endc_bkg_eff;
  
  for(int j=0; j<41; j++) {
    std::cout << "\n\n j=" << j << std::endl;
    std::cout << "counter_barr_sig_n[j]=" << counter_barr_sig_n[j] << "    counter_barr_sig_d[j]=" << counter_barr_sig_d[j] << std::endl;
    std::cout << "counter_barr_bkg_n[j]=" << counter_barr_bkg_n[j] << "    counter_barr_bkg_d[j]=" << counter_barr_bkg_d[j] << std::endl;
    std::cout << "sig_eff_barr=" << ((counter_barr_sig_n[j]*1.0)/(counter_barr_sig_d[j]*1.0)) << std::endl;
    std::cout << "bkg_eff_barr=" << ((counter_barr_bkg_n[j]*1.0)/(counter_barr_bkg_d[j]*1.0)) << std::endl;
    roc_barr_sig_eff.push_back(((counter_barr_sig_n[j]*1.0)/(counter_barr_sig_d[j]*1.0)));
    roc_barr_bkg_eff.push_back(((counter_barr_bkg_n[j]*1.0)/(counter_barr_bkg_d[j]*1.0)));
  }
  
  /*
  for(int j=0; j<41; j++) {
    std::cout << "\n\n j=" << j << std::endl;
    std::cout << "counter_endc_sig_n[j]=" << counter_endc_sig_n[j] << "    counter_endc_sig_d[j]=" << counter_endc_sig_d[j] << std::endl;
    std::cout << "counter_endc_bkg_n[j]=" << counter_endc_bkg_n[j] << "    counter_endc_bkg_d[j]=" << counter_endc_bkg_d[j] << std::endl;
    std::cout << "sig_eff_endc=" << ((counter_endc_sig_n[j]*1.0)/(counter_endc_sig_d[j]*1.0)) << std::endl;
    std::cout << "bkg_eff_endc=" << ((counter_endc_bkg_n[j]*1.0)/(counter_endc_bkg_d[j]*1.0)) << std::endl;
    roc_endc_sig_eff.push_back(((counter_endc_sig_n[j]*1.0)/(counter_endc_sig_d[j]*1.0)));
    roc_endc_bkg_eff.push_back(((counter_endc_bkg_n[j]*1.0)/(counter_endc_bkg_d[j]*1.0)));
  }
  */

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
  //  a2->Write();

  SigEff_vs_eta_deno_barr->Write();
  SigEff_vs_eta_num_barr->Write();
  BkgEff_vs_eta_deno_barr->Write();
  BkgEff_vs_eta_num_barr->Write();

  SigEff_vs_pt_deno_barr->Write();
  SigEff_vs_pt_num_barr->Write();
  BkgEff_vs_pt_deno_barr->Write();
  BkgEff_vs_pt_num_barr->Write();

  h_BDToutput_sig_barr->Write();
  h_BDToutput_bkg_barr->Write();


  outputFile->Close();
  
}
