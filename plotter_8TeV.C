#include <iostream>
#include <iomanip>
#include <vector>
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPostScript.h"
#include <TPaveStats.h>
#include <TProfile.h>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

int plotter_8TeV() {

  gROOT->SetBatch(kTRUE); 
  //  TFile *input_14TeV_cleaned  = new TFile("Hgg_StudyTMVA_timeCleaned.root");    
  //  TFile *input_14TeV_cleaned_AdhocPresel  = new TFile("Hgg_StudyTMVA_timeCleaned_AdHocPresel.root");

  TFile *input_8TeV_legacy   = new TFile("Hgg_StudyTMVA_8TeV_Legacy.root");
  TFile *input_8TeV_mine       = new TFile("Hgg_General_8TeV_ReproduceLegacy.root");
  TFile *input_8TeV_mine_gr     = new TFile("Hgg_General_8TeV_NewTrainingNewSettings.root");
  TFile *input_8TeV_mine_gr_noNewVar= new TFile("Hgg_General_8TeV_NewTrainingNewSettings_noNewVar.root");

  //  TFile *input_8TeV   = new TFile("Hgg_General_8TeV_Nov19.root");
  //  TFile *input_14TeV  = new TFile("Hgg_General_14TeV_Nov19.root");

  //////// Output File /////////
  TFile* outputFile = new TFile("Reproduce8TeV.root","RECREATE");
  Double_t x1 =.72, y1= 0.72,x2 = x1+0.25, y2= y1+0.20;
  TCanvas* b1z = new TCanvas("roc1", "ROC1");

  TGraph* roc_8_l = (TGraph*)input_8TeV_legacy->Get("Graph");
  roc_8_l->SetLineColor(kBlue);
  roc_8_l->SetMarkerColor(kBlue);
  roc_8_l->SetMarkerSize(1.0);
  roc_8_l->SetMaximum(1.2);
  roc_8_l->GetHistogram()->SetMaximum(1.2);

  TGraph* roc_8_m = (TGraph*)input_8TeV_mine->Get("Graph");
  roc_8_m->SetLineColor(kPink);
  roc_8_m->SetMarkerColor(kPink);
  roc_8_m->SetMarkerSize(1.0);
  roc_8_m->SetMaximum(1.2);
  roc_8_m->GetHistogram()->SetMaximum(1.2);

  TGraph* roc_8_m_gr = (TGraph*)input_8TeV_mine_gr->Get("Graph");
  roc_8_m_gr->SetLineColor(kGreen);
  roc_8_m_gr->SetMarkerColor(kGreen);
  roc_8_m_gr->SetMarkerSize(1.0);
  roc_8_m_gr->SetMaximum(1.2);
  roc_8_m_gr->GetHistogram()->SetMaximum(1.2);

  TGraph* roc_8_m_gr_nnv = (TGraph*)input_8TeV_mine_gr_noNewVar->Get("Graph");
  roc_8_m_gr_nnv->SetLineColor(kCyan);
  roc_8_m_gr_nnv->SetMarkerColor(kCyan);
  roc_8_m_gr_nnv->SetMarkerSize(1.0);
  roc_8_m_gr_nnv->SetMaximum(1.2);
  roc_8_m_gr_nnv->GetHistogram()->SetMaximum(1.2);

  b1z->cd();
  gPad->SetGrid();

  TMultiGraph *mg = new TMultiGraph();
  //  mg->Add(roc_14,"p");
  //  mg->Add(roc_14_adhoc,"p");
  mg->Add(roc_8_l,"p");
  mg->Add(roc_8_m,"p");
  mg->Add(roc_8_m_gr,"p");
  mg->Add(roc_8_m_gr_nnv,"p");

  mg->Draw("ap");
  b1z->Update();

  TLegend* L1z = new TLegend(x1,y1,x2,y2);
  // L1z->AddEntry(roc_14, "14TeV, 140PU (No preselection)", "L");
  L1z->AddEntry(roc_8_l,      "8TeV (Legacy Training)", "L");
  L1z->AddEntry(roc_8_m,      "8TeV (Legacy-like New Training)", "L");
  L1z->AddEntry(roc_8_m_gr_nnv,   "8TeV (New Training with diff BDT config)", "L");
  L1z->AddEntry(roc_8_m_gr,   "8TeV (New Training with 2 new input var and diff BDT config)", "L");

  //  L1z->AddEntry(roc_14_adhoc, "14TeV, 140PU (Very loose preselection)", "L");
  L1z->SetFillColor(0);
  L1z->Draw();
  b1z->Write();


  ////
  //// BDToutput_Sig_b
  ////
  std::cout << "BDT output" << std::endl;
  TCanvas* zbbb1 = new TCanvas("bdt", "bdt");

  TH1F*  bdt_legacy_sig = (TH1F*)input_8TeV_legacy->Get("BDToutput_Sig_b");
  bdt_legacy_sig->Scale(1.0/bdt_legacy_sig->Integral());
  bdt_legacy_sig->SetMarkerColor(kRed);
  bdt_legacy_sig->SetMarkerStyle(20);
  bdt_legacy_sig->SetLineColor(kRed);
  bdt_legacy_sig->SetLineWidth(2);
  bdt_legacy_sig->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_legacy_sig->DrawCopy("P");
  zbbb1->Update();
  TH1F*  bdt_legacy_bkg = (TH1F*)input_8TeV_legacy->Get("BDToutput_Bkg_b");
  bdt_legacy_bkg->Scale(1.0/bdt_legacy_bkg->Integral());
  bdt_legacy_bkg->SetMarkerColor(kBlue);
  bdt_legacy_bkg->SetMarkerStyle(20);
  bdt_legacy_bkg->SetLineColor(kBlue);
  bdt_legacy_bkg->SetLineWidth(2);
  bdt_legacy_bkg->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_legacy_bkg->DrawCopy("P:sames");
  zbbb1->Update();

  TH1F*  bdt_ReproduceLegacy_sig = (TH1F*)input_8TeV_mine->Get("BDToutput_Sig_b");
  bdt_ReproduceLegacy_sig->Scale(1.0/bdt_ReproduceLegacy_sig->Integral());
  bdt_ReproduceLegacy_sig->SetMarkerColor(kCyan);
  bdt_ReproduceLegacy_sig->SetMarkerStyle(21);
  bdt_ReproduceLegacy_sig->SetLineColor(kCyan);
  bdt_ReproduceLegacy_sig->SetLineWidth(2);
  bdt_ReproduceLegacy_sig->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_ReproduceLegacy_sig->DrawCopy("P:sames");
  zbbb1->Update();
  TH1F*  bdt_ReproduceLegacy_bkg = (TH1F*)input_8TeV_mine->Get("BDToutput_Bkg_b");
  bdt_ReproduceLegacy_bkg->Scale(1.0/bdt_ReproduceLegacy_bkg->Integral());
  bdt_ReproduceLegacy_bkg->SetLineColor(kMagenta);
  bdt_ReproduceLegacy_bkg->SetMarkerColor(kMagenta);
  bdt_ReproduceLegacy_bkg->SetMarkerStyle(21);
  bdt_ReproduceLegacy_bkg->SetLineWidth(2);
  bdt_ReproduceLegacy_bkg->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_ReproduceLegacy_bkg->DrawCopy("P:sames");
  zbbb1->Update();

  TH1F*  bdt_mineGradient_sig = (TH1F*)input_8TeV_mine_gr->Get("BDToutput_Sig_b");
  bdt_mineGradient_sig->Scale(1.0/bdt_mineGradient_sig->Integral());
  bdt_mineGradient_sig->SetLineColor(kGreen+1);
  bdt_mineGradient_sig->SetMarkerColor(kGreen+1);
  bdt_mineGradient_sig->SetMarkerStyle(24);
  bdt_mineGradient_sig->SetLineWidth(2);
  bdt_mineGradient_sig->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_mineGradient_sig->DrawCopy("P:sames");
  zbbb1->Update();
  TH1F*  bdt_mineGradient_bkg = (TH1F*)input_8TeV_mine_gr->Get("BDToutput_Bkg_b");
  bdt_mineGradient_bkg->Scale(1.0/bdt_mineGradient_bkg->Integral());
  bdt_mineGradient_bkg->SetLineColor(kOrange+1);
  bdt_mineGradient_bkg->SetMarkerColor(kOrange+1);
  bdt_mineGradient_bkg->SetMarkerStyle(24);
  bdt_mineGradient_bkg->SetLineWidth(2);
  bdt_mineGradient_bkg->SetMaximum(0.05);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_mineGradient_bkg->DrawCopy("P:sames");
  zbbb1->Update();

  TH1F*  bdt_mineGradientNoNewVar_sig = (TH1F*)input_8TeV_mine_gr_noNewVar->Get("BDToutput_Sig_b");
  bdt_mineGradientNoNewVar_sig->Scale(1.0/bdt_mineGradientNoNewVar_sig->Integral());
  bdt_mineGradientNoNewVar_sig->SetLineColor(kMagenta+2);
  bdt_mineGradientNoNewVar_sig->SetMarkerColor(kMagenta+2);
  bdt_mineGradientNoNewVar_sig->SetMarkerStyle(26);
  bdt_mineGradientNoNewVar_sig->SetLineWidth(2);
  bdt_mineGradientNoNewVar_sig->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_mineGradientNoNewVar_sig->DrawCopy("P:sames");
  zbbb1->Update();
  TH1F*  bdt_mineGradientNoNewVar_bkg = (TH1F*)input_8TeV_mine_gr_noNewVar->Get("BDToutput_Bkg_b");
  bdt_mineGradientNoNewVar_bkg->Scale(1.0/bdt_mineGradientNoNewVar_bkg->Integral());
  bdt_mineGradientNoNewVar_bkg->SetLineColor(kBlue+2);
  bdt_mineGradientNoNewVar_bkg->SetMarkerColor(kBlue+2);
  bdt_mineGradientNoNewVar_bkg->SetMarkerStyle(26);
  bdt_mineGradientNoNewVar_bkg->SetLineWidth(2);
  bdt_mineGradientNoNewVar_bkg->SetMaximum(1.0);
  zbbb1->cd();
  gPad->SetGrid();
  bdt_mineGradientNoNewVar_bkg->DrawCopy("P:sames");
  zbbb1->Update();


  ////Legend for BDT output/////
  TLegend* LL114 = new TLegend(x1,y1,x2,y2);
  LL114->AddEntry(bdt_legacy_sig, "Legacy (AdaBoost)", "LP");
  LL114->AddEntry(bdt_ReproduceLegacy_sig, "Legacy-like training (AdaBoost)", "LP");
  LL114->AddEntry(bdt_mineGradientNoNewVar_sig, "New training (GradientBoost)", "LP");
  LL114->AddEntry(bdt_mineGradient_sig, "New training (GradientBoost+2 New var)", "LP");
 
  LL114->AddEntry(bdt_legacy_bkg, "Legacy (AdaBoost)", "LP");
  LL114->AddEntry(bdt_ReproduceLegacy_bkg, "Legacy-like training (AdaBoost)", "LP");
  LL114->AddEntry(bdt_mineGradientNoNewVar_bkg, "New training (GradientBoost)", "LP");
  LL114->AddEntry(bdt_mineGradient_bkg, "New training (GradientBoost+2 New var)", "LP");

  LL114->SetFillColor(0);
  LL114->Draw();
  zbbb1->Write();



  /*
  /////////////////////////////
  //// Efficiency ////                                                                                                                                   
  /// Barrel ///                                         
  TH1F*  hist_Sig_eff_num_Eta_barr  = (TH1F*)input_8TeV_legacyTraining->Get("Sig_eff_num_vs_eta_barr");
  TH1F*  hist_Sig_eff_deno_Eta_barr = (TH1F*)input_8TeV_legacyTraining->Get("Sig_eff_deno_vs_eta_barr");
  TGraphAsymmErrors* SignalEff_vs_Eta_barr = new TGraphAsymmErrors(hist_Sig_eff_num_Eta_barr, hist_Sig_eff_deno_Eta_barr, "cp");
  SignalEff_vs_Eta_barr->SetMarkerStyle(21);
  SignalEff_vs_Eta_barr->SetMarkerColor(kRed);
  SignalEff_vs_Eta_barr->SetLineColor(kRed);
  SignalEff_vs_Eta_barr->SetLineWidth(2);
  SignalEff_vs_Eta_barr->GetHistogram()->SetMaximum(1.3);
  SignalEff_vs_Eta_barr->GetHistogram()->SetMinimum(0.0);
  SignalEff_vs_Eta_barr->SetTitle("Efficiency vs SC eta (BDT cut=0.0) BARREL");
  SignalEff_vs_Eta_barr->GetXaxis()->SetTitle("SC Eta");
  SignalEff_vs_Eta_barr->GetYaxis()->SetTitle("Efficiency");

  TH1F*  hist_Bkg_eff_num_Eta_barr  = (TH1F*)input_8TeV_legacyTraining->Get("Bkg_eff_num_vs_eta_barr");
  TH1F*  hist_Bkg_eff_deno_Eta_barr = (TH1F*)input_8TeV_legacyTraining->Get("Bkg_eff_deno_vs_eta_barr");
  TGraphAsymmErrors* BkgEff_vs_Eta_barr = new TGraphAsymmErrors(hist_Bkg_eff_num_Eta_barr, hist_Bkg_eff_deno_Eta_barr, "cp");
  BkgEff_vs_Eta_barr->SetMarkerStyle(20);
  BkgEff_vs_Eta_barr->SetMarkerColor(kBlue);
  BkgEff_vs_Eta_barr->SetLineColor(kBlue);
  BkgEff_vs_Eta_barr->SetLineWidth(2);
  BkgEff_vs_Eta_barr->GetHistogram()->SetMaximum(1.3);
  BkgEff_vs_Eta_barr->GetHistogram()->SetMinimum(0.0);
  BkgEff_vs_Eta_barr->SetTitle("Efficiency vs SC eta (BDT cut=0.0)");
  BkgEff_vs_Eta_barr->GetXaxis()->SetTitle("SC Eta");
  BkgEff_vs_Eta_barr->GetYaxis()->SetTitle("Efficiency");

  /// New Training ///
  TH1F*  hist_Sig_eff_num_Eta_barr_new  = (TH1F*)input_8TeV_myTraining->Get("Sig_eff_num_vs_eta_barr");
  TH1F*  hist_Sig_eff_deno_Eta_barr_new = (TH1F*)input_8TeV_myTraining->Get("Sig_eff_deno_vs_eta_barr");
  TGraphAsymmErrors* SignalEff_vs_Eta_barr_new = new TGraphAsymmErrors(hist_Sig_eff_num_Eta_barr_new, hist_Sig_eff_deno_Eta_barr_new, "cp");
  SignalEff_vs_Eta_barr_new->SetMarkerStyle(21);
  SignalEff_vs_Eta_barr_new->SetMarkerColor(kGreen);
  SignalEff_vs_Eta_barr_new->SetLineColor(kGreen);
  SignalEff_vs_Eta_barr_new->SetLineWidth(2);
  SignalEff_vs_Eta_barr_new->GetHistogram()->SetMaximum(1.3);
  SignalEff_vs_Eta_barr_new->GetHistogram()->SetMinimum(0.0);
  SignalEff_vs_Eta_barr_new->SetTitle("Efficiency vs SC eta (BDT cut=0.0) BARREL");
  SignalEff_vs_Eta_barr_new->GetXaxis()->SetTitle("SC Eta");
  SignalEff_vs_Eta_barr_new->GetYaxis()->SetTitle("Efficiency");

  TH1F*  hist_Bkg_eff_num_Eta_barr_new  = (TH1F*)input_8TeV_myTraining->Get("Bkg_eff_num_vs_eta_barr");
  TH1F*  hist_Bkg_eff_deno_Eta_barr_new = (TH1F*)input_8TeV_myTraining->Get("Bkg_eff_deno_vs_eta_barr");
  TGraphAsymmErrors* BkgEff_vs_Eta_barr_new = new TGraphAsymmErrors(hist_Bkg_eff_num_Eta_barr_new, hist_Bkg_eff_deno_Eta_barr_new, "cp");
  BkgEff_vs_Eta_barr_new->SetMarkerStyle(20);
  BkgEff_vs_Eta_barr_new->SetMarkerColor(kCyan);
  BkgEff_vs_Eta_barr_new->SetLineColor(kCyan);
  BkgEff_vs_Eta_barr_new->SetLineWidth(2);
  BkgEff_vs_Eta_barr_new->GetHistogram()->SetMaximum(1.3);
  BkgEff_vs_Eta_barr_new->GetHistogram()->SetMinimum(0.0);
  BkgEff_vs_Eta_barr_new->SetTitle("Efficiency vs SC eta (BDT cut=0.0)");
  BkgEff_vs_Eta_barr_new->GetXaxis()->SetTitle("SC Eta");
  BkgEff_vs_Eta_barr_new->GetYaxis()->SetTitle("Efficiency");
  ////
  TCanvas* a1 = new TCanvas("eff_eta_barr", "eff_scEta_barr");
  a1->cd();
  gPad->SetGrid();
  SignalEff_vs_Eta_barr->Draw("AP");
  BkgEff_vs_Eta_barr->Draw("P:same");
  SignalEff_vs_Eta_barr_new->Draw("P:same");
  BkgEff_vs_Eta_barr_new->Draw("P:same");
  a1->Update();

  TLegend* L01 = new TLegend(x1,y1,x2,y2);
  L01->AddEntry(SignalEff_vs_Eta_barr, "Sig Eff (Legacy training)", "LP");
  L01->AddEntry(BkgEff_vs_Eta_barr,    "Bkg Eff (Legacy training)", "LP");
  L01->AddEntry(SignalEff_vs_Eta_barr_new, "Sig Eff (New training)", "LP");
  L01->AddEntry(BkgEff_vs_Eta_barr_new,    "Bkg Eff (New training)", "LP");
  L01->SetFillColor(0);
  L01->Draw();
  a1->Update();
  a1->Write();

  //// Eff vs pt ////   
  //// Legacy training ////
  TH1F*  hist_Sig_eff_num_pt_barr  = (TH1F*)input_8TeV_legacyTraining->Get("Sig_eff_num_vs_pt_barr");
  TH1F*  hist_Sig_eff_deno_pt_barr = (TH1F*)input_8TeV_legacyTraining->Get("Sig_eff_deno_vs_pt_barr");
  TGraphAsymmErrors* SignalEff_vs_pt_barr = new TGraphAsymmErrors(hist_Sig_eff_num_pt_barr, hist_Sig_eff_deno_pt_barr, "cp");
  SignalEff_vs_pt_barr->SetMarkerStyle(21);
  SignalEff_vs_pt_barr->SetMarkerColor(kRed);
  SignalEff_vs_pt_barr->SetLineColor(kRed);
  SignalEff_vs_pt_barr->SetLineWidth(2);
  SignalEff_vs_pt_barr->GetHistogram()->SetMaximum(1.3);
  SignalEff_vs_pt_barr->GetHistogram()->SetMinimum(0.0);
  SignalEff_vs_pt_barr->SetTitle("Efficiency vs pt (BDT cut=0.0) BARREL");
  SignalEff_vs_pt_barr->GetXaxis()->SetTitle("pt");
  SignalEff_vs_pt_barr->GetYaxis()->SetTitle("Efficiency");

  TH1F*  hist_Bkg_eff_num_pt_barr  = (TH1F*)input_8TeV_legacyTraining->Get("Bkg_eff_num_vs_pt_barr");
  TH1F*  hist_Bkg_eff_deno_pt_barr = (TH1F*)input_8TeV_legacyTraining->Get("Bkg_eff_deno_vs_pt_barr");
  TGraphAsymmErrors* BkgEff_vs_pt_barr = new TGraphAsymmErrors(hist_Bkg_eff_num_pt_barr, hist_Bkg_eff_deno_pt_barr, "cp");
  BkgEff_vs_pt_barr->SetMarkerStyle(20);
  BkgEff_vs_pt_barr->SetMarkerColor(kBlue);
  BkgEff_vs_pt_barr->SetLineColor(kBlue);
  BkgEff_vs_pt_barr->SetLineWidth(2);
  BkgEff_vs_pt_barr->GetHistogram()->SetMaximum(1.3);
  BkgEff_vs_pt_barr->GetHistogram()->SetMinimum(0.0);
  BkgEff_vs_pt_barr->SetTitle("Efficiency vs pt (BDT cut=0.0) BARREL");
  BkgEff_vs_pt_barr->GetXaxis()->SetTitle("pt");
  BkgEff_vs_pt_barr->GetYaxis()->SetTitle("Efficiency");

  //// New training ////
  TH1F*  hist_Sig_eff_num_pt_barr_new  = (TH1F*)input_8TeV_myTraining->Get("Sig_eff_num_vs_pt_barr");
  TH1F*  hist_Sig_eff_deno_pt_barr_new = (TH1F*)input_8TeV_myTraining->Get("Sig_eff_deno_vs_pt_barr");
  TGraphAsymmErrors* SignalEff_vs_pt_barr_new = new TGraphAsymmErrors(hist_Sig_eff_num_pt_barr_new, hist_Sig_eff_deno_pt_barr_new, "cp");
  SignalEff_vs_pt_barr_new->SetMarkerStyle(21);
  SignalEff_vs_pt_barr_new->SetMarkerColor(kGreen);
  SignalEff_vs_pt_barr_new->SetLineColor(kGreen);
  SignalEff_vs_pt_barr_new->SetLineWidth(2);
  SignalEff_vs_pt_barr_new->GetHistogram()->SetMaximum(1.3);
  SignalEff_vs_pt_barr_new->GetHistogram()->SetMinimum(0.0);
  SignalEff_vs_pt_barr_new->SetTitle("Efficiency vs pt (BDT cut=0.0) BARREL");
  SignalEff_vs_pt_barr_new->GetXaxis()->SetTitle("pt");
  SignalEff_vs_pt_barr_new->GetYaxis()->SetTitle("Efficiency");

  TH1F*  hist_Bkg_eff_num_pt_barr_new  = (TH1F*)input_8TeV_myTraining->Get("Bkg_eff_num_vs_pt_barr");
  TH1F*  hist_Bkg_eff_deno_pt_barr_new = (TH1F*)input_8TeV_myTraining->Get("Bkg_eff_deno_vs_pt_barr");
  TGraphAsymmErrors* BkgEff_vs_pt_barr_new = new TGraphAsymmErrors(hist_Bkg_eff_num_pt_barr_new, hist_Bkg_eff_deno_pt_barr_new, "cp");
  BkgEff_vs_pt_barr_new->SetMarkerStyle(20);
  BkgEff_vs_pt_barr_new->SetMarkerColor(kCyan);
  BkgEff_vs_pt_barr_new->SetLineColor(kCyan);
  BkgEff_vs_pt_barr_new->SetLineWidth(2);
  BkgEff_vs_pt_barr_new->GetHistogram()->SetMaximum(1.3);
  BkgEff_vs_pt_barr_new->GetHistogram()->SetMinimum(0.0);
  BkgEff_vs_pt_barr_new->SetTitle("Efficiency vs pt (BDT cut=0.0) BARREL");
  BkgEff_vs_pt_barr_new->GetXaxis()->SetTitle("pt");
  BkgEff_vs_pt_barr_new->GetYaxis()->SetTitle("Efficiency");

  TCanvas* a11 = new TCanvas("eff_pt", "eff_pt");
  a11->cd();
  gPad->SetGrid();
  SignalEff_vs_pt_barr->Draw("AP");
  BkgEff_vs_pt_barr->Draw("P:same");
  SignalEff_vs_pt_barr_new->Draw("P:same");
  BkgEff_vs_pt_barr_new->Draw("P:same");
  a11->Update();

  TLegend* L021 = new TLegend(x1,y1,x2,y2);
  L021->AddEntry(SignalEff_vs_pt_barr, "Sig Eff (Legacy training)", "LP");
  L021->AddEntry(BkgEff_vs_pt_barr,    "Bkg Eff (Legacy training)", "LP");
  L021->AddEntry(SignalEff_vs_pt_barr_new, "Sig Eff (New training)", "LP");
  L021->AddEntry(BkgEff_vs_pt_barr_new, "Bkg Eff (New training)", "LP");

  L021->SetFillColor(0);
  L021->Draw();
  a11->Update();
  a11->Write();
  */

  outputFile->Close();
  return 0;
}
