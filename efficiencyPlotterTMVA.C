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
#include "TLegend.h"
#include <TProfile.h>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"


int efficiencyPlotterTMVA() {
  
  gROOT->SetBatch(kTRUE); 
  //////  TFile *input  = new TFile("Hgg_StudyTMVA_Aug_General.root");
  TFile *input  = new TFile("Hgg_StudyTMVA_timeCleaned_AdHocPresel.root");

  //////// Output File /////////
  ///////TFile* outputFile = new TFile("CompareTMVAplots_Aug_general.root","RECREATE");
  TFile* outputFile = new TFile("CompareTMVAplots_Aug_timeCleaned.root","RECREATE");
  Double_t x1 =.72, y1= 0.72,x2 = x1+0.25, y2= y1+0.20;
  
  //// Efficiency ////
  /// Barrel ///
  TH1F*  hist_Sig_eff_num_Eta_barr  = (TH1F*)input->Get("Sig_eff_num_vs_eta_barr");
  TH1F*  hist_Sig_eff_deno_Eta_barr = (TH1F*)input->Get("Sig_eff_deno_vs_eta_barr");
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

  TH1F*  hist_Bkg_eff_num_Eta_barr  = (TH1F*)input->Get("Bkg_eff_num_vs_eta_barr");
  TH1F*  hist_Bkg_eff_deno_Eta_barr = (TH1F*)input->Get("Bkg_eff_deno_vs_eta_barr");
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

  TCanvas* a1 = new TCanvas("eff_eta_barr", "eff_scEta_barr");
  a1->cd();
  gPad->SetGrid();
  SignalEff_vs_Eta_barr->Draw("AP");
  BkgEff_vs_Eta_barr->Draw("P:same");
  a1->Update();

  TLegend* L01 = new TLegend(x1,y1,x2,y2);
  L01->AddEntry(SignalEff_vs_Eta_barr, "Signal Efficiency", "LP");
  L01->AddEntry(BkgEff_vs_Eta_barr, "Background Efficiency", "LP");
  L01->SetFillColor(0);
  L01->Draw();
  a1->Update();
  a1->Write();

  
  //// Efficiency ////
  /// Endcap ///
  TH1F*  hist_Sig_eff_num_Eta_endC  = (TH1F*)input->Get("Sig_eff_num_vs_eta_endC");
  TH1F*  hist_Sig_eff_deno_Eta_endC = (TH1F*)input->Get("Sig_eff_deno_vs_eta_endC");
  TGraphAsymmErrors* SignalEff_vs_Eta_endC = new TGraphAsymmErrors(hist_Sig_eff_num_Eta_endC, hist_Sig_eff_deno_Eta_endC, "cp");
  SignalEff_vs_Eta_endC->SetMarkerStyle(21);
  SignalEff_vs_Eta_endC->SetMarkerColor(kRed);
  SignalEff_vs_Eta_endC->SetLineColor(kRed);
  SignalEff_vs_Eta_endC->SetLineWidth(2);
  SignalEff_vs_Eta_endC->GetHistogram()->SetMaximum(1.3);
  SignalEff_vs_Eta_endC->GetHistogram()->SetMinimum(0.0);
  SignalEff_vs_Eta_endC->SetTitle("Efficiency vs SC eta (BDT cut=0.0) ENDCAP");
  SignalEff_vs_Eta_endC->GetXaxis()->SetTitle("SC Eta");
  SignalEff_vs_Eta_endC->GetYaxis()->SetTitle("Efficiency");

  TH1F*  hist_Bkg_eff_num_Eta_endC  = (TH1F*)input->Get("Bkg_eff_num_vs_eta_endC");
  TH1F*  hist_Bkg_eff_deno_Eta_endC = (TH1F*)input->Get("Bkg_eff_deno_vs_eta_endC");
  TGraphAsymmErrors* BkgEff_vs_Eta_endC = new TGraphAsymmErrors(hist_Bkg_eff_num_Eta_endC, hist_Bkg_eff_deno_Eta_endC, "cp");
  BkgEff_vs_Eta_endC->SetMarkerStyle(20);
  BkgEff_vs_Eta_endC->SetMarkerColor(kBlue);
  BkgEff_vs_Eta_endC->SetLineColor(kBlue);
  BkgEff_vs_Eta_endC->SetLineWidth(2);
  BkgEff_vs_Eta_endC->GetHistogram()->SetMaximum(1.3);
  BkgEff_vs_Eta_endC->GetHistogram()->SetMinimum(0.0);
  BkgEff_vs_Eta_endC->SetTitle("Efficiency vs SC eta (BDT cut=0.0)");
  BkgEff_vs_Eta_endC->GetXaxis()->SetTitle("SC Eta");
  BkgEff_vs_Eta_endC->GetYaxis()->SetTitle("Efficiency");

  TCanvas* a1_endC = new TCanvas("eff_eta_endC", "eff_scEta_endC");
  a1_endC->cd();
  gPad->SetGrid();
  SignalEff_vs_Eta_endC->Draw("AP");
  BkgEff_vs_Eta_endC->Draw("P:same");
  a1_endC->Update();

  TLegend* L01_endC = new TLegend(x1,y1,x2,y2);
  L01_endC->AddEntry(SignalEff_vs_Eta_endC, "Signal Efficiency", "LP");
  L01_endC->AddEntry(BkgEff_vs_Eta_endC, "Background Efficiency", "LP");
  L01_endC->SetFillColor(0);
  L01_endC->Draw();
  a1_endC->Update();
  a1_endC->Write();

  //// Eff vs pt ////
  TH1F*  hist_Sig_eff_num_pt_barr  = (TH1F*)input->Get("Sig_eff_num_vs_pt_barr");
  TH1F*  hist_Sig_eff_deno_pt_barr = (TH1F*)input->Get("Sig_eff_deno_vs_pt_barr");
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


  TH1F*  hist_Bkg_eff_num_pt_barr  = (TH1F*)input->Get("Bkg_eff_num_vs_pt_barr");
  TH1F*  hist_Bkg_eff_deno_pt_barr = (TH1F*)input->Get("Bkg_eff_deno_vs_pt_barr");
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

  TCanvas* a11 = new TCanvas("eff_pt", "eff_pt");
  a11->cd();
  gPad->SetGrid();
  SignalEff_vs_pt_barr->Draw("AP");
  BkgEff_vs_pt_barr->Draw("P:same");
  a11->Update();

  TLegend* L021 = new TLegend(x1,y1,x2,y2);
  L021->AddEntry(SignalEff_vs_pt_barr, "Signal Efficiency", "LP");
  L021->AddEntry(BkgEff_vs_pt_barr, "Background Efficiency", "LP");
  L021->SetFillColor(0);
  L021->Draw();
  a11->Update();
  a11->Write();

  /////// 2D Eff ////////
  /// Barrel ///

  TH2D*  hist2D_Sig_eff_num_barr  = (TH2D*)input->Get("Sig_eff_num_vs_etaPt_barr");
  TH2D*  hist2D_Sig_eff_deno_barr = (TH2D*)input->Get("Sig_eff_deno_vs_etaPt_barr");
  TGraphAsymmErrors* SignalEff2D_barr = new TGraphAsymmErrors(hist2D_Sig_eff_num_barr, hist2D_Sig_eff_deno_barr, "cp");
  SignalEff2D_barr->SetMarkerStyle(21);
  SignalEff2D_barr->SetMarkerColor(kRed);
  SignalEff2D_barr->SetLineColor(kRed);
  SignalEff2D_barr->SetLineWidth(2);
  SignalEff2D_barr->GetHistogram()->SetMaximum(1.3);
  SignalEff2D_barr->GetHistogram()->SetMinimum(0.0);
  SignalEff2D_barr->SetTitle("Efficiency vs SC eta (BDT cut=0.0) BARREL");
  SignalEff2D_barr->GetXaxis()->SetTitle("SC Eta");
  SignalEff2D_barr->GetYaxis()->SetTitle("Efficiency");

  TH2D*  hist2D_Bkg_eff_num_barr  = (TH2D*)input->Get("Sig_eff_num_vs_etaPt_barr");
  TH2D*  hist2D_Bkg_eff_deno_barr = (TH2D*)input->Get("Bkg_eff_deno_vs_etaPt_barr");
  TGraphAsymmErrors* BkgEff2D_barr = new TGraphAsymmErrors(hist2D_Bkg_eff_num_barr, hist2D_Bkg_eff_deno_barr, "cp");
  BkgEff2D_barr->SetMarkerStyle(20);
  BkgEff2D_barr->SetMarkerColor(kBlue);
  BkgEff2D_barr->SetLineColor(kBlue);
  BkgEff2D_barr->SetLineWidth(2);
  BkgEff2D_barr->GetHistogram()->SetMaximum(1.3);
  BkgEff2D_barr->GetHistogram()->SetMinimum(0.0);
  BkgEff2D_barr->SetTitle("Efficiency vs SC eta (BDT cut=0.0)");
  BkgEff2D_barr->GetXaxis()->SetTitle("SC Eta");
  BkgEff2D_barr->GetYaxis()->SetTitle("Efficiency");

  TCanvas* a12D = new TCanvas("eff2D_barr", "eff2D_barr");
  a12D->cd();
  gPad->SetGrid();
  SignalEff2D_barr->Draw("AP");
  BkgEff2D_barr->Draw("P:same");
  a12D->Update();

  TLegend* L012D = new TLegend(x1,y1,x2,y2);
  L012D->AddEntry(SignalEff2D_barr, "Signal Efficiency", "LP");
  L012D->AddEntry(BkgEff2D_barr, "Background Efficiency", "LP");
  L012D->SetFillColor(0);
  L012D->Draw();
  a12D->Update();
  a12D->Write();

  /// pt  ///
  std::cout << "pt" << std::endl;
  TCanvas* zbbb1 = new TCanvas("pt_b", "pt_barr");

  TH1F*  pt_sig_barrel = (TH1F*)input->Get("pt_Sig_b");
  pt_sig_barrel->Scale(1.0/pt_sig_barrel->Integral());
  pt_sig_barrel->SetMarkerStyle(22);
  pt_sig_barrel->SetMarkerColor(kRed);
  pt_sig_barrel->SetLineColor(kRed);
  pt_sig_barrel->SetMarkerSize(1.0);
  pt_sig_barrel->SetLineWidth(2);
  pt_sig_barrel->SetTitle("pt (barrel)");
  pt_sig_barrel->SetMaximum(0.09);

  zbbb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pt_sig_barrel->DrawCopy("P");
  zbbb1->Update();

  TH1F*  pt_bkg_barrel = (TH1F*)input->Get("pt_Bkg_b");
  pt_bkg_barrel->Scale(1.0/pt_bkg_barrel->Integral());
  pt_bkg_barrel->SetMarkerStyle(23);
  pt_bkg_barrel->SetMarkerColor(kBlue);
  pt_bkg_barrel->SetLineColor(kBlue);
  pt_bkg_barrel->SetMarkerSize(1.0);
  pt_bkg_barrel->SetLineWidth(2);
  pt_bkg_barrel->SetTitle("pt (barrel)");
  pt_bkg_barrel->SetMaximum(0.05);

  zbbb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pt_bkg_barrel->DrawCopy("P:sames");
  zbbb1->Update();

  TLegend* LL114 = new TLegend(x1,y1,x2,y2);
  LL114->AddEntry(pt_sig_barrel, "Signal", "L");
  LL114->AddEntry(pt_bkg_barrel, "Background", "L");
  LL114->SetFillColor(0);
  LL114->Draw();
  zbbb1->Write();

  /// endcap - pt ///
  std::cout << "pt endcap" << std::endl;
  TCanvas* z1 = new TCanvas("pt_e", "pt_endc");

  TH1F*  pt_sig_endcap = (TH1F*)input->Get("pt_Sig_e");
  pt_sig_endcap->Scale(1.0/pt_sig_endcap->Integral());
  pt_sig_endcap->SetMarkerStyle(22);
  pt_sig_endcap->SetMarkerColor(kRed);
  pt_sig_endcap->SetLineColor(kRed);
  pt_sig_endcap->SetMarkerSize(1.0);
  pt_sig_endcap->SetLineWidth(2);
  pt_sig_endcap->SetTitle("pt (endcap)");
  pt_sig_endcap->SetMaximum(0.05);

  z1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pt_sig_endcap->DrawCopy("P");
  z1->Update();

  TH1F*  pt_bkg_endcap = (TH1F*)input->Get("pt_Bkg_e");
  pt_bkg_endcap->Scale(1.0/pt_bkg_endcap->Integral());
  pt_bkg_endcap->SetMarkerStyle(23);
  pt_bkg_endcap->SetMarkerColor(kBlue);
  pt_bkg_endcap->SetLineColor(kBlue);
  pt_bkg_endcap->SetMarkerSize(1.0);
  pt_bkg_endcap->SetLineWidth(2);
  pt_bkg_endcap->SetTitle("pt (endcap)");
  pt_bkg_endcap->SetMaximum(0.05);

  z1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pt_bkg_endcap->DrawCopy("P:sames");
  z1->Update();

  TLegend* LL1 = new TLegend(x1,y1,x2,y2);
  LL1->AddEntry(pt_sig_endcap, "Signal", "L");
  LL1->AddEntry(pt_bkg_endcap, "Background", "L");
  LL1->SetFillColor(0);
  LL1->Draw();
  z1->Write();

  /// etawidth  ///
  std::cout << "etawidth" << std::endl;
  TCanvas* bbb1 = new TCanvas("etawidth_b", "etawidth_barr");

  TH1F*  etawidth_sig_barrel = (TH1F*)input->Get("etawidth_Sig_b");
  etawidth_sig_barrel->Scale(1.0/etawidth_sig_barrel->Integral());
  etawidth_sig_barrel->SetMarkerStyle(22);
  etawidth_sig_barrel->SetMarkerColor(kRed);
  etawidth_sig_barrel->SetLineColor(kRed);
  etawidth_sig_barrel->SetMarkerSize(1.0);
  etawidth_sig_barrel->SetLineWidth(2);
  etawidth_sig_barrel->SetTitle("etawidth (barrel)");
  bbb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  etawidth_sig_barrel->DrawCopy("P");
  bbb1->Update();

  TH1F*  etawidth_bkg_barrel = (TH1F*)input->Get("etawidth_Bkg_b");
  etawidth_bkg_barrel->Scale(1.0/etawidth_bkg_barrel->Integral());
  etawidth_bkg_barrel->SetMarkerStyle(23);
  etawidth_bkg_barrel->SetMarkerColor(kBlue);
  etawidth_bkg_barrel->SetLineColor(kBlue);
  etawidth_bkg_barrel->SetMarkerSize(1.0);
  etawidth_bkg_barrel->SetLineWidth(2);
  etawidth_bkg_barrel->SetTitle("etawidth (barrel)");
  bbb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  etawidth_bkg_barrel->DrawCopy("P:sames");
  bbb1->Update();

  TLegend* LL11 = new TLegend(x1,y1,x2,y2);
  LL11->AddEntry(etawidth_sig_barrel, "Signal", "L");
  LL11->AddEntry(etawidth_bkg_barrel, "Background", "L");
  LL11->SetFillColor(0);
  LL11->Draw();
  bbb1->Write();

  //// etawidth - endcap ////
  std::cout << "etawidth - endcap" << std::endl;
  TCanvas* bbb12 = new TCanvas("etawidth_e", "etawidth_endc");

  TH1F*  etawidth_sig_endcap = (TH1F*)input->Get("etawidth_Sig_e");
  etawidth_sig_endcap->Scale(1.0/etawidth_sig_endcap->Integral());
  etawidth_sig_endcap->SetMarkerStyle(22);
  etawidth_sig_endcap->SetMarkerColor(kRed);
  etawidth_sig_endcap->SetLineColor(kRed);
  etawidth_sig_endcap->SetMarkerSize(1.0);
  etawidth_sig_endcap->SetLineWidth(2);
  etawidth_sig_endcap->SetTitle("etawidth (endcap)");
  bbb12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  etawidth_sig_endcap->DrawCopy("P");
  bbb12->Update();

  TH1F*  etawidth_bkg_endcap = (TH1F*)input->Get("etawidth_Bkg_e");
  etawidth_bkg_endcap->Scale(1.0/etawidth_bkg_endcap->Integral());
  etawidth_bkg_endcap->SetMarkerStyle(23);
  etawidth_bkg_endcap->SetMarkerColor(kBlue);
  etawidth_bkg_endcap->SetLineColor(kBlue);
  etawidth_bkg_endcap->SetMarkerSize(1.0);
  etawidth_bkg_endcap->SetLineWidth(2);
  etawidth_bkg_endcap->SetTitle("etawidth (endcap)");
  bbb12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  etawidth_bkg_endcap->DrawCopy("P:sames");
  bbb12->Update();

  TLegend* LL112 = new TLegend(x1,y1,x2,y2);
  LL112->AddEntry(etawidth_sig_endcap, "Signal", "L");
  LL112->AddEntry(etawidth_bkg_endcap, "Background", "L");
  LL112->SetFillColor(0);
  LL112->Draw();
  bbb12->Write();


  ///////// E2byE5  /////////
  std::cout << "E2byE5" << std::endl;
  TCanvas* cvs1 = new TCanvas("E2byE5_b", "E2byE5_barr");

  TH1F*  E2byE5_sig_barrel = (TH1F*)input->Get("E2byE5_Sig_b");
  E2byE5_sig_barrel->Scale(1.0/E2byE5_sig_barrel->Integral());
  E2byE5_sig_barrel->SetMarkerStyle(22);
  E2byE5_sig_barrel->SetMarkerColor(kRed);
  E2byE5_sig_barrel->SetLineColor(kRed);
  E2byE5_sig_barrel->SetMarkerSize(1.0);
  E2byE5_sig_barrel->SetLineWidth(2);
  E2byE5_sig_barrel->SetTitle("E2byE5 (barrel)");
  cvs1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  E2byE5_sig_barrel->DrawCopy("P");
  cvs1->Update();

  TH1F*  E2byE5_bkg_barrel = (TH1F*)input->Get("E2byE5_Bkg_b");
  E2byE5_bkg_barrel->Scale(1.0/E2byE5_bkg_barrel->Integral());
  E2byE5_bkg_barrel->SetMarkerStyle(23);
  E2byE5_bkg_barrel->SetMarkerColor(kBlue);
  E2byE5_bkg_barrel->SetLineColor(kBlue);
  E2byE5_bkg_barrel->SetMarkerSize(1.0);
  E2byE5_bkg_barrel->SetLineWidth(2);
  E2byE5_bkg_barrel->SetTitle("E2byE5 (barrel)");
  cvs1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  E2byE5_bkg_barrel->DrawCopy("P:sames");
  cvs1->Update();

  TLegend* E2byE5_L1 = new TLegend(x1,y1,x2,y2);
  E2byE5_L1->AddEntry(E2byE5_sig_barrel, "Signal", "L");
  E2byE5_L1->AddEntry(E2byE5_bkg_barrel, "Background", "L");
  E2byE5_L1->SetFillColor(0);
  E2byE5_L1->Draw();
  cvs1->Write();

  //// E2byE5 - endcap ////
  std::cout << "E2byE5 - endcap" << std::endl;
  TCanvas* cvs12 = new TCanvas("E2byE5_e", "E2byE5_endc");

  TH1F*  E2byE5_sig_endcap = (TH1F*)input->Get("E2byE5_Sig_e");
  E2byE5_sig_endcap->Scale(1.0/E2byE5_sig_endcap->Integral());
  E2byE5_sig_endcap->SetMarkerStyle(22);
  E2byE5_sig_endcap->SetMarkerColor(kRed);
  E2byE5_sig_endcap->SetLineColor(kRed);
  E2byE5_sig_endcap->SetMarkerSize(1.0);
  E2byE5_sig_endcap->SetLineWidth(2);
  E2byE5_sig_endcap->SetTitle("E2byE5 (endcap)");
  cvs12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  E2byE5_sig_endcap->DrawCopy("P");
  cvs12->Update();

  TH1F*  E2byE5_bkg_endcap = (TH1F*)input->Get("E2byE5_Bkg_e");
  E2byE5_bkg_endcap->Scale(1.0/E2byE5_bkg_endcap->Integral());
  E2byE5_bkg_endcap->SetMarkerStyle(23);
  E2byE5_bkg_endcap->SetMarkerColor(kBlue);
  E2byE5_bkg_endcap->SetLineColor(kBlue);
  E2byE5_bkg_endcap->SetMarkerSize(1.0);
  E2byE5_bkg_endcap->SetLineWidth(2);
  E2byE5_bkg_endcap->SetTitle("E2byE5 (endcap)");
  cvs12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  E2byE5_bkg_endcap->DrawCopy("P:sames");
  cvs12->Update();

  TLegend* E2byE5_L12 = new TLegend(x1,y1,x2,y2);
  E2byE5_L12->AddEntry(E2byE5_sig_endcap, "Signal", "L");
  E2byE5_L12->AddEntry(E2byE5_bkg_endcap, "Background", "L");
  E2byE5_L12->SetFillColor(0);
  E2byE5_L12->Draw();
  cvs12->Write();



  /// sc_eta_pho  ///
  std::cout << "sc_eta_pho" << std::endl;
  TCanvas* qwa1 = new TCanvas("sc_eta_pho_b", "sc_eta_pho_barr");

  TH1F*  sc_eta_pho_sig_barrel = (TH1F*)input->Get("sc_eta_pho_Sig_b");
  sc_eta_pho_sig_barrel->Scale(1.0/sc_eta_pho_sig_barrel->Integral());
  sc_eta_pho_sig_barrel->SetMarkerStyle(22);
  sc_eta_pho_sig_barrel->SetMarkerColor(kRed);
  sc_eta_pho_sig_barrel->SetLineColor(kRed);
  sc_eta_pho_sig_barrel->SetMarkerSize(1.0);
  sc_eta_pho_sig_barrel->SetLineWidth(2);
  sc_eta_pho_sig_barrel->SetMaximum(0.022);
  sc_eta_pho_sig_barrel->SetTitle("sc_eta_pho (barrel)");
  qwa1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  sc_eta_pho_sig_barrel->DrawCopy("P");
  qwa1->Update();

  TH1F*  sc_eta_pho_bkg_barrel = (TH1F*)input->Get("sc_eta_pho_Bkg_b");
  sc_eta_pho_bkg_barrel->Scale(1.0/sc_eta_pho_bkg_barrel->Integral());
  sc_eta_pho_bkg_barrel->SetMarkerStyle(23);
  sc_eta_pho_bkg_barrel->SetMarkerColor(kBlue);
  sc_eta_pho_bkg_barrel->SetLineColor(kBlue);
  sc_eta_pho_bkg_barrel->SetMarkerSize(1.0);
  sc_eta_pho_bkg_barrel->SetLineWidth(2);
  sc_eta_pho_bkg_barrel->SetMaximum(0.022);
  sc_eta_pho_bkg_barrel->SetTitle("sc_eta_pho (barrel)");
  qwa1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  sc_eta_pho_bkg_barrel->DrawCopy("P:sames");
  qwa1->Update();

  TLegend* LK13 = new TLegend(x1,y1,x2,y2);
  LK13->AddEntry(sc_eta_pho_sig_barrel, "Signal", "L");
  LK13->AddEntry(sc_eta_pho_bkg_barrel, "Background", "L");
  LK13->SetFillColor(0);
  LK13->Draw();
  qwa1->Write();


  /// sc_eta_pho  endcap ///
  std::cout << "sc_eta_pho : endcap" << std::endl;
  TCanvas* q1 = new TCanvas("sc_eta_pho_e", "sc_eta_pho_endc");

  TH1F*  sc_eta_pho_sig_endcap = (TH1F*)input->Get("sc_eta_pho_Sig_e");
  sc_eta_pho_sig_endcap->Scale(1.0/sc_eta_pho_sig_endcap->Integral());
  sc_eta_pho_sig_endcap->SetMarkerStyle(22);
  sc_eta_pho_sig_endcap->SetMarkerColor(kRed);
  sc_eta_pho_sig_endcap->SetLineColor(kRed);
  sc_eta_pho_sig_endcap->SetMarkerSize(1.0);
  sc_eta_pho_sig_endcap->SetLineWidth(2);
  sc_eta_pho_sig_endcap->SetMaximum(0.045);
  sc_eta_pho_sig_endcap->SetTitle("sc_eta_pho (endcap)");
  q1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  sc_eta_pho_sig_endcap->DrawCopy("P");
  q1->Update();

  TH1F*  sc_eta_pho_bkg_endcap = (TH1F*)input->Get("sc_eta_pho_Bkg_e");
  sc_eta_pho_bkg_endcap->Scale(1.0/sc_eta_pho_bkg_endcap->Integral());
  sc_eta_pho_bkg_endcap->SetMarkerStyle(23);
  sc_eta_pho_bkg_endcap->SetMarkerColor(kBlue);
  sc_eta_pho_bkg_endcap->SetLineColor(kBlue);
  sc_eta_pho_bkg_endcap->SetMarkerSize(1.0);
  sc_eta_pho_bkg_endcap->SetLineWidth(2);
  sc_eta_pho_bkg_endcap->SetMaximum(0.045);
  sc_eta_pho_bkg_endcap->SetTitle("sc_eta_pho (endcap)");
  q1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  sc_eta_pho_bkg_endcap->DrawCopy("P:sames");
  q1->Update();

  TLegend* K13 = new TLegend(x1,y1,x2,y2);
  K13->AddEntry(sc_eta_pho_sig_endcap, "Signal", "L");
  K13->AddEntry(sc_eta_pho_bkg_endcap, "Background", "L");
  K13->SetFillColor(0);
  K13->Draw();
  q1->Write();


  /// scRaw  ///
  std::cout << "scRAW" << std::endl;
  TCanvas* nnn1 = new TCanvas("scRaw_b", "scRaw_barr");

  TH1F*  scRaw_sig_barrel = (TH1F*)input->Get("scRaw_Sig_b");
  scRaw_sig_barrel->Scale(1.0/scRaw_sig_barrel->Integral());
  scRaw_sig_barrel->SetMarkerStyle(22);
  scRaw_sig_barrel->SetMarkerColor(kRed);
  scRaw_sig_barrel->SetLineColor(kRed);
  scRaw_sig_barrel->SetMarkerSize(1.0);
  scRaw_sig_barrel->SetLineWidth(2);
  scRaw_sig_barrel->SetTitle("scRaw (barrel)");
  nnn1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  scRaw_sig_barrel->DrawCopy("P");
  nnn1->Update();

  TH1F*  scRaw_bkg_barrel = (TH1F*)input->Get("scRaw_Bkg_b");
  scRaw_bkg_barrel->Scale(1.0/scRaw_bkg_barrel->Integral());
  scRaw_bkg_barrel->SetMarkerStyle(23);
  scRaw_bkg_barrel->SetMarkerColor(kBlue);
  scRaw_bkg_barrel->SetLineColor(kBlue);
  scRaw_bkg_barrel->SetMarkerSize(1.0);
  scRaw_bkg_barrel->SetLineWidth(2);
  scRaw_bkg_barrel->SetTitle("scRaw (barrel)");
  nnn1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  scRaw_bkg_barrel->DrawCopy("P:sames");
  nnn1->Update();

  TLegend* LK11 = new TLegend(x1,y1,x2,y2);
  LK11->AddEntry(scRaw_sig_barrel, "Signal", "L");
  LK11->AddEntry(scRaw_bkg_barrel, "Background", "L");
  LK11->SetFillColor(0);
  LK11->Draw();
  nnn1->Write();


  /// scRaw : Endcap  ///
  std::cout << "scRAW : Endcap" << std::endl;
  TCanvas* n1 = new TCanvas("scRaw_e", "scRaw_endc");

  TH1F*  scRaw_sig_endcap = (TH1F*)input->Get("scRaw_Sig_e");
  scRaw_sig_endcap->Scale(1.0/scRaw_sig_endcap->Integral());
  scRaw_sig_endcap->SetMarkerStyle(22);
  scRaw_sig_endcap->SetMarkerColor(kRed);
  scRaw_sig_endcap->SetLineColor(kRed);
  scRaw_sig_endcap->SetMarkerSize(1.0);
  scRaw_sig_endcap->SetLineWidth(2);
  scRaw_sig_endcap->SetTitle("scRaw (endcap)");
  n1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  scRaw_sig_endcap->DrawCopy("P");
  n1->Update();

  TH1F*  scRaw_bkg_endcap = (TH1F*)input->Get("scRaw_Bkg_e");
  scRaw_bkg_endcap->Scale(1.0/scRaw_bkg_endcap->Integral());
  scRaw_bkg_endcap->SetMarkerStyle(23);
  scRaw_bkg_endcap->SetMarkerColor(kBlue);
  scRaw_bkg_endcap->SetLineColor(kBlue);
  scRaw_bkg_endcap->SetMarkerSize(1.0);
  scRaw_bkg_endcap->SetLineWidth(2);
  scRaw_bkg_endcap->SetTitle("scRaw (endcap)");
  n1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  scRaw_bkg_endcap->DrawCopy("P:sames");
  n1->Update();

  TLegend* K10 = new TLegend(x1,y1,x2,y2);
  K10->AddEntry(scRaw_sig_endcap, "Signal", "L");
  K10->AddEntry(scRaw_bkg_endcap, "Background", "L");
  K10->SetFillColor(0);
  K10->Draw();
  n1->Write();



  /// R9 ///
  std::cout << "r9" << std::endl;
  TCanvas* b1 = new TCanvas("R9_b", "R9_barr");

  TH1F*  R9_sig_barrel = (TH1F*)input->Get("r9_Sig_b");
  R9_sig_barrel->Scale(1.0/R9_sig_barrel->Integral());
  R9_sig_barrel->SetMarkerStyle(22);
  R9_sig_barrel->SetMarkerColor(kRed);
  R9_sig_barrel->SetLineColor(kRed);
  R9_sig_barrel->SetMarkerSize(1.0);
  R9_sig_barrel->SetLineWidth(2);
  R9_sig_barrel->SetTitle("R9 (barrel)");
  b1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  R9_sig_barrel->DrawCopy("P");
  b1->Update();

  TH1F*  R9_bkg_barrel = (TH1F*)input->Get("r9_Bkg_b");
  R9_bkg_barrel->Scale(1.0/R9_bkg_barrel->Integral());
  R9_bkg_barrel->SetMarkerStyle(23);
  R9_bkg_barrel->SetMarkerColor(kBlue);
  R9_bkg_barrel->SetLineColor(kBlue);
  R9_bkg_barrel->SetMarkerSize(1.0);
  R9_bkg_barrel->SetLineWidth(2);
  R9_bkg_barrel->SetTitle("R9 (barrel)");
  b1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  R9_bkg_barrel->DrawCopy("P:sames");
  b1->Update();

  TLegend* L1 = new TLegend(x1,y1,x2,y2);
  L1->AddEntry(R9_sig_barrel, "Signal", "L");
  L1->AddEntry(R9_bkg_barrel, "Background", "L");
  L1->SetFillColor(0);
  L1->Draw();
  b1->Write();


  /// R9 endcap ///
  std::cout << "r9 : endcap" << std::endl;
  TCanvas* qb1 = new TCanvas("R9_e", "R9_endc");

  TH1F*  R9_sig_endcap = (TH1F*)input->Get("r9_Sig_e");
  R9_sig_endcap->Scale(1.0/R9_sig_endcap->Integral());
  R9_sig_endcap->SetMarkerStyle(22);
  R9_sig_endcap->SetMarkerColor(kRed);
  R9_sig_endcap->SetLineColor(kRed);
  R9_sig_endcap->SetMarkerSize(1.0);
  R9_sig_endcap->SetLineWidth(2);
  R9_sig_endcap->SetTitle("R9 (endcap)");
  qb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  R9_sig_endcap->DrawCopy("P");
  qb1->Update();

  TH1F*  R9_bkg_endcap = (TH1F*)input->Get("r9_Bkg_e");
  R9_bkg_endcap->Scale(1.0/R9_bkg_endcap->Integral());
  R9_bkg_endcap->SetMarkerStyle(23);
  R9_bkg_endcap->SetMarkerColor(kBlue);
  R9_bkg_endcap->SetLineColor(kBlue);
  R9_bkg_endcap->SetMarkerSize(1.0);
  R9_bkg_endcap->SetLineWidth(2);
  R9_bkg_endcap->SetTitle("R9 (endcap)");
  qb1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  R9_bkg_endcap->DrawCopy("P:sames");
  qb1->Update();

  TLegend* L19 = new TLegend(x1,y1,x2,y2);
  L19->AddEntry(R9_sig_endcap, "Signal", "L");
  L19->AddEntry(R9_bkg_endcap, "Background", "L");
  L19->SetFillColor(0);
  L19->Draw();
  qb1->Write();


  //// phiwidth ///
  std::cout << "phiwidth" << std::endl;

  TCanvas* ccc222 = new TCanvas("PHIWIDTH_b", "PHIWIDTH_barr");

  TH1F*  PHIWIDTH_sig_barrel = (TH1F*)input->Get("phiwidth_Sig_b");
  PHIWIDTH_sig_barrel->Scale(1.0/PHIWIDTH_sig_barrel->Integral());
  PHIWIDTH_sig_barrel->SetMarkerStyle(22);
  PHIWIDTH_sig_barrel->SetMarkerColor(kRed);
  PHIWIDTH_sig_barrel->SetLineColor(kRed);
  PHIWIDTH_sig_barrel->SetMarkerSize(1.0);
  PHIWIDTH_sig_barrel->SetLineWidth(2);
  PHIWIDTH_sig_barrel->SetTitle("PHIWIDTH (barrel)");
  ccc222->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  PHIWIDTH_sig_barrel->DrawCopy("P");
  ccc222->Update();

  TH1F*  PHIWIDTH_bkg_barrel = (TH1F*)input->Get("phiwidth_Bkg_b");
  PHIWIDTH_bkg_barrel->Scale(1.0/PHIWIDTH_bkg_barrel->Integral());
  PHIWIDTH_bkg_barrel->SetMarkerStyle(23);
  PHIWIDTH_bkg_barrel->SetMarkerColor(kBlue);
  PHIWIDTH_bkg_barrel->SetLineColor(kBlue);
  PHIWIDTH_bkg_barrel->SetMarkerSize(1.0);
  PHIWIDTH_bkg_barrel->SetLineWidth(2);
  PHIWIDTH_bkg_barrel->SetTitle("PHIWIDTH (barrel)");
  ccc222->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  PHIWIDTH_bkg_barrel->DrawCopy("P:sames");
  ccc222->Update();

  TLegend* LL778 = new TLegend(x1,y1,x2,y2);
  LL778->AddEntry(PHIWIDTH_sig_barrel, "Signal", "L");
  LL778->AddEntry(PHIWIDTH_bkg_barrel, "Background", "L");
  LL778->SetFillColor(0);
  LL778->Draw();
  ccc222->Write();


  ///// endcap phiwidth /////
  std::cout << "phiwidth" << std::endl;
  TCanvas* c2 = new TCanvas("PHIWIDTH_e", "PHIWIDTH_endcap");

  TH1F*  PHIWIDTH_sig_endcap = (TH1F*)input->Get("phiwidth_Sig_e");
  PHIWIDTH_sig_endcap->Scale(1.0/PHIWIDTH_sig_endcap->Integral());
  PHIWIDTH_sig_endcap->SetMarkerStyle(22);
  PHIWIDTH_sig_endcap->SetMarkerColor(kRed);
  PHIWIDTH_sig_endcap->SetLineColor(kRed);
  PHIWIDTH_sig_endcap->SetMarkerSize(1.0);
  PHIWIDTH_sig_endcap->SetLineWidth(2);
  PHIWIDTH_sig_endcap->SetTitle("PHIWIDTH (endcap)");
  c2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  PHIWIDTH_sig_endcap->DrawCopy("P");
  c2->Update();

  TH1F*  PHIWIDTH_bkg_endcap = (TH1F*)input->Get("phiwidth_Bkg_e");
  PHIWIDTH_bkg_endcap->Scale(1.0/PHIWIDTH_bkg_endcap->Integral());
  PHIWIDTH_bkg_endcap->SetMarkerStyle(23);
  PHIWIDTH_bkg_endcap->SetMarkerColor(kBlue);
  PHIWIDTH_bkg_endcap->SetLineColor(kBlue);
  PHIWIDTH_bkg_endcap->SetMarkerSize(1.0);
  PHIWIDTH_bkg_endcap->SetLineWidth(2);
  PHIWIDTH_bkg_endcap->SetTitle("PHIWIDTH (endcap)");
  c2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  PHIWIDTH_bkg_endcap->DrawCopy("P:sames");
  c2->Update();

  TLegend* L78 = new TLegend(x1,y1,x2,y2);
  L78->AddEntry(PHIWIDTH_sig_endcap, "Signal", "L");
  L78->AddEntry(PHIWIDTH_bkg_endcap, "Background", "L");
  L78->SetFillColor(0);
  L78->Draw();
  c2->Write();


  ///// pf_charged_iso_chosenvtx ///////
  std::cout << "pf_charged_iso_chosenvtx" << std::endl;
  TCanvas* b33 = new TCanvas("pf_charged_iso_chosenvtx_b", "pf_charged_iso_chosenvtx_barr");

  TH1F*  pf_charged_iso_chosenvtx_sig_barrel = (TH1F*)input->Get("pf_charged_iso_chosenvtx_Sig_b");
  pf_charged_iso_chosenvtx_sig_barrel->Scale(1.0/pf_charged_iso_chosenvtx_sig_barrel->Integral());
  pf_charged_iso_chosenvtx_sig_barrel->SetMarkerStyle(22);
  pf_charged_iso_chosenvtx_sig_barrel->SetMarkerColor(kRed);
  pf_charged_iso_chosenvtx_sig_barrel->SetLineColor(kRed);
  pf_charged_iso_chosenvtx_sig_barrel->SetMarkerSize(1.0);
  pf_charged_iso_chosenvtx_sig_barrel->SetLineWidth(2);
  pf_charged_iso_chosenvtx_sig_barrel->SetTitle("pf_charged_iso_chosenvtx (barrel)");
  b33->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_chosenvtx_sig_barrel->DrawCopy();
  b33->Update();

  TH1F*  pf_charged_iso_chosenvtx_bkg_barrel = (TH1F*)input->Get("pf_charged_iso_chosenvtx_Bkg_b");
  pf_charged_iso_chosenvtx_bkg_barrel->Scale(1.0/pf_charged_iso_chosenvtx_bkg_barrel->Integral());
  pf_charged_iso_chosenvtx_bkg_barrel->SetMarkerStyle(23);
  pf_charged_iso_chosenvtx_bkg_barrel->SetMarkerColor(kBlue);
  pf_charged_iso_chosenvtx_bkg_barrel->SetLineColor(kBlue);
  pf_charged_iso_chosenvtx_bkg_barrel->SetMarkerSize(1.0);
  pf_charged_iso_chosenvtx_bkg_barrel->SetLineWidth(2);
  pf_charged_iso_chosenvtx_bkg_barrel->SetTitle("pf_charged_iso_chosenvtx (barrel)");
  b33->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_chosenvtx_bkg_barrel->DrawCopy("sames");
  b33->Update();

  TLegend* L77 = new TLegend(x1,y1,x2,y2);
  L77->AddEntry(pf_charged_iso_chosenvtx_sig_barrel, "Signal", "L");
  L77->AddEntry(pf_charged_iso_chosenvtx_bkg_barrel, "Background", "L");
  L77->SetFillColor(0);
  L77->Draw();
  b33->Write();


  ///// pf_charged_iso_chosenvtx endcap ///////
  std::cout << "pf_charged_iso_chosenvtx -> endcap" << std::endl;
  TCanvas* b330 = new TCanvas("pf_charged_iso_chosenvtx_e", "pf_charged_iso_chosenvtx_endc");

  TH1F*  pf_charged_iso_chosenvtx_sig_endcap = (TH1F*)input->Get("pf_charged_iso_chosenvtx_Sig_e");
  pf_charged_iso_chosenvtx_sig_endcap->Scale(1.0/pf_charged_iso_chosenvtx_sig_endcap->Integral());
  pf_charged_iso_chosenvtx_sig_endcap->SetMarkerStyle(22);
  pf_charged_iso_chosenvtx_sig_endcap->SetMarkerColor(kRed);
  pf_charged_iso_chosenvtx_sig_endcap->SetLineColor(kRed);
  pf_charged_iso_chosenvtx_sig_endcap->SetMarkerSize(1.0);
  pf_charged_iso_chosenvtx_sig_endcap->SetLineWidth(2);
  pf_charged_iso_chosenvtx_sig_endcap->SetTitle("pf_charged_iso_chosenvtx (endcap)");
  b330->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_chosenvtx_sig_endcap->DrawCopy();
  b330->Update();

  TH1F*  pf_charged_iso_chosenvtx_bkg_endcap = (TH1F*)input->Get("pf_charged_iso_chosenvtx_Bkg_e");
  pf_charged_iso_chosenvtx_bkg_endcap->Scale(1.0/pf_charged_iso_chosenvtx_bkg_endcap->Integral());
  pf_charged_iso_chosenvtx_bkg_endcap->SetMarkerStyle(23);
  pf_charged_iso_chosenvtx_bkg_endcap->SetMarkerColor(kBlue);
  pf_charged_iso_chosenvtx_bkg_endcap->SetLineColor(kBlue);
  pf_charged_iso_chosenvtx_bkg_endcap->SetMarkerSize(1.0);
  pf_charged_iso_chosenvtx_bkg_endcap->SetLineWidth(2);
  pf_charged_iso_chosenvtx_bkg_endcap->SetTitle("pf_charged_iso_chosenvtx (endcap)");
  b330->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_chosenvtx_bkg_endcap->DrawCopy("sames");
  b330->Update();

  TLegend* L777 = new TLegend(x1,y1,x2,y2);
  L777->AddEntry(pf_charged_iso_chosenvtx_sig_endcap, "Signal", "L");
  L777->AddEntry(pf_charged_iso_chosenvtx_bkg_endcap, "Background", "L");
  L777->SetFillColor(0);
  L777->Draw();
  b330->Write();

  ///////////////////////// vtx1 //////////////////////////
  std::cout << "pf_charged_iso_vtx1" << std::endl;
  TCanvas* b331 = new TCanvas("pf_charged_iso_vtx1_b", "pf_charged_iso_vtx1_barr");

  TH1F*  pf_charged_iso_vtx1_sig_barrel = (TH1F*)input->Get("pf_charged_iso_vtx1_Sig_b");
  pf_charged_iso_vtx1_sig_barrel->Scale(1.0/pf_charged_iso_vtx1_sig_barrel->Integral());
  pf_charged_iso_vtx1_sig_barrel->SetMarkerStyle(22);
  pf_charged_iso_vtx1_sig_barrel->SetMarkerColor(kRed);
  pf_charged_iso_vtx1_sig_barrel->SetLineColor(kRed);
  pf_charged_iso_vtx1_sig_barrel->SetMarkerSize(1.0);
  pf_charged_iso_vtx1_sig_barrel->SetLineWidth(2);
  pf_charged_iso_vtx1_sig_barrel->SetTitle("pf_charged_iso_vtx1 (barrel)");
  b331->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_vtx1_sig_barrel->DrawCopy();
  b331->Update();

  TH1F*  pf_charged_iso_vtx1_bkg_barrel = (TH1F*)input->Get("pf_charged_iso_vtx1_Bkg_b");
  pf_charged_iso_vtx1_bkg_barrel->Scale(1.0/pf_charged_iso_vtx1_bkg_barrel->Integral());
  pf_charged_iso_vtx1_bkg_barrel->SetMarkerStyle(23);
  pf_charged_iso_vtx1_bkg_barrel->SetMarkerColor(kBlue);
  pf_charged_iso_vtx1_bkg_barrel->SetLineColor(kBlue);
  pf_charged_iso_vtx1_bkg_barrel->SetMarkerSize(1.0);
  pf_charged_iso_vtx1_bkg_barrel->SetLineWidth(2);
  pf_charged_iso_vtx1_bkg_barrel->SetTitle("pf_charged_iso_vtx1 (barrel)");
  b331->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_vtx1_bkg_barrel->DrawCopy("sames");
  b331->Update();

  TLegend* L7711 = new TLegend(x1,y1,x2,y2);
  L7711->AddEntry(pf_charged_iso_vtx1_sig_barrel, "Signal", "L");
  L7711->AddEntry(pf_charged_iso_vtx1_bkg_barrel, "Background", "L");
  L7711->SetFillColor(0);
  L7711->Draw();
  b331->Write();

  ////////// vtx2 /////////////
  std::cout << "pf_charged_iso_vtx2" << std::endl;
  TCanvas* b3311 = new TCanvas("pf_charged_iso_vtx2_b", "pf_charged_iso_vtx2_barr");

  TH1F*  pf_charged_iso_vtx2_sig_barrel = (TH1F*)input->Get("pf_charged_iso_vtx2_Sig_b");
  pf_charged_iso_vtx2_sig_barrel->Scale(1.0/pf_charged_iso_vtx2_sig_barrel->Integral());
  pf_charged_iso_vtx2_sig_barrel->SetMarkerStyle(22);
  pf_charged_iso_vtx2_sig_barrel->SetMarkerColor(kRed);
  pf_charged_iso_vtx2_sig_barrel->SetLineColor(kRed);
  pf_charged_iso_vtx2_sig_barrel->SetMarkerSize(1.0);
  pf_charged_iso_vtx2_sig_barrel->SetLineWidth(2);
  pf_charged_iso_vtx2_sig_barrel->SetTitle("pf_charged_iso_vtx2 (barrel)");
  b3311->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_vtx2_sig_barrel->DrawCopy();
  b3311->Update();

  TH1F*  pf_charged_iso_vtx2_bkg_barrel = (TH1F*)input->Get("pf_charged_iso_vtx2_Bkg_b");
  pf_charged_iso_vtx2_bkg_barrel->Scale(1.0/pf_charged_iso_vtx2_bkg_barrel->Integral());
  pf_charged_iso_vtx2_bkg_barrel->SetMarkerStyle(23);
  pf_charged_iso_vtx2_bkg_barrel->SetMarkerColor(kBlue);
  pf_charged_iso_vtx2_bkg_barrel->SetLineColor(kBlue);
  pf_charged_iso_vtx2_bkg_barrel->SetMarkerSize(1.0);
  pf_charged_iso_vtx2_bkg_barrel->SetLineWidth(2);
  pf_charged_iso_vtx2_bkg_barrel->SetTitle("pf_charged_iso_vtx2 (barrel)");
  b3311->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_vtx2_bkg_barrel->DrawCopy("sames");
  b3311->Update();

  TLegend* L771 = new TLegend(x1,y1,x2,y2);
  L771->AddEntry(pf_charged_iso_vtx2_sig_barrel, "Signal", "L");
  L771->AddEntry(pf_charged_iso_vtx2_bkg_barrel, "Background", "L");
  L771->SetFillColor(0);
  L771->Draw();
  b3311->Write();



  //////////////////////////// pf_charged_iso_badvtx ///////////////////////////
  std::cout << "pf_charged_iso_badvtx" << std::endl;

  TCanvas* b111 = new TCanvas("pf_charged_iso_badvtx_b", "pf_charged_iso_badvtx_barr");
  TH1F*  pf_charged_iso_badvtx_sig_barrel = (TH1F*)input->Get("pf_charged_iso_badvtx_Sig_b");
  pf_charged_iso_badvtx_sig_barrel->Scale(1.0/pf_charged_iso_badvtx_sig_barrel->Integral());
  pf_charged_iso_badvtx_sig_barrel->SetMarkerStyle(22);
  pf_charged_iso_badvtx_sig_barrel->SetMarkerColor(kRed);
  pf_charged_iso_badvtx_sig_barrel->SetLineColor(kRed);
  pf_charged_iso_badvtx_sig_barrel->SetMarkerSize(1.0);
  pf_charged_iso_badvtx_sig_barrel->SetLineWidth(2);
  pf_charged_iso_badvtx_sig_barrel->SetTitle("pf_charged_iso_badvtx (barrel)");
  b111->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_badvtx_sig_barrel->DrawCopy();
  b111->Update();

  TH1F*  pf_charged_iso_badvtx_bkg_barrel = (TH1F*)input->Get("pf_charged_iso_badvtx_Bkg_b");
  pf_charged_iso_badvtx_bkg_barrel->Scale(1.0/pf_charged_iso_badvtx_bkg_barrel->Integral());
  pf_charged_iso_badvtx_bkg_barrel->SetMarkerStyle(23);
  pf_charged_iso_badvtx_bkg_barrel->SetMarkerColor(kBlue);
  pf_charged_iso_badvtx_bkg_barrel->SetLineColor(kBlue);
  pf_charged_iso_badvtx_bkg_barrel->SetMarkerSize(1.0);
  pf_charged_iso_badvtx_bkg_barrel->SetLineWidth(2);
  pf_charged_iso_badvtx_bkg_barrel->SetTitle("pf_charged_iso_badvtx (barrel)");
  b111->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_badvtx_bkg_barrel->DrawCopy("sames");
  b111->Update();

  TLegend* L66 = new TLegend(x1,y1,x2,y2);
  L66->AddEntry(pf_charged_iso_badvtx_sig_barrel, "Signal", "L");
  L66->AddEntry(pf_charged_iso_badvtx_bkg_barrel, "Background", "L");
  L66->SetFillColor(0);
  L66->Draw();

  b111->Write();

  //////////////// endcap  pf_charged_iso_badvtx ///////////////////////////
  std::cout << "pf_charged_iso_badvtx -> endcap" << std::endl;

  TCanvas* g111 = new TCanvas("pf_charged_iso_badvtx_e", "pf_charged_iso_badvtx_endc");
  TH1F*  pf_charged_iso_badvtx_sig_endcap = (TH1F*)input->Get("pf_charged_iso_badvtx_Sig_e");
  pf_charged_iso_badvtx_sig_endcap->Scale(1.0/pf_charged_iso_badvtx_sig_endcap->Integral());
  pf_charged_iso_badvtx_sig_endcap->SetMarkerStyle(22);
  pf_charged_iso_badvtx_sig_endcap->SetMarkerColor(kRed);
  pf_charged_iso_badvtx_sig_endcap->SetLineColor(kRed);
  pf_charged_iso_badvtx_sig_endcap->SetMarkerSize(1.0);
  pf_charged_iso_badvtx_sig_endcap->SetLineWidth(2);
  pf_charged_iso_badvtx_sig_endcap->SetTitle("pf_charged_iso_badvtx (endcap)");
  g111->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_badvtx_sig_endcap->DrawCopy();
  g111->Update();

  TH1F*  pf_charged_iso_badvtx_bkg_endcap = (TH1F*)input->Get("pf_charged_iso_badvtx_Bkg_e");
  pf_charged_iso_badvtx_bkg_endcap->Scale(1.0/pf_charged_iso_badvtx_bkg_endcap->Integral());
  pf_charged_iso_badvtx_bkg_endcap->SetMarkerStyle(23);
  pf_charged_iso_badvtx_bkg_endcap->SetMarkerColor(kBlue);
  pf_charged_iso_badvtx_bkg_endcap->SetLineColor(kBlue);
  pf_charged_iso_badvtx_bkg_endcap->SetMarkerSize(1.0);
  pf_charged_iso_badvtx_bkg_endcap->SetLineWidth(2);
  pf_charged_iso_badvtx_bkg_endcap->SetTitle("pf_charged_iso_badvtx (endcap)");
  g111->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  pf_charged_iso_badvtx_bkg_endcap->DrawCopy("sames");
  g111->Update();

  TLegend* L666 = new TLegend(x1,y1,x2,y2);
  L666->AddEntry(pf_charged_iso_badvtx_sig_endcap, "Signal", "L");
  L666->AddEntry(pf_charged_iso_badvtx_bkg_endcap, "Background", "L");
  L666->SetFillColor(0);
  L666->Draw();

  g111->Write();


  ///// sieie /////
  std::cout << "sieie" << std::endl;

  TCanvas* d1 = new TCanvas("Sieie_b", "Sieie_barr");

  TH1F*  Sieie_sig_barrel = (TH1F*)input->Get("sieie_Sig_b");
  Sieie_sig_barrel->Scale(1.0/Sieie_sig_barrel->Integral());
  Sieie_sig_barrel->SetMarkerStyle(22);
  Sieie_sig_barrel->SetMarkerColor(kRed);
  Sieie_sig_barrel->SetLineColor(kRed);
  Sieie_sig_barrel->SetMarkerSize(1.0);
  Sieie_sig_barrel->SetLineWidth(2);
  Sieie_sig_barrel->SetTitle("Sieie (barrel)");
  d1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieie_sig_barrel->DrawCopy();
  d1->Update();

  TH1F*  Sieie_bkg_barrel = (TH1F*)input->Get("sieie_Bkg_b");
  Sieie_bkg_barrel->Scale(1.0/Sieie_bkg_barrel->Integral());
  Sieie_bkg_barrel->SetMarkerStyle(23);
  Sieie_bkg_barrel->SetMarkerColor(kBlue);
  Sieie_bkg_barrel->SetLineColor(kBlue);
  Sieie_bkg_barrel->SetMarkerSize(1.0);
  Sieie_bkg_barrel->SetLineWidth(2);
  Sieie_bkg_barrel->SetTitle("Sieie (barrel)");
  d1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieie_bkg_barrel->DrawCopy("sames");
  d1->Update();

  TLegend* L5 = new TLegend(x1,y1,x2,y2);
  L5->AddEntry(Sieie_sig_barrel, "Signal", "L");
  L5->AddEntry(Sieie_bkg_barrel, "Background", "L");
  L5->SetFillColor(0);
  L5->Draw();

  d1->Write();

  ///////// sieie endcap ////////
  std::cout << "sieie endcap" << std::endl;

  TCanvas* dm1 = new TCanvas("Sieie_e", "Sieie_endc");
  TH1F*  Sieie_sig_endcap = (TH1F*)input->Get("sieie_Sig_e");
  Sieie_sig_endcap->Scale(1.0/Sieie_sig_endcap->Integral());
  Sieie_sig_endcap->SetMarkerStyle(22);
  Sieie_sig_endcap->SetMarkerColor(kRed);
  Sieie_sig_endcap->SetLineColor(kRed);
  Sieie_sig_endcap->SetMarkerSize(1.0);
  Sieie_sig_endcap->SetLineWidth(2);
  Sieie_sig_endcap->SetTitle("Sieie (endcap)");
  dm1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieie_sig_endcap->DrawCopy();
  dm1->Update();

  TH1F*  Sieie_bkg_endcap = (TH1F*)input->Get("sieie_Bkg_e");
  Sieie_bkg_endcap->Scale(1.0/Sieie_bkg_endcap->Integral());
  Sieie_bkg_endcap->SetMarkerStyle(23);
  Sieie_bkg_endcap->SetMarkerColor(kBlue);
  Sieie_bkg_endcap->SetLineColor(kBlue);
  Sieie_bkg_endcap->SetMarkerSize(1.0);
  Sieie_bkg_endcap->SetLineWidth(2);
  Sieie_bkg_endcap->SetTitle("Sieie (endcap)");
  dm1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieie_bkg_endcap->DrawCopy("sames");
  dm1->Update();

  TLegend* L52 = new TLegend(x1,y1,x2,y2);
  L52->AddEntry(Sieie_sig_endcap, "Signal", "L");
  L52->AddEntry(Sieie_bkg_endcap, "Background", "L");
  L52->SetFillColor(0);
  L52->Draw();

  dm1->Write();
  /////////


  /// sieip ///
  std::cout << "sieip" << std::endl;

  TCanvas* dd1 = new TCanvas("Sieip_b", "Sieip_barr");

  TH1F*  Sieip_sig_barrel = (TH1F*)input->Get("sieip_pho_Sig_b");
  Sieip_sig_barrel->Scale(1.0/Sieip_sig_barrel->Integral());
  Sieip_sig_barrel->SetMarkerStyle(22);
  Sieip_sig_barrel->SetMarkerColor(kRed);
  Sieip_sig_barrel->SetLineColor(kRed);
  Sieip_sig_barrel->SetMarkerSize(1.0);
  Sieip_sig_barrel->SetLineWidth(2);
  Sieip_sig_barrel->SetTitle("Sieip (barrel)");
  dd1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieip_sig_barrel->DrawCopy();
  dd1->Update();

  TH1F*  Sieip_bkg_barrel = (TH1F*)input->Get("sieip_pho_Bkg_b");
  Sieip_bkg_barrel->Scale(1.0/Sieip_bkg_barrel->Integral());
  Sieip_bkg_barrel->SetMarkerStyle(23);
  Sieip_bkg_barrel->SetMarkerColor(kBlue);
  Sieip_bkg_barrel->SetLineColor(kBlue);
  Sieip_bkg_barrel->SetMarkerSize(1.0);
  Sieip_bkg_barrel->SetLineWidth(2);
  Sieip_bkg_barrel->SetTitle("Sieip (barrel)");
  dd1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieip_bkg_barrel->DrawCopy("sames");
  dd1->Update();

  TLegend* L25 = new TLegend(x1,y1,x2,y2);
  L25->AddEntry(Sieip_sig_barrel, "Signal", "L");
  L25->AddEntry(Sieip_bkg_barrel, "Background", "L");
  L25->SetFillColor(0);
  L25->Draw();

  dd1->Write();
  /////////

  ///// sieip -> endcap /////
  std::cout << "sieip endcap" << std::endl;

  TCanvas* dd14 = new TCanvas("Sieip_e", "Sieip_endc");

  TH1F*  Sieip_sig_endcap = (TH1F*)input->Get("sieip_pho_Sig_e");
  Sieip_sig_endcap->Scale(1.0/Sieip_sig_endcap->Integral());
  Sieip_sig_endcap->SetMarkerStyle(22);
  Sieip_sig_endcap->SetMarkerColor(kRed);
  Sieip_sig_endcap->SetLineColor(kRed);
  Sieip_sig_endcap->SetMarkerSize(1.0);
  Sieip_sig_endcap->SetLineWidth(2);
  Sieip_sig_endcap->SetTitle("Sieip (endcap)");
  dd14->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieip_sig_endcap->DrawCopy();
  dd1->Update();

  TH1F*  Sieip_bkg_endcap = (TH1F*)input->Get("sieip_pho_Bkg_e");
  Sieip_bkg_endcap->Scale(1.0/Sieip_bkg_endcap->Integral());
  Sieip_bkg_endcap->SetMarkerStyle(23);
  Sieip_bkg_endcap->SetMarkerColor(kBlue);
  Sieip_bkg_endcap->SetLineColor(kBlue);
  Sieip_bkg_endcap->SetMarkerSize(1.0);
  Sieip_bkg_endcap->SetLineWidth(2);
  Sieip_bkg_endcap->SetTitle("Sieip (endcap)");
  dd14->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Sieip_bkg_endcap->DrawCopy("sames");
  dd14->Update();

  TLegend* L251 = new TLegend(x1,y1,x2,y2);
  L251->AddEntry(Sieip_sig_endcap, "Signal", "L");
  L251->AddEntry(Sieip_bkg_endcap, "Background", "L");
  L251->SetFillColor(0);
  L251->Draw();

  dd14->Write();
  /////////


  /// pf_pho_iso ///
  std::cout << "pf_pho_iso" << std::endl;

  TCanvas* pp1 = new TCanvas("Pf_Pho_Iso_b", "Pf_Pho_Iso_barr");

  TH1F*  Pf_Pho_Iso_sig_barrel = (TH1F*)input->Get("pf_pho_iso_Sig_b");
  Pf_Pho_Iso_sig_barrel->Scale(1.0/Pf_Pho_Iso_sig_barrel->Integral());
  Pf_Pho_Iso_sig_barrel->SetMarkerStyle(22);
  Pf_Pho_Iso_sig_barrel->SetMarkerColor(kRed);
  Pf_Pho_Iso_sig_barrel->SetLineColor(kRed);
  Pf_Pho_Iso_sig_barrel->SetMarkerSize(1.0);
  Pf_Pho_Iso_sig_barrel->SetLineWidth(2);
  Pf_Pho_Iso_sig_barrel->SetTitle("Pf_Pho_Iso (barrel)");
  pp1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_sig_barrel->DrawCopy();
  pp1->Update();

  TH1F*  Pf_Pho_Iso_bkg_barrel = (TH1F*)input->Get("pf_pho_iso_Bkg_b");
  Pf_Pho_Iso_bkg_barrel->Scale(1.0/Pf_Pho_Iso_bkg_barrel->Integral());
  Pf_Pho_Iso_bkg_barrel->SetMarkerStyle(23);
  Pf_Pho_Iso_bkg_barrel->SetMarkerColor(kBlue);
  Pf_Pho_Iso_bkg_barrel->SetLineColor(kBlue);
  Pf_Pho_Iso_bkg_barrel->SetMarkerSize(1.0);
  Pf_Pho_Iso_bkg_barrel->SetLineWidth(2);
  Pf_Pho_Iso_bkg_barrel->SetTitle("Pf_Pho_Iso (barrel)");
  pp1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_bkg_barrel->DrawCopy("sames");
  pp1->Update();

  TLegend* L254 = new TLegend(x1,y1,x2,y2);
  L254->AddEntry(Pf_Pho_Iso_sig_barrel, "Signal", "L");
  L254->AddEntry(Pf_Pho_Iso_bkg_barrel, "Background", "L");
  L254->SetFillColor(0);
  L254->Draw();

  pp1->Write();
  /////////

  ///// endcap  pf_pho_iso ///
  std::cout << "pf_pho_iso endcap" << std::endl;
  TCanvas* pp2 = new TCanvas("Pf_Pho_Iso_e", "Pf_Pho_Iso_endc");

  TH1F*  Pf_Pho_Iso_sig_endcap = (TH1F*)input->Get("pf_pho_iso_Sig_e");
  Pf_Pho_Iso_sig_endcap->Scale(1.0/Pf_Pho_Iso_sig_endcap->Integral());
  Pf_Pho_Iso_sig_endcap->SetMarkerStyle(22);
  Pf_Pho_Iso_sig_endcap->SetMarkerColor(kRed);
  Pf_Pho_Iso_sig_endcap->SetLineColor(kRed);
  Pf_Pho_Iso_sig_endcap->SetMarkerSize(1.0);
  Pf_Pho_Iso_sig_endcap->SetLineWidth(2);
  Pf_Pho_Iso_sig_endcap->SetTitle("Pf_Pho_Iso (endcap)");
  pp2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_sig_endcap->DrawCopy();
  pp2->Update();

  TH1F*  Pf_Pho_Iso_bkg_endcap = (TH1F*)input->Get("pf_pho_iso_Bkg_e");
  Pf_Pho_Iso_bkg_endcap->Scale(1.0/Pf_Pho_Iso_bkg_endcap->Integral());
  Pf_Pho_Iso_bkg_endcap->SetMarkerStyle(23);
  Pf_Pho_Iso_bkg_endcap->SetMarkerColor(kBlue);
  Pf_Pho_Iso_bkg_endcap->SetLineColor(kBlue);
  Pf_Pho_Iso_bkg_endcap->SetMarkerSize(1.0);
  Pf_Pho_Iso_bkg_endcap->SetLineWidth(2);
  Pf_Pho_Iso_bkg_endcap->SetTitle("Pf_Pho_Iso (endcap)");
  pp2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_bkg_endcap->DrawCopy("sames");
  pp2->Update();

  TLegend* L2543 = new TLegend(x1,y1,x2,y2);
  L2543->AddEntry(Pf_Pho_Iso_sig_endcap, "Signal", "L");
  L2543->AddEntry(Pf_Pho_Iso_bkg_endcap, "Background", "L");
  L2543->SetFillColor(0);
  L2543->Draw();

  pp2->Write();
  /////////


  //////
  /// pf_pho_iso/pt ///
  std::cout << "pf_pho_iso / pt" << std::endl;

  TCanvas* ppp1 = new TCanvas("Pf_Pho_Iso_pt_b", "Pf_Pho_Iso_pt_barr");

  TH1F*  Pf_Pho_Iso_sig_barrel1 = (TH1F*)input->Get("pf_pho_iso_Pt_Sig_b");
  Pf_Pho_Iso_sig_barrel1->Scale(1.0/Pf_Pho_Iso_sig_barrel1->Integral());
  Pf_Pho_Iso_sig_barrel1->SetMarkerStyle(22);
  Pf_Pho_Iso_sig_barrel1->SetMarkerColor(kRed);
  Pf_Pho_Iso_sig_barrel1->SetLineColor(kRed);
  Pf_Pho_Iso_sig_barrel1->SetMarkerSize(1.0);
  Pf_Pho_Iso_sig_barrel1->SetLineWidth(2);
  Pf_Pho_Iso_sig_barrel1->SetTitle("Pf_Pho_Iso / pt (barrel)");
  ppp1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_sig_barrel1->DrawCopy();
  ppp1->Update();

  TH1F*  Pf_Pho_Iso_bkg_barrel1 = (TH1F*)input->Get("pf_pho_iso_Pt_Bkg_b");
  Pf_Pho_Iso_bkg_barrel1->Scale(1.0/Pf_Pho_Iso_bkg_barrel1->Integral());
  Pf_Pho_Iso_bkg_barrel1->SetMarkerStyle(23);
  Pf_Pho_Iso_bkg_barrel1->SetMarkerColor(kBlue);
  Pf_Pho_Iso_bkg_barrel1->SetLineColor(kBlue);
  Pf_Pho_Iso_bkg_barrel1->SetMarkerSize(1.0);
  Pf_Pho_Iso_bkg_barrel1->SetLineWidth(2);
  Pf_Pho_Iso_bkg_barrel1->SetTitle("Pf_Pho_Iso / pt (barrel)");
  ppp1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_bkg_barrel1->DrawCopy("sames");
  ppp1->Update();

  TLegend* L2541 = new TLegend(x1,y1,x2,y2);
  L2541->AddEntry(Pf_Pho_Iso_sig_barrel1, "Signal", "L");
  L2541->AddEntry(Pf_Pho_Iso_bkg_barrel1, "Background", "L");
  L2541->SetFillColor(0);
  L2541->Draw();

  ppp1->Write();
  /////////

  ///// endcap  pf_pho_iso ///
  std::cout << "pf_pho_iso endcap pt " << std::endl;
  TCanvas* ppp2 = new TCanvas("Pf_Pho_Iso_pt_e", "Pf_Pho_Iso_pt_endc");

  TH1F*  Pf_Pho_Iso_sig_endcap1 = (TH1F*)input->Get("pf_pho_iso_Pt_Sig_e");
  Pf_Pho_Iso_sig_endcap1->Scale(1.0/Pf_Pho_Iso_sig_endcap1->Integral());
  Pf_Pho_Iso_sig_endcap1->SetMarkerStyle(22);
  Pf_Pho_Iso_sig_endcap1->SetMarkerColor(kRed);
  Pf_Pho_Iso_sig_endcap1->SetLineColor(kRed);
  Pf_Pho_Iso_sig_endcap1->SetMarkerSize(1.0);
  Pf_Pho_Iso_sig_endcap1->SetLineWidth(2);
  Pf_Pho_Iso_sig_endcap1->SetTitle("Pf_Pho_Iso / pt (endcap)");
  ppp2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_sig_endcap1->DrawCopy();
  ppp2->Update();

  TH1F*  Pf_Pho_Iso_bkg_endcap1 = (TH1F*)input->Get("pf_pho_iso_Pt_Bkg_e");
  Pf_Pho_Iso_bkg_endcap1->Scale(1.0/Pf_Pho_Iso_bkg_endcap1->Integral());
  Pf_Pho_Iso_bkg_endcap1->SetMarkerStyle(23);
  Pf_Pho_Iso_bkg_endcap1->SetMarkerColor(kBlue);
  Pf_Pho_Iso_bkg_endcap1->SetLineColor(kBlue);
  Pf_Pho_Iso_bkg_endcap1->SetMarkerSize(1.0);
  Pf_Pho_Iso_bkg_endcap1->SetLineWidth(2);
  Pf_Pho_Iso_bkg_endcap1->SetTitle("Pf_Pho_Iso / pt (endcap)");
  ppp2->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Pf_Pho_Iso_bkg_endcap1->DrawCopy("sames");
  ppp2->Update();

  TLegend* L25431 = new TLegend(x1,y1,x2,y2);
  L25431->AddEntry(Pf_Pho_Iso_sig_endcap1, "Signal", "L");
  L25431->AddEntry(Pf_Pho_Iso_bkg_endcap1, "Background", "L");
  L25431->SetFillColor(0);
  L25431->Draw();

  ppp2->Write();
  /////////

  /// rho ///
  std::cout << "rho" << std::endl;

  TCanvas* mm1 = new TCanvas("Rho_b", "Rho_barr");

  TH1F*  Rho_sig_barrel = (TH1F*)input->Get("rho_Sig_b");
  Rho_sig_barrel->Scale(1.0/Rho_sig_barrel->Integral());
  Rho_sig_barrel->SetMarkerStyle(22);
  Rho_sig_barrel->SetMarkerColor(kRed);
  Rho_sig_barrel->SetLineColor(kRed);
  Rho_sig_barrel->SetMarkerSize(1.0);
  Rho_sig_barrel->SetLineWidth(2);
  Rho_sig_barrel->SetTitle("Rho (barrel)");
  mm1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Rho_sig_barrel->DrawCopy();
  mm1->Update();

  TH1F*  Rho_bkg_barrel = (TH1F*)input->Get("rho_Bkg_b");
  Rho_bkg_barrel->Scale(1.0/Rho_bkg_barrel->Integral());
  Rho_bkg_barrel->SetMarkerStyle(23);
  Rho_bkg_barrel->SetMarkerColor(kBlue);
  Rho_bkg_barrel->SetLineColor(kBlue);
  Rho_bkg_barrel->SetMarkerSize(1.0);
  Rho_bkg_barrel->SetLineWidth(2);
  Rho_bkg_barrel->SetTitle("Rho (barrel)");
  mm1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Rho_bkg_barrel->DrawCopy("sames");
  mm1->Update();

  TLegend* LPO1 = new TLegend(x1,y1,x2,y2);
  LPO1->AddEntry(Rho_sig_barrel, "Signal", "L");
  LPO1->AddEntry(Rho_bkg_barrel, "Background", "L");
  LPO1->SetFillColor(0);
  LPO1->Draw();

  mm1->Write();

  ///// endcap rho ///
  std::cout << "rho endcap" << std::endl;
  TCanvas* ww1 = new TCanvas("Rho_e", "Rho_endc");

  TH1F*  Rho_sig_endcap = (TH1F*)input->Get("rho_Sig_e");
  Rho_sig_endcap->Scale(1.0/Rho_sig_endcap->Integral());
  Rho_sig_endcap->SetMarkerStyle(22);
  Rho_sig_endcap->SetMarkerColor(kRed);
  Rho_sig_endcap->SetLineColor(kRed);
  Rho_sig_endcap->SetMarkerSize(1.0);
  Rho_sig_endcap->SetLineWidth(2);
  Rho_sig_endcap->SetTitle("Rho (endcap)");
  ww1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Rho_sig_endcap->DrawCopy();
  ww1->Update();

  TH1F*  Rho_bkg_endcap = (TH1F*)input->Get("rho_Bkg_e");
  Rho_bkg_endcap->Scale(1.0/Rho_bkg_endcap->Integral());
  Rho_bkg_endcap->SetMarkerStyle(23);
  Rho_bkg_endcap->SetMarkerColor(kBlue);
  Rho_bkg_endcap->SetLineColor(kBlue);
  Rho_bkg_endcap->SetMarkerSize(1.0);
  Rho_bkg_endcap->SetLineWidth(2);
  Rho_bkg_endcap->SetTitle("Rho (endcap)");
  ww1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  Rho_bkg_endcap->DrawCopy("sames");
  ww1->Update();

  TLegend* LPO12 = new TLegend(x1,y1,x2,y2);
  LPO12->AddEntry(Rho_sig_endcap, "Signal", "L");
  LPO12->AddEntry(Rho_bkg_endcap, "Background", "L");
  LPO12->SetFillColor(0);
  LPO12->Draw();

  ww1->Write();

  /////////

  ////// Gradient //////
  TCanvas* h1 = new TCanvas("BDT_b", "BDT_barr");

  TH1F*  BDT_sig_barrel = (TH1F*)input->Get("BDToutput_Sig_b");
  BDT_sig_barrel->Scale(1.0/BDT_sig_barrel->Integral());
  BDT_sig_barrel->SetMarkerStyle(22);
  BDT_sig_barrel->SetMarkerColor(kRed);
  BDT_sig_barrel->SetLineColor(kRed);
  BDT_sig_barrel->SetMarkerSize(1.0);
  BDT_sig_barrel->SetLineWidth(2);
  BDT_sig_barrel->SetTitle("BDT (barrel)");
  BDT_sig_barrel->SetMaximum(1.0);

  h1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  BDT_sig_barrel->DrawCopy();
  h1->Update();

  TH1F*  BDT_bkg_barrel = (TH1F*)input->Get("BDToutput_Bkg_b");
  BDT_bkg_barrel->Scale(1.0/BDT_bkg_barrel->Integral());
  BDT_bkg_barrel->SetMarkerStyle(23);
  BDT_bkg_barrel->SetMarkerColor(kBlue);
  BDT_bkg_barrel->SetLineColor(kBlue);
  BDT_bkg_barrel->SetMarkerSize(1.0);
  BDT_bkg_barrel->SetLineWidth(2);
  BDT_bkg_barrel->SetTitle("BDT (barrel)");
  BDT_bkg_barrel->SetMaximum(1.0);

  h1->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  BDT_bkg_barrel->DrawCopy("sames");
  h1->Update();

  TLegend* L13 = new TLegend(x1,y1,x2,y2);
  L13->AddEntry(BDT_sig_barrel, "Signal", "L");
  L13->AddEntry(BDT_bkg_barrel, "Background", "L");
  L13->SetFillColor(0);
  L13->Draw();

  h1->Write();

  //// endcap  Gradient //////
  TCanvas* h12 = new TCanvas("BDT_e", "BDT_endc");

  TH1F*  BDT_sig_endcap = (TH1F*)input->Get("BDToutput_Sig_e");
  BDT_sig_endcap->Scale(1.0/BDT_sig_endcap->Integral());
  BDT_sig_endcap->SetMarkerStyle(22);
  BDT_sig_endcap->SetMarkerColor(kRed);
  BDT_sig_endcap->SetLineColor(kRed);
  BDT_sig_endcap->SetMarkerSize(1.0);
  BDT_sig_endcap->SetLineWidth(2);
  BDT_sig_endcap->SetTitle("BDT (endcap)");
  BDT_sig_endcap->SetMaximum(1.0);

  h12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  BDT_sig_endcap->DrawCopy();
  h12->Update();

  TH1F*  BDT_bkg_endcap = (TH1F*)input->Get("BDToutput_Bkg_e");
  BDT_bkg_endcap->Scale(1.0/BDT_bkg_endcap->Integral());
  BDT_bkg_endcap->SetMarkerStyle(23);
  BDT_bkg_endcap->SetMarkerColor(kBlue);
  BDT_bkg_endcap->SetLineColor(kBlue);
  BDT_bkg_endcap->SetMarkerSize(1.0);
  BDT_bkg_endcap->SetLineWidth(2);
  BDT_bkg_endcap->SetTitle("BDT (endcap)");
  BDT_bkg_endcap->SetMaximum(1.0);

  h12->cd();
  gPad->SetGrid();
  //  gPad->SetLogy();
  BDT_bkg_endcap->DrawCopy("sames");
  h12->Update();

  TLegend* L135 = new TLegend(x1,y1,x2,y2);
  L135->AddEntry(BDT_sig_endcap, "Signal", "L");
  L135->AddEntry(BDT_bkg_endcap, "Background", "L");
  L135->SetFillColor(0);
  L135->Draw();

  h12->Write();


  outputFile->Close();
  return 0;
}
