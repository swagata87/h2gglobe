#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"

#include <iostream>

Int_t pho_n, gp_n;
Short_t gp_mother[10000], gp_status[10000], gp_pdgid[10000];
TClonesArray* pho_p4 = new TClonesArray("TLorentzVector");
TClonesArray* gp_p4 = new TClonesArray("TLorentzVector");

TH1F* pt_prompt  = new TH1F("pt_prompt",  "", 200, 0, 200);
TH1F* eta_prompt = new TH1F("eta_prompt", "", 100, -1.5, 1.5);
TH1F* pt_fake    = new TH1F("pt_fake",    "", 200, 0, 200);
TH1F* eta_fake   = new TH1F("eta_fake",   "", 100, -1.5, 1.5);

Bool_t GenMatchedPhoton(int ipho) { //, float& energy) {

  Bool_t is_prompt = false;
  TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(ipho);
  
  for(int ip=0; ip<gp_n; ++ip) {
    if( gp_status[ip] != 1 || gp_pdgid[ip] != 22 ) {
      continue;
    }
    TLorentzVector * p4 = (TLorentzVector*)gp_p4->At(ip);
    if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
    int mother_id = abs( gp_pdgid[ gp_mother[ip] ] );
   
    if( mother_id <= 25 ) {
      
      float dr = phop4->DeltaR(*p4);
      if (dr<0.3 && fabs((p4->Pt()-phop4->Pt())/p4->Pt()) < 0.5) {
	//std::cout << "Reco: " << phop4->Eta() << " " << phop4->Phi() << " " << dr << " " << std::endl;
	
	is_prompt = true;
	break;
      }
    }
  }
  
  return is_prompt;
}

void test() {

  TChain* chain = new TChain("event");
  chain->Add("/afs/cern.ch/user/s/sani/eos/cms/store/group/phys_higgs/cmshgg/processed/FixPFIsolationIssue/NewClustering/GJet_Pt40_doubleEMEnriched_TuneZ2star_14TeV-pythia6_TP2023SHCALDR-SHCALJan23_PU140BX25_PH2_1K_FB_V6-v2/*.root");

  chain->SetBranchStatus("*", 0);

  chain->SetBranchStatus("pho_n", 1);
  chain->SetBranchStatus("pho_p4", 1);
  chain->SetBranchStatus("gp_n", 1);
  chain->SetBranchStatus("gp_p4", 1);
  chain->SetBranchStatus("gp_mother", 1);
  chain->SetBranchStatus("gp_status", 1);
  chain->SetBranchStatus("gp_pdgid", 1);

  chain->SetBranchAddress("pho_n", &pho_n);
  chain->SetBranchAddress("pho_p4", &pho_p4);
  chain->SetBranchAddress("gp_n", &gp_n);
  chain->SetBranchAddress("gp_p4", &gp_p4);
  chain->SetBranchAddress("gp_mother", &gp_mother);
  chain->SetBranchAddress("gp_status", &gp_status);
  chain->SetBranchAddress("gp_pdgid", &gp_pdgid);
  
  Int_t entries = chain->GetEntries();
  for (int z=0; z<100000; z++) {
    if (z%1000 == 0)
      std::cout << z << std::endl;
    chain->GetEntry(z);
    for (int i=0; i<pho_n; i++) {
      TLorentzVector* phop4 = (TLorentzVector*)pho_p4->At(i);
      
      if (phop4->Pt() > 25. && fabs(phop4->Eta())<1.479) {
	bool is_prompt = false;
	//bool is_prompt = GenMatchedPhoton(i);
	//std::cout << is_prompt << std::endl;
	
	//TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(ipho);
  
	for(int ip=0; ip<gp_n; ++ip) {
	  //std::cout << gp_status[ip] << " " <<  gp_pdgid[ip] << std::endl;
	  if( gp_status[ip] != 1 || gp_pdgid[ip] != 22 ) {
	    continue;
	  }
	  TLorentzVector * p4 = (TLorentzVector*)gp_p4->At(ip);
	  if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
	  int mother_id = abs( gp_pdgid[ gp_mother[ip] ] );
	  
	  if( mother_id <= 25 ) {
	    
	    float dr = phop4->DeltaR(*p4);
	    if (dr<0.3 && fabs((p4->Pt()-phop4->Pt())/p4->Pt()) < 0.5) {
	      //std::cout << "Reco: " << phop4->Eta() << " " << phop4->Phi() << " " << dr << " " << std::endl;
	      
	      is_prompt = true;
	      break;
	    }
	  }
	}
	

	if (is_prompt) {
	  pt_prompt->Fill(phop4->Pt());
	  eta_prompt->Fill(phop4->Eta());
	} else {
	  pt_fake->Fill(phop4->Pt());
	  eta_fake->Fill(phop4->Eta());
	}
      }
    }
  }

  TCanvas* c = new TCanvas("c", "c");
  c->Divide(2, 1);
  c->cd(1);
  pt_prompt->Draw();
  pt_fake->Draw("SAME");
  c->cd(2);
  eta_prompt->Draw();
  eta_fake->Draw("SAME");

  c->SaveAs("phase2.root");
  
}
