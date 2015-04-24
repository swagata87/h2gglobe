#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TString.h"
#include "TMath.h"

#include <iostream>

float weight = 1.;
bool passID = true;
Int_t pho_n, gp_n, genjet_algo1_n;
Int_t dipho_gensel, dipho_genfakeind, dipho_genpromptind;
Float_t dipho_genprompte, dipho_genfakee;
Float_t dipho_genprompteta, dipho_genfakeeta;
Float_t dipho_genpromptphi, dipho_genfakephi;

Short_t gp_mother[10000], gp_status[10000], gp_pdgid[10000];
Float_t genjet_algo1_em[150];
TClonesArray* pho_p4 = new TClonesArray("TLorentzVector");
TClonesArray* gp_p4 = new TClonesArray("TLorentzVector");
TClonesArray* genjet_algo1_p4 = new TClonesArray("TLorentzVector");

TH1F* recoEffPromptPt_den = new TH1F("recoEffPromptPt_den", "", 10000, 0, 1000);
TH1F* recoEffPromptPt_num = new TH1F("recoEffPromptPt_num", "", 10000, 0, 1000);
TH1F* idPromptPt_num      = new TH1F("idPromptPt_num", "", 10000, 0, 1000);
TH1F* recoEffFakePt_den   = new TH1F("recoEffFakePt_den", "", 10000, 0, 1000);
TH1F* recoEffFakePt_num   = new TH1F("recoEffFakePt_num", "", 10000, 0, 1000);
TH1F* idFakePt_num        = new TH1F("idFakePt_num", "", 10000, 0, 1000);

TH1F* recoEffPromptEta_den = new TH1F("recoEffPromptEta_den", "", 60, -3, 3);
TH1F* recoEffPromptEta_num = new TH1F("recoEffPromptEta_num", "", 60, -3, 3);
TH1F* idPromptEta_num      = new TH1F("idPromptEta_num", "",      60, -3, 3);
TH1F* recoEffFakeEta_den   = new TH1F("recoEffFakeEta_den", "",   60, -3, 3);
TH1F* recoEffFakeEta_num   = new TH1F("recoEffFakeEta_num", "",   60, -3, 3);
TH1F* idFakeEta_num        = new TH1F("idFakeEta_num", "",        60, -3, 3);
    
TH2F* nGenFakeablePt  = new TH2F("nGenFakeablePt", "", 10000, 0, 1000, 10000, 0, 1000);
TH2F* nGenFakeableEta = new TH2F("nGenFakeableEta", "", 60, -3, 3, 60, -3, 3);

bool isPrompt(int ip) {

  if( gp_status[ip] != 1 || gp_pdgid[ip] != 22 ) {
    return false;
  }
  
  TLorentzVector * p4 = (TLorentzVector*)gp_p4->At(ip);
  if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { 
    return false; 
  }

  int mother_id = abs( gp_pdgid[ gp_mother[ip] ] );
  
  if( mother_id <= 25 ) {
    return true;
  }
  
  return false;  
}

bool isEMGenJet(int ip) {
  
  TLorentzVector * p4 = (TLorentzVector*)genjet_algo1_p4->At(ip);
  if(p4->Pt() < 20. || fabs(p4->Eta()) > 3.) { 
    return false; 
  }
  
  if((genjet_algo1_em[ip]/p4->E()) >.5) {
    return true;
  }
  
  return false;  
}

bool selectGenEvents(int& iPrompt, int& iFake) {
  
  bool promptPho = false;
  
  for(int ip=0; ip<gp_n; ++ip) {
    if (isPrompt(ip)) {
      promptPho = true;
      iPrompt = ip;
      break;
    }
  }
  
  if (promptPho) {
    for (int i=0; i<genjet_algo1_n; i++) {
      TLorentzVector * p4 = (TLorentzVector*)genjet_algo1_p4->At(i);
      if (p4->Pt() > 20 and (genjet_algo1_em[i]/p4->E()) > .5) {
	iFake = i;
	return true;
      }
    }
  }
  
  return false;
}

void promptRecoEfficiency() {
  
  for(int ip=0; ip<gp_n; ++ip) {
    TLorentzVector * p4 = (TLorentzVector*)gp_p4->At(ip);
    
    if (isPrompt(ip)) {
      for (int ipho=0; ipho<pho_n; ipho++) {
	TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(ipho);
	
	if (phop4->Pt() > 25. && fabs(phop4->Eta())<1.479) {  

	  recoEffPromptPt_den->Fill(p4->Pt(), weight);
	  recoEffPromptEta_den->Fill(p4->Eta(), weight);

	  float dr = phop4->DeltaR(*p4);
	  if (dr<0.3 && fabs((p4->Pt()-phop4->Pt())/p4->Pt()) < 0.5) {
	    recoEffPromptPt_num->Fill(p4->Pt(), weight);
	    recoEffPromptEta_num->Fill(p4->Eta(), weight);
	  }
	}
      }
    }
  }
}

float dR(float e1, float e2, float p1, float p2) {

  float dp = std::abs(p1-p2); 
  if (dp > TMath::Pi()) 
    dp -= 2.*TMath::Pi(); 
 
  return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

void nPromptID(int iPho) {
  
  TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(iPho);

  if (phop4->Pt() > 25. && fabs(phop4->Eta())<1.479) {
    float dr = dR(phop4->Eta(), dipho_genprompteta, phop4->Phi(), dipho_genpromptphi);
 
    Float_t promptTheta = 2*atan(exp(-dipho_genprompteta));
    Float_t promptet = dipho_genprompte*sin(promptTheta);
    if (dr<0.3 && fabs((promptet-phop4->Et())/promptet) < 0.5) {
      if (passID) {
	idPromptPt_num->Fill(promptet, weight);
	idPromptEta_num->Fill(dipho_genprompteta, weight);
      }
    }
  }
}

void fakeRecoEfficiency() {
  
  for(int ip=0; ip<genjet_algo1_n; ++ip) {
    TLorentzVector * p4 = (TLorentzVector*)genjet_algo1_p4->At(ip);
    
    if (isEMGenJet(ip)) {
      for (int ipho=0; ipho<pho_n; ipho++) {
	TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(ipho);
	
	if (phop4->Pt() > 25. && fabs(phop4->Eta())<1.479) { 

	  recoEffFakePt_den->Fill(p4->Pt(), weight);
	  recoEffFakeEta_den->Fill(p4->Eta(), weight);
    
	  float dr = phop4->DeltaR(*p4);
	  if (dr<0.3 && fabs((p4->Pt()-phop4->Pt())/p4->Pt()) < 0.5) {	    
	    recoEffFakePt_num->Fill(p4->Pt(), weight);
	    recoEffFakeEta_num->Fill(p4->Eta(), weight);
	  }
	}
      }
    }
  }
}

void nFakeID(int ipho) {

  TLorentzVector* phop4 = (TLorentzVector*) pho_p4->At(ipho);
	
  if (phop4->Pt() > 25. && fabs(phop4->Eta())<1.479) { 

    float dr = dR(phop4->Eta(), dipho_genfakeeta, phop4->Phi(), dipho_genfakephi);
    Float_t fakeTheta = 2*atan(exp(-dipho_genfakeeta));
    Float_t fakeet = dipho_genfakee*sin(fakeTheta);
  
    if (dr<0.3 && fabs((fakeet-phop4->Et())/fakeet) < 0.5) {	    
      if (passID) {
	idFakePt_num->Fill(fakeet, weight);
	idFakeEta_num->Fill(dipho_genfakeeta, weight);
      }
    }
  }
}

void nGenDiphoton() {
  
  Float_t promptTheta = 2*atan(exp(-dipho_genprompteta));
  Float_t fakeTheta = 2*atan(exp(-dipho_genfakeeta));
  Float_t promptet = dipho_genprompte*sin(promptTheta);
  Float_t fakeet   = dipho_genfakee*sin(fakeTheta);
  nGenFakeablePt->Fill(promptet, fakeet, weight);
  nGenFakeableEta->Fill(dipho_genprompteta, dipho_genfakeeta, weight);
}

Bool_t GenMatchedPhoton(int ipho, int& index) {
  
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
	index = ip;
	is_prompt = true;
	break;
      }
    }
  }
  
  return is_prompt;
}

void efficiencyCalculator(bool recoEff=false) {

  TChain* chain = new TChain("event");
  //chain->Add("/afs/cern.ch/user/s/sani/work/globe_for_Shashlik/src/h2gglobe/AnalysisScripts/my_reduced_tuple_foo.root");
  chain->Add("/afs/cern.ch/user/s/sani/eos/cms/store/group/phys_higgs/cmshgg/processed/FixPFIsolationIssue/22AprFix/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_14TeV-pythia6_TP2023SHCALDR-SHCALMar26_PU140BX25_PH2_1K_FB_V6-v1/*.root");

  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("pho_n", 1);
  chain->SetBranchStatus("pho_p4", 1);
  if (recoEff) {
    chain->SetBranchStatus("gp_n", 1);
    chain->SetBranchStatus("gp_p4", 1);
    chain->SetBranchStatus("gp_mother", 1);
    chain->SetBranchStatus("gp_status", 1);
    chain->SetBranchStatus("gp_pdgid", 1);
    chain->SetBranchStatus("genjet_algo1_n",  1);
    chain->SetBranchStatus("genjet_algo1_p4", 1);
    chain->SetBranchStatus("genjet_algo1_em", 1);

    chain->SetBranchAddress("gp_n", &gp_n);
    chain->SetBranchAddress("gp_p4", &gp_p4);
    chain->SetBranchAddress("gp_mother", &gp_mother);
    chain->SetBranchAddress("gp_status", &gp_status);
    chain->SetBranchAddress("gp_pdgid", &gp_pdgid);
    
    chain->SetBranchAddress("genjet_algo1_n",  &genjet_algo1_n);
    chain->SetBranchAddress("genjet_algo1_p4", &genjet_algo1_p4);
    chain->SetBranchAddress("genjet_algo1_em", &genjet_algo1_em);    
  } else {
    chain->SetBranchStatus("dipho_gensel", 1);
    chain->SetBranchStatus("dipho_genfakeind", 1);
    chain->SetBranchStatus("dipho_genpromptind", 1);
    chain->SetBranchStatus("dipho_genfakee", 1);
    chain->SetBranchStatus("dipho_genprompte", 1);
    chain->SetBranchStatus("dipho_genfakeeta", 1);
    chain->SetBranchStatus("dipho_genprompteta", 1);
    chain->SetBranchStatus("dipho_genfakephi", 1);
    chain->SetBranchStatus("dipho_genpromptphi", 1);
    
    chain->SetBranchAddress("dipho_gensel", &dipho_gensel);
    chain->SetBranchAddress("dipho_genpromptind", &dipho_genpromptind);
    chain->SetBranchAddress("dipho_genfakeind", &dipho_genfakeind);
    chain->SetBranchAddress("dipho_genprompte", &dipho_genprompte);
    chain->SetBranchAddress("dipho_genfakee", &dipho_genfakee);
    chain->SetBranchAddress("dipho_genprompteta", &dipho_genprompteta);
    chain->SetBranchAddress("dipho_genfakeeta", &dipho_genfakeeta);
    chain->SetBranchAddress("dipho_genpromptphi", &dipho_genpromptphi);
    chain->SetBranchAddress("dipho_genfakephi", &dipho_genfakephi);
  }
 
  chain->SetBranchAddress("pho_n", &pho_n);
  chain->SetBranchAddress("pho_p4", &pho_p4);
  
  Int_t entries = chain->GetEntries();
  for (int z=0; z<entries; z++) {
    if (z%1000 == 0)
      std::cout << z << std::endl;
    chain->GetEntry(z);

    if (recoEff) {
      promptRecoEfficiency();
      fakeRecoEfficiency();
    } else {
      if (dipho_gensel) {
	nGenDiphoton();
	
	for (int ipho=0; ipho<pho_n; ipho++) {
	  nPromptID(ipho);
	  nFakeID(ipho);
	}
      }
    }
  }

  TString filename = "output_eff.root";
  if (recoEff)
    filename = "output_eff_reco.root";

  TFile* file = new TFile(filename, "recreate");
  
  if (recoEff) {
    recoEffPromptPt_den ->Write();
    recoEffPromptPt_num ->Write();
    recoEffFakePt_den   ->Write();
    recoEffFakePt_num   ->Write();  
    recoEffPromptEta_den ->Write();
    recoEffPromptEta_num ->Write();
    recoEffFakeEta_den   ->Write();
    recoEffFakeEta_num   ->Write();
  } else {
    nGenFakeablePt  ->Write();
    nGenFakeableEta ->Write();
    idFakeEta_num        ->Write();
    idFakePt_num        ->Write();
    idPromptPt_num      ->Write();
    idPromptEta_num      ->Write();
  }

  file->Close();
  
}
