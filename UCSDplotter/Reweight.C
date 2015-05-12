#define Reweight_cxx
#include "Reweight.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>


void Reweight::Loop()
{
  gROOT->SetBatch();
  int counter_sieie_clean=0;
  int counter_sieie_general=0;

  float new_wt=0.0;

  TFile* outputFile = new TFile("2DReweighted_opttree_Ph2_Shashlik.root","RECREATE");
  TTree* myTree = new TTree("opttree_2Dwt","opttree_2Dwt");

  myTree->Branch("run", &run, "run/I");
  myTree->Branch("lumis", &lumis, "lumis/I");
  myTree->Branch("event", &event, "event/I");
  myTree->Branch("itype", &itype, "itype/I");
  myTree->Branch("nvtx", &nvtx, "nvtx/I");
  myTree->Branch("nPU", &nPU, "nPU/I");
  myTree->Branch("weight", &new_wt, "new_wt/F");
  myTree->Branch("rho", &rho, "rho/F");
  myTree->Branch("nPho", &nPho, "nPho/I");
  myTree->Branch("mva_legacy", &mva_legacy, "mva_legacy/F");

  myTree->Branch("mydipho_gensel", &mydipho_gensel, "mydipho_gensel/F");
  myTree->Branch("mydipho_genfakeind", &mydipho_genfakeind, "mydipho_genfakeind/F");
  myTree->Branch("mydipho_genpromptind", &mydipho_genpromptind, "mydipho_genpromptind/F");
  myTree->Branch("mydipho_genfakee", &mydipho_genfakee, "mydipho_genfakee/F");
  myTree->Branch("mydipho_genprompte", &mydipho_genprompte, "mydipho_genprompte/F");
  myTree->Branch("mydipho_genfakeeta", &mydipho_genfakeeta, "mydipho_genfakeeta/F");
  myTree->Branch("mydipho_genprompteta", &mydipho_genprompteta, "mydipho_genprompteta/F");
  myTree->Branch("mydipho_genfakephi", &mydipho_genfakephi, "mydipho_genfakephi/F");
  myTree->Branch("mydipho_genpromptphi", &mydipho_genpromptphi, "mydipho_genpromptphi/F");

  myTree->Branch("pT_pho", &pT_pho, "pT_pho/F");
  myTree->Branch("p4_eta_pho", &p4_eta_pho, "p4_eta_pho/F");
  myTree->Branch("sc_eta_pho", &sc_eta_pho, "sc_eta_pho/F");
  myTree->Branch("Et_pho", &Et_pho, "Et_pho/F");
  myTree->Branch("E2byE5_pho", &E2byE5_pho, "E2byE5_pho/F");
  myTree->Branch("E2byE5_cleaned_pho", &E2byE5_cleaned_pho, "E2byE5_cleaned_pho/F");
  myTree->Branch("r9_pho", &r9_pho, "r9_pho/F");
  myTree->Branch("r9_cleaned_pho", &r9_cleaned_pho, "r9_cleaned_pho/F");
  myTree->Branch("sieie_pho", &sieie_pho, "sieie_pho/F");
  myTree->Branch("sieie_cleaned_pho", &sieie_cleaned_pho, "sieie_cleaned_pho/F");
  myTree->Branch("hoe_pho", &hoe_pho, "hoe_pho/F");
  myTree->Branch("isEB_pho", &isEB_pho, "isEB_pho/F");
  myTree->Branch("trkiso03_pho", &trkiso03_pho, "trkiso03_pho/F");
  myTree->Branch("ecaliso03_pho", &ecaliso03_pho, "ecaliso03_pho/F");
  myTree->Branch("hcaliso03_pho", &hcaliso03_pho, "hcaliso03_pho/F");
  myTree->Branch("pf_pho_iso", &pf_pho_iso, "pf_pho_iso/F");
  myTree->Branch("pf_pho_iso_eta030", &pf_pho_iso_eta030, "pf_pho_iso_eta030/F");
  myTree->Branch("pf_pho_iso_eta045", &pf_pho_iso_eta045, "pf_pho_iso_eta045/F");
  myTree->Branch("pf_pho_iso_eta060", &pf_pho_iso_eta060, "pf_pho_iso_eta060/F");
  myTree->Branch("pf_pho_iso_eta075", &pf_pho_iso_eta075, "pf_pho_iso_eta075/F");
  myTree->Branch("pf_pho_iso_eta090", &pf_pho_iso_eta090, "pf_pho_iso_eta090/F");
  myTree->Branch("pf_pho_iso_dR070",  &pf_pho_iso_dR070,  "pf_pho_iso_dR070/F");
  myTree->Branch("pfPhoIso_cleaned_t15", &pfPhoIso_cleaned_t15, "pfPhoIso_cleaned_t15/F");
  myTree->Branch("pfPhoIso_cleaned_t10", &pfPhoIso_cleaned_t10, "pfPhoIso_cleaned_t10/F");
  myTree->Branch("ESEffSigmaRR_pho", &ESEffSigmaRR_pho, "ESEffSigmaRR_pho/F");
  myTree->Branch("genmatched_pho", &genmatched_pho, "genmatched_pho/I");
  myTree->Branch("etawidth_pho", &etawidth_pho, "etawidth_pho/F");
  myTree->Branch("phiwidth_pho", &phiwidth_pho, "phiwidth_pho/F");
  myTree->Branch("sieip_pho", &sieip_pho, "sieip_pho/F");
  myTree->Branch("sieip_cleaned_pho", &sieip_cleaned_pho, "sieip_cleaned_pho/F");
  myTree->Branch("scRaw_pho", &scRaw_pho, "scRaw_pho/F");
  myTree->Branch("hasPixelSeed_pho", &hasPixelSeed_pho, "hasPixelSeed_pho/I");
  myTree->Branch("pf_charged_iso_chosenvtx", &pf_charged_iso_chosenvtx, "pf_charged_iso_chosenvtx/F");
  myTree->Branch("pf_charged_iso_vtx1", &pf_charged_iso_vtx1, "pf_charged_iso_vtx1/F");
  myTree->Branch("pf_charged_iso_vtx2",&pf_charged_iso_vtx2, "pf_charged_iso_vtx2/F");
  myTree->Branch("pf_charged_iso_badvtx", &pf_charged_iso_badvtx, "pf_charged_iso_badvtx/F");
  myTree->Branch("sipip_pho", &sipip_pho, "sipip_pho/F");

  TH2F* hist2D_Pt_AbsEta_sig = new TH2F("pt_AbsEta_Sig", "pt_absEta_Sig", 250, 0.0, 500.0, 5, 0, 2.5);
  TH2F* hist2D_Pt_AbsEta_bkg = new TH2F("pt_AbsEta_Bkg", "pt_absEta_Bkg", 250, 0.0, 500.0, 5, 0, 2.5);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    float absEta=fabs(sc_eta_pho);
    
    if(genmatched_pho==1.0) {
      hist2D_Pt_AbsEta_sig->Fill(pT_pho,absEta,weight);
    }
    if(genmatched_pho==0.0) {
      hist2D_Pt_AbsEta_bkg->Fill(pT_pho,absEta,weight);
    }
  }
  
  hist2D_Pt_AbsEta_sig->Scale(1.0/hist2D_Pt_AbsEta_sig->Integral());
  hist2D_Pt_AbsEta_bkg->Scale(1.0/hist2D_Pt_AbsEta_bkg->Integral());
  
  TH2F ratio = ((*hist2D_Pt_AbsEta_sig)/(*hist2D_Pt_AbsEta_bkg));
  ratio.SetNameTitle("ratio_2D", "ratio");


  //////......///////
  for (Long64_t kentry=0; kentry<nentries;kentry++) {
    //   std::cout << "\nEvt " << kentry << std::endl;
    Long64_t ientry = LoadTree(kentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(kentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    float myAbsEta=fabs(sc_eta_pho);

    if(genmatched_pho==1.0) {
      new_wt=weight;
    }
    if(genmatched_pho==0.0) {
      Int_t binx_pt = ratio.GetXaxis()->FindBin(pT_pho);
      Int_t biny_absEta = ratio.GetYaxis()->FindBin(myAbsEta);
      Int_t bin = ratio.GetBin(binx_pt,biny_absEta,0);
      float mywt = ratio.GetBinContent(bin);
      if (mywt==0.0)  mywt=1.0;
      new_wt=weight*mywt;
    }
    if (sieie_cleaned_pho>=1.0)  counter_sieie_clean++;
    if (sieie_pho>=1.0)  counter_sieie_general++;
    
    if (sieie_cleaned_pho>=1.0)  continue;
    if (sieie_pho>=1.0)  continue;

    myTree->Fill();
  }  
  
  std::cout << "counter_sieie_clean=" << counter_sieie_clean << "  counter_sieie_general=" << counter_sieie_general << std::endl;  
  outputFile->cd();
  ratio.Write();
  myTree->Write();
  outputFile->Close();
  
}
