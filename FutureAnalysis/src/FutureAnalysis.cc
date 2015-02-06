#include "../interface/FutureAnalysis.h"
#include "Sorters.h"
#include "../../LoopAll.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
FutureAnalysis::FutureAnalysis()  
{
    name_ = "FutureAnalysis";
}

// ----------------------------------------------------------------------------------------------------
FutureAnalysis::~FutureAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::Term(LoopAll& l) 
{
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::Init(LoopAll& l) 
{
  std::cout << "Init -> FutureAnalysis "<< std::endl;
  l.SetAllMVA();
  if (bdtTrainingType == "")
    bdtTrainingType = bdtTrainingPhilosophy;

  // 2013 ID MVA                                                                                                                       
  if( photonLevel2013IDMVA_EB != "" && photonLevel2013IDMVA_EE != "" ) {
    l.tmvaReaderID_2013_Barrel->BookMVA("AdaBoost",photonLevel2013IDMVA_EB.c_str());
    l.tmvaReaderID_2013_Endcap->BookMVA("AdaBoost",photonLevel2013IDMVA_EE.c_str());
  } else if( photonLevel2012IDMVA_EB != "" && photonLevel2012IDMVA_EE != "" ) {
    l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevel2012IDMVA_EB.c_str());
    l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevel2012IDMVA_EE.c_str());
    assert( bdtTrainingType == "Moriond2013" );
  } else if (photonLevel2013_7TeV_IDMVA_EB != "" && photonLevel2013_7TeV_IDMVA_EE != "" ) {
    l.tmvaReaderID_2013_7TeV_MIT_Barrel->BookMVA("AdaBoost",photonLevel2013_7TeV_IDMVA_EB.c_str());
    l.tmvaReaderID_2013_7TeV_MIT_Endcap->BookMVA("AdaBoost",photonLevel2013_7TeV_IDMVA_EE.c_str());
  } else {
    assert( run7TeV4Xanalysis );
  }
}

void FutureAnalysis::ResetAnalysis(){}

void FutureAnalysis::GetBranches(TTree * outputTree, std::set<TBranch *>& ) 
{}


// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
  bool isCorrectVertex;
  TVector3* vtx;
  float etaStrip_run1 = 0.015;
  float etaStrip_030  = 0.030;
  float etaStrip_045  = 0.045;
  float etaStrip_060  = 0.060;
  float etaStrip_075  = 0.075;
  float etaStrip_090  = 0.090;

  float dRmax = 0.3;
  float dRVeto_run1 = 0.0;
  float dRVeto_070  = 0.070;

  int type = l.itype[l.current];
  float weight       = l.sampleContainer[l.current_sample_index].weight();
  //std::cout << "pho_n=" << l.pho_n << " ecalhit_n=" << l.ecalhit_n << std::endl;
  for (int ipho=0; ipho<l.pho_n; ipho++){
    //std::cout << "ipho=" << ipho << std::endl;
    TLorentzVector* p4_pho =  (TLorentzVector*) l.pho_p4->At(ipho);
    if (p4_pho->Pt() < 40.0) continue; 
    if (l.pho_genmatched==0) continue;
    float sc_eta_ipho =  ((TVector3*)l.sc_xyz->At(l.pho_scind[ipho]))->Eta();
    float phi_ipho = p4_pho->Phi();
    if (fabs(sc_eta_ipho) >1.479 ) continue;
    l.FillHist("pf_pho_iso_run1parameters", 0, l.pho_pfiso_myphoton03[ipho], weight);
    l.FillHist("pf_pho_iso_eta030", 0, l.pho_pfiso_myphoton03_eta030[ipho], weight);
    l.FillHist("pf_pho_iso_eta045", 0, l.pho_pfiso_myphoton03_eta045[ipho], weight);
    l.FillHist("pf_pho_iso_eta060", 0, l.pho_pfiso_myphoton03_eta060[ipho], weight);
    l.FillHist("pf_pho_iso_eta075", 0, l.pho_pfiso_myphoton03_eta075[ipho], weight);
    l.FillHist("pf_pho_iso_eta090", 0, l.pho_pfiso_myphoton03_eta090[ipho], weight);

    //    float iso=l.pho_pfiso_myphoton03[ipho];
    //  if (iso<60.0 || iso>100.0) continue;
    l.FillHist2D("my_pho_position", 0, sc_eta_ipho, phi_ipho, weight);
   
    //// PFcandidate loop ////
    for (int ipf=0; ipf<l.pfcand_n; ipf++) {
      if (l.pfcand_pdgid[ipf] != 4) continue;
      TLorentzVector* p4_pf = (TLorentzVector*) l.pfcand_p4->At(ipf);
      TVector3* pfvtx = (TVector3*)l.pfcand_posvtx->At(ipf);
      TVector3* phoEcalPos = (TVector3*)l.sc_xyz->At(l.pho_scind[ipho]);
      TVector3 photonDirectionWrtVtx = TVector3(phoEcalPos->X() - pfvtx->X(),
						phoEcalPos->Y() - pfvtx->Y(),
						phoEcalPos->Z() - pfvtx->Z());
      float dPhi = (photonDirectionWrtVtx.Phi() - p4_pf->Phi());
      ///////float dPhi = (p4_pho->Phi() - p4_pf->Phi());
      if(dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
      if(dPhi < (-1.0*(TMath::Pi()))) dPhi = (-1.0*(TMath::TwoPi())) - dPhi;
      ///   float deta = ((sc_eta_ipho - p4_pf->Eta()));
      float deta = (photonDirectionWrtVtx.Eta() - p4_pf->Eta());
      float dR = photonDirectionWrtVtx.DeltaR(p4_pf->Vect());
      if(dR > dRmax) continue ; 
      l.FillHist2D("pf_detadphi_noveto", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_run1) continue;
      if (dR < dRVeto_run1)           continue;
      l.FillHist2D("my_pfetaphi_Run1veto", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_030)  continue;
      l.FillHist2D("my_pfetaphi_deta030", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_045)  continue;
      l.FillHist2D("my_pfetaphi_deta045", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_060)  continue;
      l.FillHist2D("my_pfetaphi_deta060", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_075)  continue;
      l.FillHist2D("my_pfetaphi_deta075", 0, deta, dPhi, weight);

      if (fabs(deta) < etaStrip_090)  continue;
      l.FillHist2D("my_pfetaphi_deta090", 0, deta, dPhi, weight);

    }

    /*
    //// ecalhit loop ////
    for (int iecal=0; iecal<l.ecalhit_n; iecal++) {
      ////std::cout << "iecal=" << iecal << std::endl;
      TLorentzVector* p4_ecal = (TLorentzVector*) l.ecalhit_p4->At(iecal);
      float dPhi = (p4_pho->Phi() - p4_ecal->Phi());
      if(dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
      if(dPhi < (-1.0*(TMath::Pi()))) dPhi = (-1.0*(TMath::TwoPi())) - dPhi;
      float deta = ((sc_eta_ipho - p4_ecal->Eta()));
      l.FillHist2D("my_etaphi", 0, deta, dPhi, weight);

      //      float abs_dPhi = fabs(p4_pho->Phi() - p4_ecal->Phi());
      // if(abs_dPhi > TMath::Pi()) abs_dPhi = TMath::TwoPi() - dPhi;
      //  float abs_deta = fabs((sc_eta_ipho - p4_ecal->Eta()));

      //std::cout << "deta=" << deta << " dPhi=" << dPhi << std::endl;
      //  l.FillHist2D("my_absetaphi", 0, abs_deta, abs_dPhi, weight);
      
      }*/
  }
   
  for (int ipho=0;ipho<l.pho_n;ipho++){
    l.pho_s4ratio[ipho] = l.pho_e2x2[ipho]/l.bc_s25[l.sc_bcseedind[l.pho_scind[ipho]]];
    float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
    l.pho_ESEffSigmaRR[ipho] = 0.0;
    if(rr2>0. && rr2<999999.) {
      l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
    }
  }
  //  std::cout << "\n Evt----" << jentry << std::endl;                                                              
                
 
  //  int type = l.itype[l.current];
  //  float weight       = l.sampleContainer[l.current_sample_index].weight();
  float sampleweight = l.sampleContainer[l.current_sample_index].weight();

  int ivtx = l.dipho_vtxind[0];
  // std::cout << "Vertex Index = " << ivtx << std::endl;
  float nVtx = l.vtx_std_n ;
  vtx = (TVector3*)l.vtx_std_xyz->At(ivtx);
  isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
  int category = -1;
  if( isCorrectVertex ) { l.FillHist("my_nvtx_rv", category+1, nVtx, weight); }
  l.FillHist("my_nvtx", category+1, nVtx, weight);

  int lead = l.dipho_leadind[0];
  int sublead = l.dipho_subleadind[0];
  TLorentzVector* p4_pho_lead = (TLorentzVector*) l.pho_p4->At(lead);
  TLorentzVector* p4_pho_sublead = (TLorentzVector*) l.pho_p4->At(sublead);

  float phoIDMVA_lead = l.photonIDMVA(lead, ivtx, *p4_pho_lead, bdtTrainingType.c_str() );
  float phoIDMVA_sublead = l.photonIDMVA(sublead, ivtx, *p4_pho_sublead, bdtTrainingType.c_str() );

  TLorentzVector* ph_lead = dynamic_cast<TLorentzVector*>(l.pho_p4->At(lead));
  float pt_lead = ph_lead->Pt();
  float AbsEta_lead = fabs(ph_lead->Eta());

  TLorentzVector* ph_sublead = dynamic_cast<TLorentzVector*>(l.pho_p4->At(sublead));
  float pt_sublead = ph_sublead->Pt();
  float AbsEta_sublead = fabs(ph_sublead->Eta());
  ///  if ( (l.pho_isEB[lead])&&(l.pho_isEB[sublead]) || (!(l.pho_isEB[lead]))&&(!(l.pho_isEB[sublead])) ) {
  // if (pt_lead>33.0 && pt_sublead>25.0) {
  if (pt_lead>15.0 && pt_sublead>15.0) {  
    if ( (PhotonPreSelection(l,lead,ivtx)==true) && (PhotonPreSelection(l,sublead,ivtx)==true) ) {
      if (makeOptTree) {
	//// std::cout << "\nLeadAbseta->" << AbsEta_lead << " LeadIsEB->"  << l.pho_isEB[lead] << std::endl;                      
	//// std::cout << "SubLeadAbseta->" << AbsEta_sublead << " SubLeadIsEB->"  << l.pho_isEB[sublead] << std::endl;            
	FillFlatTree(l, type, weight, lead,    ivtx, phoIDMVA_lead);
	FillFlatTree(l, type, weight, sublead, ivtx, phoIDMVA_sublead);
      }
      //    }
    }
  }
  return false;
}

// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SelectEvents(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}
   
// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SelectEventsReduction(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::SkimEvents(LoopAll&, int)
{
  return true;
}

// ----------------------------------------------------------------------------------------------------
void FutureAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
}

void FutureAnalysis::FillFlatTree(LoopAll& l, Int_t type, Float_t weight, Int_t pho_indx, Int_t vtx, Float_t mva) {
  l.FillTree("run",        l.run,         "Future_myTry");
  l.FillTree("lumis",      l.lumis,       "Future_myTry");
  l.FillTree("event",      l.event,       "Future_myTry");
  l.FillTree("weight",     weight,        "Future_myTry");
  l.FillTree("itype",      type,          "Future_myTry");
  l.FillTree("nvtx",       l.vtx_std_n,   "Future_myTry");
  l.FillTree("nPU",        l.pu_n,        "Future_myTry");
  l.FillTree("rho",        l.rho_algo1,   "Future_myTry");
  l.FillTree("nPho",       l.pho_n,       "Future_myTry");

  float pho_E2byE5 = 0.0;
  float pho_E2byE5_cleaned = 0.0;
  if (pho_indx != -1) {
    TLorentzVector* p4_photon = (TLorentzVector*) l.pho_p4->At(pho_indx);
    //std::cout << "E2x2=" << l.pho_e2x2[pho_indx] << "  E5x5=" << l.pho_e5x5[pho_indx] << std::endl;                                
    pho_E2byE5 = ((l.pho_e2x2[pho_indx])/(l.pho_e5x5[pho_indx]));
      //l.pho_e2x2[pho_indx]/l.bc_s25[l.sc_bcseedind[l.pho_scind[pho_indx]]];
    //    std::cout << "\n\nsc_eta_pho=" << ((TVector3*)l.sc_xyz->At(l.pho_scind[pho_indx]))->Eta() << "  pho_E2byE5=" << pho_E2byE5 << "  pho_e2x2=" << l.pho_e2x2[pho_indx] << "  bc_s25=" << l.bc_s25[l.sc_bcseedind[l.pho_scind[pho_indx]]] << "  sc_bcseedind=" << l.sc_bcseedind[l.pho_scind[pho_indx]]  << std::endl;
    //l.pho_e2x2[pho_indx]/l.pho_e5x5[pho_indx];
    pho_E2byE5_cleaned = l.pho_e2x2_cleaned[pho_indx]/l.pho_e5x5_cleaned[pho_indx];
      //l.pho_e2x2_cleaned[pho_indx]/l.bc_s25[l.sc_bcseedind[l.pho_scind[pho_indx]]];
    ///l.pho_e2x2_cleaned[pho_indx]/l.pho_e5x5_cleaned[pho_indx];
    
    if (isnan(pho_E2byE5)) pho_E2byE5=0.0;
    if (isnan(pho_E2byE5_cleaned)) pho_E2byE5_cleaned=0.0;
    if (isinf(pho_E2byE5)) pho_E2byE5=0.0;
    if (isinf(pho_E2byE5_cleaned)) pho_E2byE5_cleaned=0.0;


    if ( isnan(pho_E2byE5) ) std::cout << "NAN ******* pho_E2byE5 " << pho_E2byE5 << std::endl;
    if ( isnan(pho_E2byE5_cleaned) ) std::cout << "NAN *******pho_E2byE5_cleaned "<< pho_E2byE5_cleaned << std::endl;

    if ( isinf(pho_E2byE5) ) std::cout << "INF ********* pho_E2byE5 " << pho_E2byE5 << std::endl;
    if ( isinf(pho_E2byE5_cleaned) ) std::cout << "INF ********* pho_E2byE5_cleaned "<< pho_E2byE5_cleaned << std::endl;

    float pf_charged_iso03_chosenvtx =  (*l.pho_pfiso_mycharged03)[pho_indx][vtx];
    
    int ivtx_1 = l.dipho_vtxind[1];
    float pf_charged_iso03_vtx1 =  (*l.pho_pfiso_mycharged03)[pho_indx][ivtx_1];

    int ivtx_2 = l.dipho_vtxind[2];
    float pf_charged_iso03_vtx2 =  (*l.pho_pfiso_mycharged03)[pho_indx][ivtx_2];

    float val_Pf_charged03_iso_badVtx = -999.0;
    for(int iv=0; iv < l.vtx_std_n; iv++) {
      if((*l.pho_pfiso_mycharged03)[pho_indx][iv] > val_Pf_charged03_iso_badVtx) {
        val_Pf_charged03_iso_badVtx = (*l.pho_pfiso_mycharged03)[pho_indx][iv];
      }
    }
    //    std::cout << "pho_E2byE5=" << pho_E2byE5 << std::endl;                                                                     
    float pfPhoIso_cleaned_t15 = pfEcalIso(l, pho_indx, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, 4, 15.0);
    float pfPhoIso_cleaned_t10 = pfEcalIso(l, pho_indx, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, 4, 10.0);

    TVector3* phoEcalPos = (TVector3*)l.sc_xyz->At(l.pho_scind[pho_indx]);
    float phoEcalPosX = phoEcalPos->X() ;
    float phoEcalPosY = phoEcalPos->Y() ;
    float phoEcalPosZ = phoEcalPos->Z() ; 
    //    std::cout << "phoEcalPosX = " << phoEcalPosX << std::endl;

    //    std::cout << "mva = " << mva << std::endl;                       
    l.FillTree("mva_legacy",               mva,                                                      "Future_myTry");
    l.FillTree("pho_EcalPosX",             phoEcalPosX,                                              "Future_myTry");
    l.FillTree("pho_EcalPosY",             phoEcalPosY,                                              "Future_myTry");
    l.FillTree("pho_EcalPosZ",             phoEcalPosZ,                                              "Future_myTry");
    l.FillTree("pT_pho",                   p4_photon->Pt(),                                          "Future_myTry");
    l.FillTree("p4_eta_pho",               p4_photon->Eta(),                                         "Future_myTry");
    l.FillTree("sc_eta_pho",               ((TVector3*)l.sc_xyz->At(l.pho_scind[pho_indx]))->Eta(),  "Future_myTry" );
    l.FillTree("Et_pho",                   p4_photon->Et(),                                          "Future_myTry");
    l.FillTree("E2byE5_pho",               pho_E2byE5,                                               "Future_myTry");
    l.FillTree("E2byE5_cleaned_pho",       pho_E2byE5_cleaned,                                       "Future_myTry");
    //   if (l.pho_r9[pho_indx]<0.5) std::cout << "r9<0.5" << std::endl;
    l.FillTree("r9_pho",                   l.pho_r9[pho_indx],                                       "Future_myTry");
    l.FillTree("r9_cleaned_pho",           l.pho_r9_cleaned[pho_indx],                               "Future_myTry");

    l.FillTree("sieie_pho",                l.pho_sieie[pho_indx],                                    "Future_myTry");
    l.FillTree("sieie_cleaned_pho",        l.pho_sieie_cleaned[pho_indx],                            "Future_myTry");
    l.FillTree("sipip_pho",                l.pho_sipip[pho_indx],                                    "Future_myTry");
    l.FillTree("hoe_pho",                  l.pho_hoe[pho_indx],                                      "Future_myTry");
    l.FillTree("isEB_pho",                 l.pho_isEB[pho_indx],                                     "Future_myTry");
    l.FillTree("trkiso03_pho",             l.pho_trksumpthollowconedr03[pho_indx],                   "Future_myTry");
    l.FillTree("ecaliso03_pho",            l.pho_ecalsumetconedr03[pho_indx],                        "Future_myTry");
    l.FillTree("hcaliso03_pho",            l.pho_hcalsumetconedr03[pho_indx],                        "Future_myTry");
    l.FillTree("pf_pho_iso",               l.pho_pfiso_myphoton03[pho_indx],                         "Future_myTry");
    l.FillTree("pf_pho_iso_eta030",        l.pho_pfiso_myphoton03_eta030[pho_indx],                  "Future_myTry");
    l.FillTree("pf_pho_iso_eta045",        l.pho_pfiso_myphoton03_eta045[pho_indx],                  "Future_myTry");
    l.FillTree("pf_pho_iso_eta060",        l.pho_pfiso_myphoton03_eta060[pho_indx],                  "Future_myTry");
    l.FillTree("pf_pho_iso_eta075",        l.pho_pfiso_myphoton03_eta075[pho_indx],                  "Future_myTry");
    l.FillTree("pf_pho_iso_eta090",        l.pho_pfiso_myphoton03_eta090[pho_indx],                  "Future_myTry");
    l.FillTree("pf_pho_iso_dR070",         l.pho_pfiso_myphoton03_dR070[pho_indx],                   "Future_myTry");
    // if (pf_charged_iso03_chosenvtx>14.0) std::cout<< "Pf_chIso03_chosenvtx>14.0" << std::endl; 
    l.FillTree("pf_charged_iso_chosenvtx", pf_charged_iso03_chosenvtx,                               "Future_myTry");
    l.FillTree("pf_charged_iso_vtx1",      pf_charged_iso03_vtx1,                                    "Future_myTry");
    l.FillTree("pf_charged_iso_vtx2",      pf_charged_iso03_vtx2,                                    "Future_myTry");
    l.FillTree("pf_charged_iso_badvtx",    val_Pf_charged03_iso_badVtx,                              "Future_myTry");
    l.FillTree("genmatched_pho",           (int)l.pho_genmatched[pho_indx],                          "Future_myTry");
    l.FillTree("etawidth_pho",             l.pho_etawidth[pho_indx],                                 "Future_myTry");
    l.FillTree("phiwidth_pho",             l.sc_sphi[l.pho_scind[pho_indx]],                         "Future_myTry");
    l.FillTree("sieip_pho",                l.pho_sieip[pho_indx],                                    "Future_myTry");
    l.FillTree("sieip_cleaned_pho",        l.pho_sieip_cleaned[pho_indx],                            "Future_myTry");
    l.FillTree("scRaw_pho",                l.sc_raw[l.pho_scind[pho_indx]],                          "Future_myTry");
    l.FillTree("hasPixelSeed_pho",         l.pho_haspixseed[pho_indx],                               "Future_myTry");
    l.FillTree("pfPhoIso_cleaned_t15",     pfPhoIso_cleaned_t15,                                     "Future_myTry");
    l.FillTree("pfPhoIso_cleaned_t10",     pfPhoIso_cleaned_t10,                                     "Future_myTry");
  }

  else {
    l.FillTree("mva_legacy",                9999.0,       "Future_myTry");
    l.FillTree("pho_EcalPosX",              9999.0,       "Future_myTry");
    l.FillTree("pho_EcalPosY",              9999.0,       "Future_myTry");
    l.FillTree("pho_EcalPosZ",              9999.0,       "Future_myTry");
    l.FillTree("pT_pho",                    9999.0,       "Future_myTry");
    l.FillTree("p4_eta_pho",                9999.0,       "Future_myTry");
    l.FillTree("sc_eta_pho",                9999.0,       "Future_myTry");
    l.FillTree("Et_pho",                    9999.0,       "Future_myTry");
    l.FillTree("E2byE5_pho",                9999.0,       "Future_myTry");
    l.FillTree("E2byE5_cleaned_pho",        9999.0,       "Future_myTry");
    l.FillTree("r9_pho",                    9999.0,       "Future_myTry");
    l.FillTree("r9_cleaned_pho",            9999.0,       "Future_myTry");
    l.FillTree("sieie_pho",                 9999.0,       "Future_myTry");
    l.FillTree("sieie_cleaned_pho",         9999.0,       "Future_myTry");
    l.FillTree("sipip_pho",                 9999.0,       "Future_myTry");
    l.FillTree("hoe_pho",                   9999.0,       "Future_myTry");
    l.FillTree("isEB_pho",                  9999.0,       "Future_myTry");
    l.FillTree("trkiso03_pho",              9999.0,       "Future_myTry");
    l.FillTree("ecaliso03_pho",             9999.0,       "Future_myTry");
    l.FillTree("hcaliso03_pho",             9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso",                9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_eta030",         9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_eta045",         9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_eta060",         9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_eta075",         9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_eta090",         9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso_dR070",          9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_chosenvtx",  9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_vtx1",       9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_vtx2",       9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_badvtx",     9999.0,       "Future_myTry");
    l.FillTree("genmatched_pho",            9999.0,       "Future_myTry");
    l.FillTree("etawidth_pho",              9999.0,       "Future_myTry");
    l.FillTree("phiwidth_pho",              9999.0,       "Future_myTry");
    l.FillTree("sieip_pho",                 9999.0,       "Future_myTry");
    l.FillTree("sieip_cleaned_pho",         9999.0,       "Future_myTry");
    l.FillTree("scRaw_pho",                 9999.0,       "Future_myTry");
    l.FillTree("hasPixelSeed_pho",          9999.0,       "Future_myTry");
    l.FillTree("pfPhoIso_cleaned_t15",      9999.0,       "Future_myTry");
    l.FillTree("pfPhoIso_cleaned_t10",      9999.0,       "Future_myTry");

  }

  l.FillTreeContainer();
}

float FutureAnalysis::pfEcalIso(LoopAll& l, int phoindex, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel,
				float etaStripEndcap, float thrBarrel, float thrEndcaps, int pfToUse, float time_cut) {

  float dRVeto, etaStrip, thr;
  if (l.pho_isEB[phoindex]) {
    dRVeto = dRVetoBarrel;
    etaStrip = etaStripBarrel;
    thr = thrBarrel;
  } else {
    dRVeto = dRVetoEndcap;
    etaStrip = etaStripEndcap;
    thr = thrEndcaps;
  }
  float sum = 0;
  for(unsigned i=0; i<l.pfcand_n; i++) {

    if (l.pfcand_pdgid[i] == pfToUse) {

      TVector3* pfvtx = (TVector3*)l.pfcand_posvtx->At(i);
      TVector3* phoEcalPos = (TVector3*)l.sc_xyz->At(l.pho_scind[phoindex]);

      TVector3 photonDirectionWrtVtx = TVector3(phoEcalPos->X() - pfvtx->X(),
                                                phoEcalPos->Y() - pfvtx->Y(),
                                                phoEcalPos->Z() - pfvtx->Z());

      TLorentzVector* pfc = (TLorentzVector*)l.pfcand_p4->At(i);

      if( pfc->Pt() < thr )
        continue;

      float dEta = fabs(photonDirectionWrtVtx.Eta() - pfc->Eta());
      float dR = photonDirectionWrtVtx.DeltaR(pfc->Vect());

      if (dEta < etaStrip) continue;
      if(dR > dRmax || dR < dRVeto) continue;
      if ((fabs(l.pfcand_time[i])) > time_cut) continue;

      sum += pfc->Pt();
    }
  }
  return sum;
}

bool FutureAnalysis::PhotonPreSelection(LoopAll &l, Int_t photon_index, Int_t vtx) {
  //  std::cout << "inside PhotonPreSelection() " << std::endl;                                                                      
  TLorentzVector* ph = dynamic_cast<TLorentzVector*>(l.pho_p4->At(photon_index));
  float pt = ph->Pt();
  float fabs_eta = fabs(ph->Eta());
  if (l.pho_sieie[photon_index] > 1.0) return false;
  if (l.pho_sieie_cleaned[photon_index] > 1.0) return false;

  //// **** AD HOC PRESELECTION **** ////
  ///// if (l.pho_r9[photon_index] <= 0.5) return false;
  ////  float charged_iso03_chosen =  (*l.pho_pfiso_mycharged03)[photon_index][vtx];
  //// if (charged_iso03_chosen >= 14.0) return false;

  if(fabs_eta>2.5) return false;
  if ( (fabs_eta>1.4442) && (fabs_eta<1.566) ) return false;

  int r9_category = (int) (l.pho_r9[photon_index] <= 0.9);
  int photon_category = r9_category + 2*l.PhotonEtaCategory(photon_index,2);
 
  ///// MIT Preselection //////
  float mitCuts_hoe[4]                 = {0.082,0.075,0.075,0.075};
  float mitCuts_sieie[4]               = {0.014,0.014,0.034,0.034};
  float mitCuts_ecaliso[4]             = {50,4,50,4};
  float mitCuts_hcaliso[4]             = {50,4,50,4};
  float mitCuts_trkiso[4]              = {50,4,50,4};
  float mitCuts_pfiso[4]               = {4,4,4,4}; // WARN if depends on category should change below

  float val_hoe        = l.pho_hoe[photon_index];
  float val_sieie      = l.pho_sieie[photon_index];
  float val_ecaliso = l.pho_ecalsumetconedr03[photon_index] - 0.012*ph->Et();
  float val_hcaliso = l.pho_hcalsumetconedr03[photon_index] - 0.005*ph->Et();
  float val_trkiso  = l.pho_trksumpthollowconedr03[photon_index] - 0.002*ph->Et(); 
  int val_pho_isconv = l.pho_isconv[photon_index];
  float val_pfiso02 = (*l.pho_pfiso_mycharged02)[photon_index][vtx];
  
  //  if (val_hoe             >= mitCuts_hoe[photon_category]         ) return false;
  //  if (val_sieie           >= mitCuts_sieie[photon_category]       ) return false;
  //  if (val_hcaliso         >= mitCuts_hcaliso[photon_category]     ) return false;
  //  if (val_trkiso          >= mitCuts_trkiso[photon_category]      ) return false;
  //  if ((!val_pho_isconv)) return false; // Electron Rejection based Conversion Safe Veto 
  //  if ( val_pfiso02 >= mitCuts_pfiso[photon_category]) return false;
  
  
  return true;
}


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
