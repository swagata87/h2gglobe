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

}

void FutureAnalysis::ResetAnalysis(){}

void FutureAnalysis::GetBranches(TTree * outputTree, std::set<TBranch *>& ) 
{}


// ----------------------------------------------------------------------------------------------------
bool FutureAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
  bool isCorrectVertex;
  TVector3* vtx;

  int type = l.itype[l.current];
  float weight       = l.sampleContainer[l.current_sample_index].weight();
  ////std::cout << "weight = " << weight << std::endl;
  for (int ipho=0;ipho<l.pho_n;ipho++){
    l.pho_s4ratio[ipho] = l.pho_e2x2[ipho]/l.bc_s25[l.sc_bcseedind[l.pho_scind[ipho]]];
    float rr2=l.pho_eseffsixix[ipho]*l.pho_eseffsixix[ipho]+l.pho_eseffsiyiy[ipho]*l.pho_eseffsiyiy[ipho];
    l.pho_ESEffSigmaRR[ipho] = 0.0;
    if(rr2>0. && rr2<999999.) {
      l.pho_ESEffSigmaRR[ipho] = sqrt(rr2);
    }
  }
  float sampleweight = l.sampleContainer[l.current_sample_index].weight();
  int ivtx = l.dipho_vtxind[0];
  // std::cout << "Vertex Index = " << ivtx << std::endl;
  float nVtx = l.vtx_std_n ;
  //// std::cout << "dipho_n=" << l.dipho_n << std::endl;
  vtx = (TVector3*)l.vtx_std_xyz->At(ivtx);
  isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
  int category = -1;
  if( isCorrectVertex ) { l.FillHist("my_nvtx_rv", category+1, nVtx, weight); }
  l.FillHist("my_nvtx", category+1, nVtx, weight);

  int lead = l.dipho_leadind[0];
  int sublead = l.dipho_subleadind[0];
  //// std::cout << "lead=" << lead << " sublead=" << sublead << std::endl;

  //  TLorentzVector* p4_pho_lead = (TLorentzVector*) l.pho_p4->At(lead);
  //  TLorentzVector* p4_pho_sublead = (TLorentzVector*) l.pho_p4->At(sublead);

  TLorentzVector* ph_lead = dynamic_cast<TLorentzVector*>(l.pho_p4->At(lead));
  ///std::cout << "Fine upto this, could access l.pho_p4->At(lead) " << std::endl;
  float pt_lead = ph_lead->Pt();
  float AbsEta_lead = fabs(ph_lead->Eta());
  ////std::cout << "Going to access l.pho_p4->At(sublead) " << std::endl;

  TLorentzVector* ph_sublead = dynamic_cast<TLorentzVector*>(l.pho_p4->At(sublead));
  //// std::cout << "Could access l.pho_p4->At(sublead) " << std::endl;

  float pt_sublead = ph_sublead->Pt();
  float AbsEta_sublead = fabs(ph_sublead->Eta());

  if (PADEBUG) std::cout << "Before pT cut and preselection cut" << std::endl;
  if (pt_lead>25.0 && pt_sublead>25.0) {  
    if ( (PhotonPreSelection(l,lead,ivtx)==true) && (PhotonPreSelection(l,sublead,ivtx)==true) ) {
      if (makeOptTree) {
	FillFlatTree(l, type, weight, lead,    ivtx);
	FillFlatTree(l, type, weight, sublead, ivtx);
      }
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

void FutureAnalysis::FillFlatTree(LoopAll& l, Int_t type, Float_t weight, Int_t pho_indx, Int_t vtx) {
  l.FillTree("run",                l.run,                  "Future_myTry");
  l.FillTree("lumis",              l.lumis,                "Future_myTry");
  l.FillTree("event",              l.event,                "Future_myTry");
  l.FillTree("weight",             weight,                 "Future_myTry");
  l.FillTree("itype",              type,                   "Future_myTry");
  l.FillTree("nvtx",               l.vtx_std_n,            "Future_myTry");
  l.FillTree("nPU",                l.pu_n,                 "Future_myTry");
  l.FillTree("rho",                l.rho_algo1,            "Future_myTry");
  l.FillTree("nPho",               l.pho_n,                "Future_myTry");
  l.FillTree("mydipho_gensel",       l.dipho_gensel,         "Future_myTry");
  l.FillTree("mydipho_genfakeind",   l.dipho_genfakeind,     "Future_myTry");
  l.FillTree("mydipho_genpromptind", l.dipho_genpromptind,   "Future_myTry");
  l.FillTree("mydipho_genfakee",     l.dipho_genfakee,       "Future_myTry");
  l.FillTree("mydipho_genprompte",   l.dipho_genprompte,     "Future_myTry");
  l.FillTree("mydipho_genfakeeta",   l.dipho_genfakeeta,     "Future_myTry");
  l.FillTree("mydipho_genprompteta", l.dipho_genprompteta,   "Future_myTry");
  l.FillTree("mydipho_genfakephi",   l.dipho_genfakephi,     "Future_myTry");
  l.FillTree("mydipho_genpromptphi", l.dipho_genpromptphi,   "Future_myTry");
  
  
  if (PADEBUG) std::cout  << "Global evt variables filled" << std::endl;
  float pho_E2byE5 = 0.0;
  if (pho_indx != -1) {
    TLorentzVector* p4_photon = (TLorentzVector*) l.pho_p4->At(pho_indx);
    pho_E2byE5 = ((l.pho_e2x2[pho_indx])/(l.pho_e5x5[pho_indx]));

    if (isnan(pho_E2byE5)) pho_E2byE5=0.0;
    if (isinf(pho_E2byE5)) pho_E2byE5=0.0;

    if ( isnan(pho_E2byE5) ) std::cout << "NAN ******* pho_E2byE5 " << pho_E2byE5 << std::endl;
    if ( isinf(pho_E2byE5) ) std::cout << "INF ********* pho_E2byE5 " << pho_E2byE5 << std::endl;
    if (PADEBUG) std::cout  << "Will access pho_pfiso_mycharged03 branch " << pho_indx << " " << vtx  << std::endl;
    float pf_charged_iso03_chosenvtx =  (*l.pho_pfiso_mycharged03)[pho_indx][vtx];
    if (PADEBUG) std::cout  << "Have accessed pho_pfiso_mycharged03 branch" << std::endl;
    if (PADEBUG) std::cout << "pf_charged_iso03_chosenvtx " << pf_charged_iso03_chosenvtx << std::endl;
    int ivtx_1 = l.dipho_vtxind[1];
    float pf_charged_iso03_vtx1 =  (*l.pho_pfiso_mycharged03)[pho_indx][ivtx_1];

    int ivtx_2 = l.dipho_vtxind[2];
    float pf_charged_iso03_vtx2 =  (*l.pho_pfiso_mycharged03)[pho_indx][ivtx_2];
    
    float val_Pf_charged03_iso_badVtx = -999.0;
    if (PADEBUG) std::cout  << "Will find val_Pf_charged03_iso_badVtx" << std::endl;
    for(int iv=0; iv < l.vtx_std_n; iv++) {
      if((*l.pho_pfiso_mycharged03)[pho_indx][iv] > val_Pf_charged03_iso_badVtx) {
        val_Pf_charged03_iso_badVtx = (*l.pho_pfiso_mycharged03)[pho_indx][iv];
      }
    }
    //    std::cout << "pho_E2byE5=" << pho_E2byE5 << std::endl;                                                                     
    
    TVector3* phoEcalPos = (TVector3*)l.sc_xyz->At(l.pho_scind[pho_indx]);
    float phoEcalPosX = phoEcalPos->X() ;
    float phoEcalPosY = phoEcalPos->Y() ;
    float phoEcalPosZ = phoEcalPos->Z() ; 
    //    std::cout << "phoEcalPosX = " << phoEcalPosX << std::endl;
    if (PADEBUG) std::cout  << "Will fill photon variables" << std::endl;

    l.FillTree("pho_EcalPosX",             phoEcalPosX,                                              "Future_myTry");
    l.FillTree("pho_EcalPosY",             phoEcalPosY,                                              "Future_myTry");
    l.FillTree("pho_EcalPosZ",             phoEcalPosZ,                                              "Future_myTry");
    l.FillTree("pT_pho",                   p4_photon->Pt(),                                          "Future_myTry");
    l.FillTree("p4_eta_pho",               p4_photon->Eta(),                                         "Future_myTry");
    l.FillTree("sc_eta_pho",               ((TVector3*)l.sc_xyz->At(l.pho_scind[pho_indx]))->Eta(),  "Future_myTry" );
    l.FillTree("Et_pho",                   p4_photon->Et(),                                          "Future_myTry");
    l.FillTree("E2byE5_pho",               pho_E2byE5,                                               "Future_myTry");
    l.FillTree("r9_pho",                   l.pho_r9[pho_indx],                                       "Future_myTry");
    l.FillTree("sieie_pho",                l.pho_sieie[pho_indx],                                    "Future_myTry");
    l.FillTree("sipip_pho",                l.pho_sipip[pho_indx],                                    "Future_myTry");
    l.FillTree("hoe_pho",                  l.pho_hoe[pho_indx],                                      "Future_myTry");
    l.FillTree("isEB_pho",                 l.pho_isEB[pho_indx],                                     "Future_myTry");
    l.FillTree("trkiso03_pho",             l.pho_trksumpthollowconedr03[pho_indx],                   "Future_myTry");
    l.FillTree("ecaliso03_pho",            l.pho_ecalsumetconedr03[pho_indx],                        "Future_myTry");
    l.FillTree("hcaliso03_pho",            l.pho_hcalsumetconedr03[pho_indx],                        "Future_myTry");
    l.FillTree("pf_pho_iso",               l.pho_pfiso_myphoton03[pho_indx],                         "Future_myTry");
    if (PADEBUG) std::cout << "Will fill pf_charged_iso_chosenvtx" << std::endl;
    l.FillTree("pf_charged_iso_chosenvtx", pf_charged_iso03_chosenvtx,                               "Future_myTry");
    if (PADEBUG) std::cout << "Will fill pf_charged_iso_vtx1" << std::endl;
    l.FillTree("pf_charged_iso_vtx1",      pf_charged_iso03_vtx1,                                    "Future_myTry");
    if (PADEBUG) std::cout << "Will fill pf_charged_iso_vtx2" << std::endl;
    l.FillTree("pf_charged_iso_vtx2",      pf_charged_iso03_vtx2,                                    "Future_myTry");
    if (PADEBUG) std::cout << "Will fill pf_charged_iso_badvtx" << std::endl; 
    l.FillTree("pf_charged_iso_badvtx",    val_Pf_charged03_iso_badVtx,                              "Future_myTry");
    l.FillTree("genmatched_pho",           (int)l.pho_genmatched[pho_indx],                          "Future_myTry");
    l.FillTree("etawidth_pho",             l.pho_etawidth[pho_indx],                                 "Future_myTry");
    l.FillTree("phiwidth_pho",             l.sc_sphi[l.pho_scind[pho_indx]],                         "Future_myTry");
    l.FillTree("sieip_pho",                l.pho_sieip[pho_indx],                                    "Future_myTry");
    l.FillTree("scRaw_pho",                l.sc_raw[l.pho_scind[pho_indx]],                          "Future_myTry");
    l.FillTree("hasPixelSeed_pho",         l.pho_haspixseed[pho_indx],                               "Future_myTry");
    float rr2=l.pho_eseffsixix[pho_indx]*l.pho_eseffsixix[pho_indx]+l.pho_eseffsiyiy[pho_indx]*l.pho_eseffsiyiy[pho_indx];
    float pho_ESEffSigmaRR = 0.0;
    if(rr2>0. && rr2<999999.) {
      pho_ESEffSigmaRR = sqrt(rr2);
    }
    l.FillTree("ESEffSigmaRR_pho",         pho_ESEffSigmaRR,                               "Future_myTry");
  }

  else {
    l.FillTree("pho_EcalPosX",              9999.0,       "Future_myTry");
    l.FillTree("pho_EcalPosY",              9999.0,       "Future_myTry");
    l.FillTree("pho_EcalPosZ",              9999.0,       "Future_myTry");
    l.FillTree("pT_pho",                    9999.0,       "Future_myTry");
    l.FillTree("p4_eta_pho",                9999.0,       "Future_myTry");
    l.FillTree("sc_eta_pho",                9999.0,       "Future_myTry");
    l.FillTree("Et_pho",                    9999.0,       "Future_myTry");
    l.FillTree("E2byE5_pho",                9999.0,       "Future_myTry");
    l.FillTree("r9_pho",                    9999.0,       "Future_myTry");
    l.FillTree("sieie_pho",                 9999.0,       "Future_myTry");
    l.FillTree("sipip_pho",                 9999.0,       "Future_myTry");
    l.FillTree("hoe_pho",                   9999.0,       "Future_myTry");
    l.FillTree("isEB_pho",                  9999.0,       "Future_myTry");
    l.FillTree("trkiso03_pho",              9999.0,       "Future_myTry");
    l.FillTree("ecaliso03_pho",             9999.0,       "Future_myTry");
    l.FillTree("hcaliso03_pho",             9999.0,       "Future_myTry");
    l.FillTree("pf_pho_iso",                9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_chosenvtx",  9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_vtx1",       9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_vtx2",       9999.0,       "Future_myTry");
    l.FillTree("pf_charged_iso_badvtx",     9999.0,       "Future_myTry");
    l.FillTree("genmatched_pho",            9999.0,       "Future_myTry");
    l.FillTree("etawidth_pho",              9999.0,       "Future_myTry");
    l.FillTree("phiwidth_pho",              9999.0,       "Future_myTry");
    l.FillTree("sieip_pho",                 9999.0,       "Future_myTry");
    l.FillTree("scRaw_pho",                 9999.0,       "Future_myTry");
    l.FillTree("hasPixelSeed_pho",          9999.0,       "Future_myTry");
    l.FillTree("ESEffSigmaRR_pho",          9999.0,       "Future_myTry");
  }

  l.FillTreeContainer();
}


bool FutureAnalysis::PhotonPreSelection(LoopAll &l, Int_t photon_index, Int_t vtx) {
  //  std::cout << "inside PhotonPreSelection() " << std::endl;                                                                      
  TLorentzVector* ph = dynamic_cast<TLorentzVector*>(l.pho_p4->At(photon_index));
  float fabs_eta = fabs(ph->Eta());
  if (l.pho_sieie[photon_index] > 1.0) return false;
  /// if (fabs_eta>2.5) std::cout << "fabs_eta=" << fabs_eta << std::endl;
  if(fabs_eta>3.0) return false;
  if ( (fabs_eta>1.4442) && (fabs_eta<1.566) ) return false;

  return true;
}


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
