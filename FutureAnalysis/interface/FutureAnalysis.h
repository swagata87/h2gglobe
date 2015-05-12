#ifndef __FUTUREANALYSIS__
#define __FUTUREANALYSIS__

#include "BaseAnalysis.h"
#include "RooContainer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class FutureAnalysis : public BaseAnalysis
{
public:
  
  FutureAnalysis();
  virtual ~FutureAnalysis();
  
  virtual const std::string & name() const { return name_; };
  
  // LoopAll analysis interface implementation
  void Init(LoopAll&);
  void Term(LoopAll&);
  
  virtual void FillReductionVariables(LoopAll& l, int jentry);   
  virtual bool SelectEventsReduction(LoopAll&, int);
  virtual void ReducedOutputTree(LoopAll &l, TTree * outputTree); 
  virtual void ResetAnalysis();
  virtual void GetBranches(TTree * outputTree, std::set<TBranch *>& );
  
  virtual bool SkimEvents(LoopAll&, int);
  virtual bool SelectEvents(LoopAll&, int);
  virtual bool Analysis(LoopAll&, Int_t);
  void FillFlatTree(LoopAll&, Int_t, Float_t, Int_t, Int_t);
  bool PhotonPreSelection(LoopAll&, Int_t, Int_t);

  bool makeOptTree;
  
protected:
  std::string name_;
};

#endif


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
