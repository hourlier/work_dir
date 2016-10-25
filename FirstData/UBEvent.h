#ifndef DEF_UBEVENT
#define DEF_UBEVENT

#include "TH1D.h"
#include "TH2D.h"

class uBEvent
{
 public:
  
  uBEvent();
  uBEvent(Int_t run, Int_t subrun, Int_t evtnum, TH2D* hreadOut);
  ~uBEvent();

  //
  // Setters
  //

  void SetFullReadOut(TH2D *h);
  void SetRunNumber(Int_t run);
  void SetSubRunNumber(Int_t subrun);
  void SetEvtNumber(Int_t evtnum);

  //
  // Getters
  //

  TH1D* GetWaveForm(Int_t chnum) const;
  TH1D* GetCorrectedWF(Int_t chnum);
  TH2D* GetFullReadOut() const;
  Int_t GetEvtNumber() const;
  Int_t GetRunNumber() const;
  Int_t GetSubRunNumber() const;
  Int_t GetPMTNumer(Int_t chnum);   // Get the correspondancy #PMT / #channel
  Int_t GetChlNumber(Int_t pmtnum); //
  TH1D* GetWFderivative(Int_t chnum) const;// Get the deriavtive of a given WF
  Double_t GetBaseline(Int_t chnum) const;// retrieve the baseline for a given channel w/o re-eavluating it
  
  //
  // Others 
  //

  void DrawFullReadOut();              // Draws the 2D histogram of the radout and the superimposed WF
  double EvalBaseline(Int_t chnum);    // Eveluates the baseline for a given chanel
  double ShowEvalBaseline(Int_t chnum);// Shows how the baseline is evauated
  void GetBaselines();                 // Fills the baselines vector for the event
  void LocatePulses();                 // Locate pulses present on multiple WF at the same time;
  TH1D* SumWF();

 private:
  Int_t Run;
  Int_t SubRun;
  Int_t EvtNumber;
  Int_t N_PMT;
  Int_t N_samples;
  TH2D* hReadOut;
  std::vector<Double_t> baselines;
};

#endif
