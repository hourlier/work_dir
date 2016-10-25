#include "UBEvent.cxx"
#include "util.h"

TH1D *hSigma = new TH1D("hSigma","hSigma",100,0,10);
TH1D *hAmplitude = new TH1D("hAmplitude","hAmplitude",1000,0,2000);
TH1D *hTime =new TH1D("hTime","hTime",15000,0,1500);
TH1D *hT0 = new TH1D("hT0","hT0",15000,0,1500);

void LookEvt(){
  cout << "Hello World" << endl;
}

TTree* GetOpticalInfo(){
  TTree *T;
  TFile *fIN = TFile::Open("rawdigits.root","READ");
  TDirectory *dir = fIN->GetDirectory("rawdigitwriter");
  T = (TTree*)dir->Get("OpDetWaveforms");
  return T;
}

uBEvent* Get1Evt(int evtnum = 17801){

  TTree *T = GetOpticalInfo();

  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Int_t           opcrate;
  Int_t           opslot;
  Int_t           opfemch;
  Int_t           frame;
  Int_t           sample;
  Int_t           readoutch;
  Int_t           category;
  Int_t           gaintype;
  Double_t        timestamp;
  Double_t        trig_timestamp;
  Double_t        beam_timestamp;
  std::vector<short>   *adcs = 0;

  // Set branch addresses.
  T->SetBranchAddress("run",&run);
  T->SetBranchAddress("subrun",&subrun);
  T->SetBranchAddress("event",&event);
  T->SetBranchAddress("opcrate",&opcrate);
  T->SetBranchAddress("opslot",&opslot);
  T->SetBranchAddress("opfemch",&opfemch);
  T->SetBranchAddress("frame",&frame);
  T->SetBranchAddress("sample",&sample);
  T->SetBranchAddress("readoutch",&readoutch);
  T->SetBranchAddress("category",&category);
  T->SetBranchAddress("gaintype",&gaintype);
  T->SetBranchAddress("timestamp",&timestamp);
  T->SetBranchAddress("trig_timestamp",&trig_timestamp);
  T->SetBranchAddress("beam_timestamp",&beam_timestamp);
  T->SetBranchAddress("adcs", &adcs);

  TH2D *hWF = new TH2D("hWF","hWF;t;fem ch",1500,0,1500,32,0,32);
  for(int i = 0;i<T->GetEntries();i++){
    T->GetEntry(i);
    if(event !=  evtnum) continue;
    if(adcs->size()<60) continue;

    for(int j = 0;j<adcs->size();j++){
      hWF->SetBinContent(j+1,opfemch,adcs->at(j));
    }
  }
  uBEvent* evt = new uBEvent(run, subrun, evtnum, hWF);
  return evt;
}

std::vector<int> GetListOfEvt(){
  TTree *T = GetOpticalInfo();
  int event;
  int prevEvt;
  T->SetBranchAddress("event",&event);
  std::vector<int> EvtList;
  for(int i = 0;i<T->GetEntries();i++){
    T->GetEntry(i);
    if(i>0 && event == prevEvt)continue;
    else if(i>0 && event != prevEvt){
      //cout << event << endl;
      EvtList.push_back(event);
      prevEvt = event;
    }
    else{
      //cout << event << endl;
      EvtList.push_back(event);
      prevEvt =event;
    }
  }
  return EvtList;
}

uBEvent* Load1stEvent(){
  std::vector<int> EvtList = GetListOfEvt();
  uBEvent *evt = Get1Evt(EvtList.at(0));
  evt->DrawFullReadOut();
  return evt;
}

uBEvent* LoadNextEvent(uBEvent* evt){
  std::vector<int> EvtList = GetListOfEvt();
  int currentRunNum = evt->GetRunNumber();
  int currentSubRunNum = evt->GetSubRunNumber();
  int currentEvtNum = evt->GetEvtNumber();
  
  int i = 0;
  while(currentEvtNum != EvtList.at(i) && i < EvtList.size()){
    i++;
  }
  cout << "current evt :  " << currentEvtNum << " \t / \t" << EvtList.at(i) << endl;
  cout << "next evt : " << EvtList.at(i+1) << endl;
  uBEvent *NewEvt = Get1Evt(EvtList.at(i+1));
  NewEvt->DrawFullReadOut();
  return NewEvt;
}

void FindPulse(uBEvent *evt){
  TH1D *hWF;
  Int_t maxPeak;
  double t0  =1500;
  TF1 *f = new TF1("f","gaus(0)",0,1500);
  for(int i = 0;i<32;i++){
    hWF = evt->GetCorrectedWF(i);
    maxPeak = hWF->GetMaximumBin();
    if(hWF->GetBinContent(maxPeak) > 5){
      cout << endl;
      cout << "WF #" << i << endl;
      cout << endl;
      f->SetParameters(hWF->GetBinContent(maxPeak),maxPeak,2);
      hWF->Fit("f","rl","",maxPeak-5,maxPeak+2);
      hAmplitude->Fill(f->GetParameter(0));
      hTime->Fill(f->GetParameter(1));
      hSigma->Fill(f->GetParameter(2));
      if(f->GetParameter(1) < t0)t0 = f->GetParameter(1);
    }
    hT0->Fill(t0);
  }
  TCanvas *cResults = new TCanvas("cResults","cResults",800,600);
  cResults->Divide(2,2);
  cResults->cd(1);
  hAmplitude->Draw();
  cResults->cd(2);
  hTime->Draw();
  cResults->cd(3);
  hSigma->Draw();
  cResults->cd(4);
  hT0->Draw();
}

void FindAllPulses(){
  int Nevt = GetListOfEvt().size();
  uBEvent *evt = Load1stEvent();
  FindPulse(evt);
  for(int i = 1;i<Nevt-1;i++){
    evt = LoadNextEvent(evt);
    FindPulse(evt);
  }
}
