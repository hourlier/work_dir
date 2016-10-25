#include "UBEvent.h"
#include "TH1D.h"
#include "TH2D.h"

//
// constructor/destructor
//

uBEvent::uBEvent(){
  N_PMT = 32;
  N_samples = 1500;
  hReadOut = new TH2D("hReadOut","hReadOut;ticks;PMT",N_samples,0,N_samples,N_PMT,0,N_PMT);
  Run = -1;
  SubRun = -1;
  EvtNumber = -1;
}

uBEvent::uBEvent(Int_t run, Int_t subrun, Int_t evtnum, TH2D* hreadOut){
  N_PMT = 32;
  N_samples = 1500;
  Run = run;
  SubRun = subrun;
  EvtNumber = evtnum;
  hReadOut = (TH2D*)hreadOut->Clone("hReadOut");
  hReadOut->SetNameTitle("hReadOut",Form("hReadOut_%d_%d_%d;ticks;PMT",Run, SubRun, EvtNumber));
  GetBaselines();
}

uBEvent::~uBEvent(){hReadOut->Delete();}

//
// Setters
//

void uBEvent::SetFullReadOut(TH2D *h){hReadOut = (TH2D*)h->Clone("hReadOut");}

void uBEvent::SetRunNumber(Int_t run){Run = run;}

void uBEvent::SetSubRunNumber(Int_t subrun){SubRun = subrun;}

void uBEvent::SetEvtNumber(Int_t evtnum){EvtNumber = evtnum;}

//
// Getters
//


TH1D* uBEvent::GetWaveForm(Int_t chnum) const {
  TH1D *hwaveform = hReadOut->ProjectionX(Form("hWF_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),chnum+1,chnum+1);
  //hwaveform->SetMaximum(4095);
  //hwaveform->SetMinimum(2000);
  return hwaveform;
}

TH2D* uBEvent::GetFullReadOut() const {TH2D *hFull = (TH2D*)hReadOut->Clone("hFull");return hFull;}

Int_t uBEvent::GetEvtNumber() const { Int_t num = EvtNumber; return num;}

Int_t uBEvent::GetRunNumber() const {Int_t run = Run; return run;}

Int_t uBEvent::GetSubRunNumber() const {Int_t subrun = SubRun; return subrun;}

Int_t uBEvent::GetPMTNumer(Int_t chnum){return chnum+1;}

Int_t uBEvent::GetChlNumber(Int_t pmtnum){return pmtnum-1;}

TH1D* uBEvent::GetWFderivative(Int_t chnum) const {
  TH1D *hWF = GetWaveForm(chnum);
  TH1D *hDerivative = new TH1D("hDerivative",Form("hDerivative_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),hWF->GetNbinsX(),hWF->GetXaxis()->GetXmin(),hWF->GetXaxis()->GetXmax());
  hDerivative->SetNameTitle("hDerivative",Form("hDerivative_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum));
  hDerivative->SetBinContent(1,0);
  for(int i = 1;i< hWF->GetNbinsX();i++){
    hDerivative->SetBinContent(i+1,0.5*(hWF->GetBinContent(i+1)-hWF->GetBinContent(i)));
  }
  return hDerivative;
}

Double_t uBEvent::GetBaseline(Int_t chnum) const {double BL = baselines.at(chnum); return BL;}

//
// Others
//

void uBEvent::DrawFullReadOut(){
  TCanvas *cReadOut = new TCanvas("cReadOut","cReadOut",1000,500);
  cReadOut->Divide(2,1);
  cReadOut->cd(1);
  hReadOut->Draw("colz");
  TBox *beam = new TBox(210,0,340,32);
  beam->SetLineColor(2);
  beam->SetLineWidth(1);
  beam->SetFillStyle(0);
  beam->Draw("same");
  cReadOut->cd(2);
  TH1D *hWF = GetCorrectedWF(0);
  hWF->Draw();
  hWF->GetYaxis()->SetRangeUser(-200,2048);
  for(int i = 1;i<N_PMT;i++){
    GetCorrectedWF(i)->Draw("same");
  }
}

double uBEvent::EvalBaseline(Int_t chnum){
  //TCanvas *cBaseline = new TCanvas("cBaseline","cBaseline",1000,1000);
  //cBaseline->Divide(2,2);
  //cBaseline->cd(1);
  TH1D *hWF = GetWaveForm(chnum);
  //hWF->Draw();

  //cBaseline->cd(2);
  TH1D *hDerivative = GetWFderivative(chnum);
  //hDerivative->Draw();

  int LasHighTick = 0;
  int DT = 50;
  double threshold = 1.25;
  int TimeSinceHigh = DT-5;
  std::vector<TBox*> box;
  bool boxopen = false;
  

  TH1D *hBaseline = new TH1D("hBaseline",Form("hBaseline_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),100,2000,2100);
  TH1D *hFirstBL = new TH1D("hFirstBL",Form("hFirstBL_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),100,2000,2100);

  for(int i = 0;i<hDerivative->GetNbinsX();i++){
    if(TMath::Abs(hDerivative->GetBinContent(i+1)) < threshold){
      if(TimeSinceHigh >= DT && boxopen == false){
	TBox *b = new TBox();
	b->SetY1(-2);
	b->SetY2(2);
	b->SetLineColor(2);
	b->SetLineWidth(2);
	b->SetFillStyle(0);
	b->SetX1(i);
	box.push_back(b);
	boxopen = true;
	TimeSinceHigh = i-LasHighTick;
	//cout << "order to open box " << i << endl;
      }
      TimeSinceHigh =i-LasHighTick;
    }
    else{
      if(boxopen == true){
	boxopen = false;
	box.at(box.size()-1)->SetX2(i-DT);
	//box.at(box.size()-1)->Draw();
	//cout << "order to close box " << i << endl;
      }
      boxopen = false;
      TimeSinceHigh = 0;
      LasHighTick = i;
    }
    if(i == hDerivative->GetNbinsX()-1 && boxopen == true){
      box.at(box.size()-1)->SetX2(i);
      //box.at(box.size()-1)->Draw();
      //cout << "order to close box " << i << endl;
      boxopen = false;
    }
    if(boxopen){
      hBaseline->Fill(hWF->GetBinContent(i+1));
      if(box.size() == 1){hFirstBL->Fill(hWF->GetBinContent(i+1));}
    }
  }
  
  //cBaseline->cd(3);
  //cBaseline->cd(3)->SetLogy();
  //hBaseline->Draw();

  //TF1 *fBaseline = new TF1("fBaseline","gaus(0)",2000,2100);
  //fBaseline->SetParameters(hBaseline->GetBinContent(hBaseline->GetMaximumBin()),hBaseline->GetMean(),hBaseline->GetRMS());
  //hBaseline->Fit("fBaseline","rlon");
  hFirstBL->SetLineColor(4);
  hFirstBL->SetLineWidth(2);
  //hFirstBL->Draw("same");

  double baseline;
  if(box.size()>1){
    hBaseline->Add(hFirstBL,-1);
  }

  if(TMath::Abs((hFirstBL->GetMean()-hBaseline->GetMean()))/sqrt(pow(hFirstBL->GetRMS(),2)+pow(hBaseline->GetRMS(),2)) > 1.5){
    baseline = hFirstBL->GetMean();
  }
  else{
    hBaseline->Add(hFirstBL,1);
    baseline = hBaseline->GetMean();
  }

  hBaseline->Delete();
  hFirstBL->Delete();
  hDerivative->Delete();

  return baseline;
}

double uBEvent::ShowEvalBaseline(Int_t chnum){
  TCanvas *cBaseline = new TCanvas("cBaseline","cBaseline",1000,1000);
  cBaseline->Divide(2,2);
  cBaseline->cd(1);
  TH1D *hWF = GetWaveForm(chnum);
  hWF->Draw();
  
  cBaseline->cd(2);
  TH1D *hDerivative = GetWFderivative(chnum);
  hDerivative->Draw();
  int LasHighTick = 0;
  int DT = 50;
  double threshold = 1.25;
  int TimeSinceHigh = DT-5;
  std::vector<TBox*> box;
  bool boxopen = false;


  TH1D *hBaseline = new TH1D("hBaseline",Form("hBaseline_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),100,2000,2100);
  TH1D *hFirstBL = new TH1D("hFirstBL",Form("hFirstBL_%d_%d_%d_%d",Run,SubRun,EvtNumber,chnum),100,2000,2100);

  for(int i = 0;i<hDerivative->GetNbinsX();i++){
    if(TMath::Abs(hDerivative->GetBinContent(i+1)) < threshold){
      if(TimeSinceHigh >= DT && boxopen == false){
        TBox *b = new TBox();
        b->SetY1(-2);
        b->SetY2(2);
        b->SetLineColor(2);
        b->SetLineWidth(2);
        b->SetFillStyle(0);
        b->SetX1(i);
        box.push_back(b);
        boxopen = true;
        TimeSinceHigh = i-LasHighTick;
        //cout << "order to open box " << i << endl;
      }
      TimeSinceHigh =i-LasHighTick;
    }
    else{
      if(boxopen == true){
        boxopen = false;
        box.at(box.size()-1)->SetX2(i-DT);
        box.at(box.size()-1)->Draw();
        //cout << "order to close box " << i << endl;
      }
      boxopen = false;
      TimeSinceHigh = 0;
      LasHighTick = i;
    }
    if(i == hDerivative->GetNbinsX()-1 && boxopen == true){
      box.at(box.size()-1)->SetX2(i);
      box.at(box.size()-1)->Draw();
      //cout << "order to close box " << i << endl;
      boxopen = false;
    }
    if(boxopen){
      hBaseline->Fill(hWF->GetBinContent(i+1));
      if(box.size() == 1){hFirstBL->Fill(hWF->GetBinContent(i+1));}
    }
  }

  cBaseline->cd(3);
  cBaseline->cd(3)->SetLogy();
  hBaseline->Draw();

  hFirstBL->SetLineColor(4);
  hFirstBL->SetLineWidth(2);
  hFirstBL->Draw("same");

  double baseline;
  if(box.size()>1){
    hBaseline->Add(hFirstBL,-1);
  }

  if(TMath::Abs((hFirstBL->GetMean()-hBaseline->GetMean()))/sqrt(pow(hFirstBL->GetRMS(),2)+pow(hBaseline->GetRMS(),2)) > 1.5){
    baseline = hFirstBL->GetMean();
    cout << "potentially bad baseline" << endl;
  }
  else{
    hBaseline->Add(hFirstBL,1);
    baseline = hBaseline->GetMean();
  }


  return baseline;
}

void uBEvent::GetBaselines(){
  baselines.clear();
  double BL;
  for(int i = 0;i<N_PMT;i++){
    BL = EvalBaseline(i);
    baselines.push_back(BL);
  }
}



TH1D* uBEvent::GetCorrectedWF(Int_t chnum){
  TH1D *hWF = GetWaveForm(chnum);
  TH1D *hBL = (TH1D*)hWF->Clone("hBL");
  double baseline = GetBaseline(chnum);
  for(int i = 0;i<hBL->GetNbinsX();i++){
    hBL->SetBinContent(i+1,baseline);
  }
  hWF->Add(hBL,-1);
  return hWF;
}

TH1D* uBEvent::SumWF(){
  TH1D *hWF = GetCorrectedWF(0);
  for(int i = 1;i< N_PMT;i++){
    hWF->Add(GetCorrectedWF(i),1);
  }
  return hWF;
}
