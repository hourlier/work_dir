#include "UBEvent.cxx"


TH2D* GenerateReadOut(){
  TH2D *hreadout = new TH2D("hreadout","hreadout",1500,0,1500,32,0,32);
  TF1 *f = new TF1("f","gaus(0)",0,1500);
  for(int i = 0;i<32;i++){
    f->SetParameters(20,10*i+100,5);
    for(int j = 0;j<1500;j++){
      hreadout->SetBinContent(j+1,i+1,f->Eval(j)+2050);
    }
  }
  return hreadout;
}

void test(){
  cout << "Hello World"<< endl;
  TH2D *hreadout = GenerateReadOut();
  hreadout->Draw("colz");
  uBEvent evt(1,1,1,hreadout);
  TCanvas *c2 = new TCanvas();
  evt.GetFullReadOut()->Draw("colz");
  TCanvas *c3 = new TCanvas();
  evt.GetWaveForm(1)->Draw();
  evt.GetWaveForm(10)->Draw("same");
}
