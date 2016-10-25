double SPEshape(double *x, double *par){
  double baseline = par[0];
  double amplitude = par[1];
  double t0 = par[2];
  double sigma_1 = par[3];
  double sigma_2 = par[4];

  if(x[0] < t0){
    return baseline+amplitude*TMath::Gaus(x[0],t0,sigma_1);
  }
  else{
    return baseline+amplitude*TMath::Gaus(x[0],t0,sigma_2);
  }
}

double SPEshape(double x, double amplitude, double T0, double sigma_rise, double sigma_fall){
  if(x<T0){
    return amplitude*TMath::Gaus(x,T0,sigma_rise);
  }
  else{
    return amplitude*TMath::Gaus(x,T0,sigma_fall);
  }
}

double WFshape(double *x, double *par){
  int NPE = par[0];
  double baseline = par[1];
  double value = baseline;
  for(int i = 0;i<NPE;i++){
    value += SPEshape(x[0],par[2+4*i], par[3+4*i],par[4+4*i], par[5+4*i]);
  }
  return value;
}

TF1* GetfWF(int NPE){
  TF1 *fWF = new TF1("fWF",WFshape, 0,1500,4*NPE+2);
  fWF->FixParameter(0,NPE);
  return fWF;
}

TH1D* GetDerivative(TH1D *h){
  TH1D *hDerivative = new TH1D("hDerivative","hDerivative",h->GetNbinsX(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  hDerivative->SetBinContent(1,0);
  for(int i = 1;i<h->GetNbinsX();i++){
    hDerivative->SetBinContent(i+1,0.5*(h->GetBinContent(i+1)-h->GetBinContent(i)));
  }
  return hDerivative;
}

TH1D* EvaluateBaseline(TH1D *h){
  TH1D *hBaseline = new TH1D("hBaseline","hBaseline",4096,0,4096);
  TH1D *hDerivative = GetDerivative(h);
  
  for(int i = 0;i<h->GetNbinsX();i++){
    if(TMath::Abs(hDerivative->GetBinContent(i+1)) < 1){
      hBaseline->Fill(h->GetBinContent(i+1));
    }
  }
  return hBaseline;
}
