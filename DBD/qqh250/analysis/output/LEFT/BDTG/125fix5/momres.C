void momres(){
  Double_t xmax = 1e-3, xmin = 1e-6, ymax = 100, ymin = 0;
  TCanvas *c1 = new TCanvas("c1","momres",0,0,600,400);
  TH1F *frame = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  c1->SetLogx();
  frame->SetTitle("qqh250-L; Momentum Resolution (GeV^{-1}); #Delta(#sigma#timesBR)/(#sigma#timesBR) (%)");

  //benchmark numbers
  Double_t x[13]  = {1e-3 , 5e-4 , 3e-4 , 2e-4 , 1e-4 , 5e-5 , 3e-5 , 2e-5 , 1e-5  , 5e-6  , 3e-6  , 2e-6  , 1e-6};
  Double_t y[13]  = {51.38, 46.55, 42.43, 40.32, 36.04, 33.44, 31.38, 29.84, 27.970, 26.594, 25.573, 24.808, 23.952};
  Double_t ex[13] = {0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0     , 0     , 0     , 0     , 0};
  Double_t ey[13] = {0.20 , 0.18 , 0.16 , 0.15 , 0.13 , 0.12 , 0.11 , 0.10 , 0.095 , 0.090 , 0.086 , 0.083 , 0.080};
  
  TGraphErrors *gr1 = new TGraphErrors(13,x,y,ex,ey);
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerStyle(22);
  gr1->Draw("P");

  //result of qqh250-L
  Double_t x2[1] = {2e-5}; 
  Double_t y2[1] = {32.78};
  Double_t ex2[1] = {0};
  Double_t ey2[1] = {0.11};

  TGraphErrors *gr2 = new TGraphErrors(1,x2,y2,ex2,ey2);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerSize(1);
  gr2->SetMarkerStyle(20);
  gr2->Draw("P");

  //legend
  TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
  leg->AddEntry(gr2, "full", "p");
  leg->AddEntry(gr1, "benchmark", "p");
  leg->Draw();
}
