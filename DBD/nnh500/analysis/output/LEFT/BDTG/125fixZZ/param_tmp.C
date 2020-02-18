#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TStyle.h"
using namespace RooFit;

void param()
{
  gStyle->SetOptFit(1111);
  TFile *file = TFile::Open("alldata_after.root");
  TTree *tree = (TTree*) file->Get("datatest");

  RooRealVar processid("processid","processid",0,300000);
  RooRealVar type("type","type",0,100);
  RooRealVar mumu_mass("mumu_mass","mumu_mass",120,130);
  RooRealVar BDTGoutput("BDTGoutput","BDTGoutput",-1,1);
  RooRealVar weight2("weight2","weight2",0,200);
  RooDataSet *data = new RooDataSet( "data", "template data", RooArgSet(processid,type,mumu_mass,BDTGoutput,weight2), Import(*tree) );
  RooDataSet *wdata = new RooDataSet( data->GetName(), data->GetTitle(), data, *data->get(), 0, weight2.GetName());

  //prepare canvas
  TCanvas* c = new TCanvas("practice","pratice",1200,800) ;
  c->Divide(2,2);
  c->cd(1);

  //choose only signal
  RooDataSet *sdata = wdata->reduce( Cut("((processid>=108161&&processid<=108164)&&(type==12||type==14||type==16))&&mumu_mass>120&&BDTGoutput>BDTGBDTG") );
  RooDataSet *sdata_reduced = sdata->reduce(SelectVars(mumu_mass));
  //RooDataHist *sdata_binned = sdata_reduced->binnedClone();
  RooPlot* xframe = mumu_mass.frame();
  sdata_reduced->plotOn(xframe);
  //sdata_binned->plotOn(xframe);
  xframe->Draw();

  c->cd(2);

  //use Crystall Ball function for signal modelling
  RooRealVar mean("mean","mean",125);
  RooRealVar sigma("sigma","sigma",0.1,0,5);
  RooRealVar alpha("alpha","alpha",1,-10,10);
  RooRealVar n("n","n",1,-10,10);
  RooCBShape CB("CB","CB function",mumu_mass,mean,sigma,alpha,n);
  RooRealVar gmean("gmean","gmean",125,120,130);
  RooRealVar gsigma("gsigma","gsigma",0.1,0,5);
  RooGaussian gauss("gauss","gauss",mumu_mass,gmean,gsigma);
  RooRealVar CBfrac("CBfrac","CBfrac",0.5,0,1);
  RooAddPdf CBgauss("CBgauss","CBgauss",RooArgSet(CB,gauss),RooArgSet(CBfrac));

  RooPlot* xframe2 = mumu_mass.frame(Title("Data and CB+gauss"));
  sdata_reduced->plotOn(xframe2,DataError(RooAbsData::SumW2));
  //sdata_binned->plotOn(xframe2);
  RooFitResult *fitresult;
  fitresult = CBgauss.fitTo(*sdata_reduced,Range(120,130),Minos(kTRUE),SumW2Error(kTRUE),Save());
  //fitresult = CBgauss.fitTo(*sdata_binned,Range(120,130),Minos(kTRUE),Save());
  CBgauss.plotOn(xframe2,LineColor(kRed));
  xframe2->Draw();
  cout << "signal fitting result with BDTGcut BDTGBDTG" << endl;

  c->cd(3);

  //choose only background
  mumu_mass.setBins(NBINSBKGNBINSBKG);
  RooDataSet *bdata = wdata->reduce( Cut("!((processid>=108161&&processid<=108164)&&(type==12||type==14||type==16))&&mumu_mass>120&&BDTGoutput>BDTGBDTG") );
  RooDataSet *bdata_reduced = bdata->reduce(SelectVars(mumu_mass));
  RooDataHist *bdata_binned = bdata_reduced->binnedClone();
  RooPlot* xframe3 = mumu_mass.frame(Bins(NBINSBKGNBINSBKG));
  bdata_binned->plotOn(xframe3,Bins(NBINSBKGNBINSBKG));
  xframe3->Draw();

  c->cd(4);

  //use pol1 for background modelling
  RooRealVar p0("p0","p0",0,-100,100);
  RooChebychev poly("poly","poly",mumu_mass,RooArgSet(p0));

  RooPlot* xframe4 = mumu_mass.frame(Title("Data and Linear Function"),Bins(NBINSBKGNBINSBKG)) ;
  bdata_binned->plotOn(xframe4,Bins(NBINSBKGNBINSBKG),DataError(RooAbsData::SumW2));
  RooFitResult *fitresult_b;
  fitresult_b = poly.fitTo(*bdata_binned,Range(120,130),Minos(kTRUE),SumW2Error(kTRUE),Save());
  poly.plotOn(xframe4);
  xframe4->Draw();
  cout << "background fitting result with BDTGcut BDTGBDTG" << endl;
  p0.Print();

  cout << "==========final numbers==========" << endl;
  cout << "# sig events: " << sdata->sumEntries() << endl;
  cout << "# bkg events: " << bdata->sumEntries() << endl;
  cout << "CB mean " << mean.Print() << endl;
  cout << "CB sigma " << sigma.Print() << endl;
  cout << "CB alpha " << alpha.Print() << endl;
  cout << "CB n " << n.Print() << endl;
  cout << "CB frac " << CBfrac.Print() << endl;
  cout << "gauss mean " << gmean.Print() << endl;
  cout << "gauss sigma " << gsigma.Print() << endl;
  cout << "CBgauss chi2/ndf " << xframe2->chiSquare() << endl;
  cout << "p0 " << p0.Print() << endl;
  cout << "pol1 chi2/ndf " << xframe4->chiSquare() << endl;

  ofstream output("output_BDTGBDTG.txt");
  TString No;
  if(BDTGBDTG==-0.95) No = "01";
  else if(BDTGBDTG==-0.9 ) No = "02";
  else if(BDTGBDTG==-0.85) No = "03";
  else if(BDTGBDTG==-0.8 ) No = "04";
  else if(BDTGBDTG==-0.75) No = "05";
  else if(BDTGBDTG==-0.7 ) No = "06";
  else if(BDTGBDTG==-0.65) No = "07";
  else if(BDTGBDTG==-0.6 ) No = "08";
  else if(BDTGBDTG==-0.55) No = "09";
  else if(BDTGBDTG==-0.5 ) No = "10";
  else if(BDTGBDTG==-0.45) No = "11";
  else if(BDTGBDTG==-0.4 ) No = "12";
  else if(BDTGBDTG==-0.35) No = "13";
  else if(BDTGBDTG==-0.3 ) No = "14";
  else if(BDTGBDTG==-0.25) No = "15";
  else if(BDTGBDTG==-0.2 ) No = "16";
  else if(BDTGBDTG==-0.15) No = "17";
  else if(BDTGBDTG==-0.1 ) No = "18";
  else if(BDTGBDTG==-0.05) No = "19";
  else if(BDTGBDTG==0    ) No = "20";
  else if(BDTGBDTG==0.05 ) No = "21";
  else if(BDTGBDTG==0.1  ) No = "22";
  else if(BDTGBDTG==0.15 ) No = "23";
  else if(BDTGBDTG==0.2  ) No = "24";
  else if(BDTGBDTG==0.25 ) No = "25";
  else if(BDTGBDTG==0.3  ) No = "26";
  else if(BDTGBDTG==0.35 ) No = "27";
  else if(BDTGBDTG==0.4  ) No = "28";
  else if(BDTGBDTG==0.45 ) No = "29";
  else if(BDTGBDTG==0.5  ) No = "30";
  else if(BDTGBDTG==0.55 ) No = "31";
  else if(BDTGBDTG==0.6  ) No = "32";
  else if(BDTGBDTG==0.65 ) No = "33";
  else if(BDTGBDTG==0.7  ) No = "34";
  else if(BDTGBDTG==0.75 ) No = "35";
  else if(BDTGBDTG==0.8  ) No = "36";
  else if(BDTGBDTG==0.85 ) No = "37";
  else if(BDTGBDTG==0.9  ) No = "38";
  else if(BDTGBDTG==0.95 ) No = "39";
  output << No << " BDTGBDTG "  << sigma.getVal() << " " << alpha.getVal() << " " << n.getVal() << " "
	 << gmean.getVal() << " " << gsigma.getVal() << " " << CBfrac.getVal() << " " << p0.getVal() << endl;

  cout << "finished with BDTGcut BDTGBDTG" << endl;
  cerr << "finished with BDTGcut BDTGBDTG" << endl;
}
