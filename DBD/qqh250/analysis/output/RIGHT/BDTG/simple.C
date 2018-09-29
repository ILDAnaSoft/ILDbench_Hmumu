/*
stackwithcut.C

write histogram and evaluate significance
this macro does not use hev, instead of calculated weight
*/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtupleD.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

using namespace std;

struct procinfo{
  int procid;
  double crosssection;
  int epol;
  int ppol;
  int prockind;
};

procinfo g_procinfo[300000];    // map is not usable in ROOT??

void MakeProcinfo( const char *file ){

  memset(g_procinfo,0, sizeof(g_procinfo));
  
  ifstream in(file);
  
  //TString str;                                                          
  TString str_proc;
  TString str_pol;
  procinfo tmpinfo;
  while(!in.eof()){
    //in >> tmpinfo.procid >> tmpinfo.crosssection >> str >> tmpinfo.prockind;
    //in >> tmpinfo.procid >> str_proc >> str_pol >> tmpinfo.crosssection; //procid500.txt specific
    in >> str_proc >> str_pol >> tmpinfo.procid >> tmpinfo.crosssection; //procid250.txt specific
    
    if(str_pol == "eL.pL"){tmpinfo.epol = -1; tmpinfo.ppol = -1;}
    else if(str_pol == "eL.pW" || str_pol == "eL.pB"){tmpinfo.epol = -1; tmpinfo.ppol =  0;}
    else if(str_pol == "eL.pR"){tmpinfo.epol = -1; tmpinfo.ppol =  1;}
    else if(str_pol == "eW.pL" || str_pol == "eB.pL"){tmpinfo.epol =  0; tmpinfo.ppol = -1;}
    else if(str_pol == "eW.pW" || str_pol == "eW.pB" || str_pol == "eB.pW" || str_pol == "eB.pB"){tmpinfo.epol =  0; tmpinfo.ppol =  0;}
    else if(str_pol == "eW.pR" || str_pol == "eB.pR"){tmpinfo.epol =  0; tmpinfo.ppol =  1;}
    else if(str_pol == "eR.pL"){tmpinfo.epol =  1; tmpinfo.ppol = -1;}
    else if(str_pol == "eR.pW" || str_pol == "eR.pB"){tmpinfo.epol =  1; tmpinfo.ppol =  0;}
    else if(str_pol == "eR.pR"){tmpinfo.epol =  1; tmpinfo.ppol =  1;}

    if( str_proc == "ffh_mumu" ) tmpinfo.prockind = 1;
    else if( str_proc == "e1e1h" || str_proc == "e2e2h" || str_proc == "e3e3h" ||
	     str_proc == "nnh"   || str_proc == "qqh" ) tmpinfo.prockind = 2;
    else if( str_proc.BeginsWith("2f_") == true ) tmpinfo.prockind = 3;
    //else if( tmpinfo.procid >= 37473 && tmpinfo.procid <= 37488 ) tmpinfo.prockind = 4;//aa_2f
    //else if( tmpinfo.procid >= 37489 && tmpinfo.procid <= 37544 ) tmpinfo.prockind = 5;//3f
    else if( str_proc.BeginsWith("4f_") == true ) tmpinfo.prockind = 4;
    else if( str_proc.BeginsWith("aa_") == true ) tmpinfo.prockind = 5; //aa_2f
    else if( str_proc.BeginsWith("ae_") == true || str_proc.BeginsWith("ea_") == true ) tmpinfo.prockind = 6;//1f_3f
    else tmpinfo.prockind = 7;

    g_procinfo[tmpinfo.procid] = tmpinfo;
  }

}

void MakeAllWithWeight(
		       TTree *nt, int nbin, double nlow, double nhigh,
		       const char *varexp, const char *selection,
		       double epol = 0.9, double ppol = 0.35,
		       const char *idname = "processid",
		       bool isiddouble = false, bool tau = false, bool taupol = false )
{
  if( selection == NULL ) selection = "1"; // select all events

  // Make TreeFormulae
  TTreeFormula tfexp( "tfexp", varexp, nt );
  TTreeFormula tfsel( "tfsel", selection, nt );
  double idd;
  int idi;

  //Make histograms
  TH1F *h[3];
  h[ 0] = new TH1F("all","all",nbin,nlow,nhigh);
  h[ 1] = new TH1F("sig","sig",nbin,nlow,nhigh);
  h[ 2] = new TH1F("bkg","bkg",nbin,nlow,nhigh);

  int higgsdecay1;
  nt->SetBranchAddress("higgsdecay1", &higgsdecay1);
  int type;
  nt->SetBranchAddress("type", &type);
  double weight;
  nt->SetBranchAddress("weight", &weight);

  // stats                                                                
  double entryKind[3];
  memset(entryKind,0,sizeof(entryKind));
  int counter[3];
  memset(counter,0,sizeof(counter));

  if(isiddouble)
    nt->SetBranchAddress(idname, &idd);
  else
    nt->SetBranchAddress(idname, &idi);

  for( int i = 0; i < nt->GetEntries(); i++ ){
    nt->GetEntry(i);
    if( i % 100000 == 0 && i > 0 ) cout << i << " entries processed." << endl;

    if( tfsel.EvalInstance(i) == 0 ) continue;

    // obtain event info
    int procid = (isiddouble ? (int)idd : idi);
    procinfo &info = g_procinfo[procid];

    // fill stats
    if( info.prockind == 1 && ( type >= 1 && type <= 5 ) ){//sig
      counter[0]++;
      entryKind[0] += 2*weight;
    }
    else{//bkg
      counter[1]++;
      entryKind[1] += 2*weight;
    }

    h[ 0]->Fill(tfexp.EvalInstance(i),2*weight); //all
    if( info.prockind == 1 && ( type >= 1 && type <= 5 ) ) h[ 1]->Fill(tfexp.EvalInstance(i),2*weight);//sig
    else h[ 2]->Fill(tfexp.EvalInstance(i),2*weight);//bkg

  }

  // write histograms and stats
  TLegend *leg = new TLegend(0.7,0.7,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  h[0]->SetLineColor(kBlack);h[0]->SetMarkerColor(kBlack);h[0]->Draw();//ALL
  leg->AddEntry(h[0], "All", "lp");
  if( counter[0] != 0 ){//sig
    h[1]->SetLineColor(kBlue);h[1]->SetMarkerColor(kBlue);h[1]->SetLineWidth(3);
    h[1]->Draw("same");
    leg->AddEntry(h[1], "signal", "lp");
  }
  if( counter[1] != 0 ){//bkg
    h[2]->SetLineColor(kRed);h[2]->SetMarkerColor(kRed);
    h[2]->Draw("same");
    leg->AddEntry(h[2], "background", "lp");
  }
  leg->Draw();

  /*
  //normalize to 1
  double norm_sig = h[1]->Integral();
  double norm_bkg = h[2]->Integral();
  h[1]->Scale(1./norm_sig); h[2]->Scale(1./norm_bkg);
  h[1]->Draw(); h[2]->Draw("same");
  */

  //show statistics
  cout << setw(12) << "signal" << setw(12) << "background" << endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << endl;

  cout << "signal = " << entryKind[0] << ", background = " << entryKind[1] << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + entryKind[1] );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;

}
