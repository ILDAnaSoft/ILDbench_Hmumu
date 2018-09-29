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
    in >> tmpinfo.procid >> str_proc >> str_pol >> tmpinfo.crosssection; //procid500.txt specific
    //in >> str_proc >> str_pol >> tmpinfo.procid >> tmpinfo.crosssection;
    
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
    else if( str_proc.BeginsWith("aa_") == true &&
	     !( tmpinfo.procid >= 37473 && tmpinfo.procid <= 37488 ) ) tmpinfo.prockind = 5;//aa_4f
    else if( ( str_proc.BeginsWith("ae_") == true || str_proc.BeginsWith("ea_") == true ) &&
	     !( tmpinfo.procid >= 37489 && tmpinfo.procid <= 37544 ) ) tmpinfo.prockind = 6;//5f
    else tmpinfo.prockind = 7;

    g_procinfo[tmpinfo.procid] = tmpinfo;
  }

}

void MakeAllWithWeight(
		       ofstream& ofs, ofstream& ofs2, TTree *nt,
		       int nbin, double nlow, double nhigh,
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
  TH1F *h[9];
  h[ 0] = new TH1F("all","all",nbin,nlow,nhigh);
  h[ 1] = new TH1F("nnh_mumu","nnh_mumu",nbin,nlow,nhigh);
  h[ 2] = new TH1F("other_h_mumu","other_h_mumu",nbin,nlow,nhigh);
  h[ 3] = new TH1F("other_h_decay","other_h_decay",nbin,nlow,nhigh);
  h[ 4] = new TH1F("2f","2f",nbin,nlow,nhigh);
  //h[ 5] = new TH1F("aa_2f","aa_2f",nbin,nlow,nhigh);
  //h[ 6] = new TH1F("ea_3f","ea_3f",nbin,nlow,nhigh);
  h[ 5] = new TH1F("4f","4f",nbin,nlow,nhigh);
  h[ 6] = new TH1F("aa_4f","aa_4f",nbin,nlow,nhigh);
  h[ 7] = new TH1F("ea_5f","ea_5f",nbin,nlow,nhigh);
  h[ 8] = new TH1F("other","other",nbin,nlow,nhigh);

  int higgsdecay1;
  nt->SetBranchAddress("higgsdecay1", &higgsdecay1);
  int type;
  nt->SetBranchAddress("type", &type);
  double weight;
  nt->SetBranchAddress("weight", &weight);

  // stats                                                                
  double entryKind[15];
  memset(entryKind,0,sizeof(entryKind));
  int counter[15];
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
    if( info.prockind == 1 && ( type == 12 || type == 14 || type == 16 ) ){//nnh, h->mumu
      counter[0]++;
      entryKind[0] += 2.*weight;
    }
    else if( info.prockind == 1 && !( type == 12 || type == 14 || type == 16 ) ){//qqh/llh, h->mumu
      counter[1]++;
      entryKind[1] += 2.*weight;
    }
    else if( info.prockind == 2 ){//ffh, h->other
      counter[2]++;
      entryKind[2] += 2.*weight;
    }
    else if( info.prockind == 3 ){//2f
      counter[3]++;
      entryKind[3] += 2.*weight;
    }
    else if( info.prockind == 4 ){//4f
      counter[4]++;
      entryKind[4] += 2.*weight;
    }
    else if( info.prockind == 5 ){//aa_4f
      counter[5]++;
      entryKind[5] += 2.*weight;
    }
    else if( info.prockind == 6 ){//ea_5f
      counter[6]++;
      entryKind[6] += 2.*weight;
    }
    else{//other
      counter[7]++;
      entryKind[7] += 2.*weight;
    }

    h[ 0]->Fill(tfexp.EvalInstance(i),2.*weight); //all
    if( info.prockind == 1 && ( type == 12 || type == 14 || type == 16 ) ) h[ 1]->Fill(tfexp.EvalInstance(i),2.*weight);//nnh, h->mumu
    else if( info.prockind == 1 && !( type == 12 || type == 14 || type == 16 ) ) h[ 2]->Fill(tfexp.EvalInstance(i),2.*weight);//qqh/llh, h->mumu
    else if( info.prockind == 2 ) h[ 3]->Fill(tfexp.EvalInstance(i),2.*weight);//ffh, h->other
    else if( info.prockind == 3 ) h[ 4]->Fill(tfexp.EvalInstance(i),2.*weight);//2f
    //else if( info.prockind == 4 ) h[ 5]->Fill(tfexp.EvalInstance(i),2.*weight);//aa_2f
    //else if( info.prockind == 5 ) h[ 6]->Fill(tfexp.EvalInstance(i),2.*weight);//ea_3f
    else if( info.prockind == 4 ) h[ 5]->Fill(tfexp.EvalInstance(i),2.*weight);//4f
    else if( info.prockind == 5 ) h[ 6]->Fill(tfexp.EvalInstance(i),2.*weight);//aa_4f
    else if( info.prockind == 6 ) h[ 7]->Fill(tfexp.EvalInstance(i),2.*weight);//ea_5f
    else h[ 8]->Fill(tfexp.EvalInstance(i),2.*weight); //other

  }

  // write histograms and stats
  TLegend *leg = new TLegend(0.7,0.5,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  h[0]->SetLineColor(kBlack);h[0]->SetMarkerColor(kBlack);h[0]->Draw();//ALL
  leg->AddEntry(h[0], "All", "lp");
  if( counter[0] != 0 ){//nnh, h->mumu
    h[1]->SetLineColor(kBlue);h[1]->SetMarkerColor(kBlue);h[1]->SetLineWidth(3);
    h[1]->Draw("same");
    leg->AddEntry(h[1], "#nu#nuh, h->#mu#mu", "lp");
  }
  if( counter[1] != 0 ){//qqh/llh, h->mumu
    h[2]->SetLineColor(kBlue);h[2]->SetMarkerColor(kBlue);h[2]->SetLineStyle(2);
    h[2]->Draw("same");
    leg->AddEntry(h[2], "qqh/llh, h->#mu#mu", "lp");
  }
  if( counter[2] != 0 ){//ffh, h->others
    h[3]->SetLineColor(kGray);h[3]->SetMarkerColor(kGray);
    h[3]->Draw("same");
    leg->AddEntry(h[3], "ffh, h->other", "lp");
  }
  if( counter[3] != 0 ){//2f
    h[4]->SetLineColor(kRed);h[4]->SetMarkerColor(kRed);
    h[4]->Draw("same");
    leg->AddEntry(h[4], "2f", "lp");
  }
  if( counter[4] != 0 ){//4f
    h[5]->SetLineColor(kTeal+4);h[5]->SetMarkerColor(kTeal+4);
    h[5]->Draw("same");
    leg->AddEntry(h[5], "4f", "lp");
  }
  if( counter[5] != 0 ){//aa_4f
    h[6]->SetLineColor(kTeal+4);h[6]->SetMarkerColor(kTeal+4);h[6]->SetLineStyle(2);
    h[6]->Draw("same");
    leg->AddEntry(h[6], "#gamma#gamma->4f", "lp");
  }
  if( counter[6] != 0 ){//ea_5f
    h[7]->SetLineColor(kViolet);h[7]->SetMarkerColor(kViolet);
    h[7]->Draw("same");
    leg->AddEntry(h[7], "e#gamma->5f", "lp");
  }
  /*
  if( counter[7] != 0 ){//other
    h[8]->SetLineColor(kYellow);h[8]->SetMarkerColor(kYellow);
    h[8]->Draw("same");
    leg->AddEntry(h[8], "strange", "lp");
  }
  */
  leg->Draw();

  //show statistics
  cout << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
       << setw(12) << "2f" << setw(12) << "4f" << setw(12) << "aa_4f" << setw(12) << "5f" 
       << setw(12) << "strange" << endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
       << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5]
       << setw(12) << entryKind[6] << setw(12) << entryKind[7] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
       << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5]
       << setw(12) << counter[6] << setw(12) << counter[7] << endl;
  ofs << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
      << setw(12) << "2f" << setw(12) << "4f" << setw(12) << "aa_4f" << setw(12) << "5f" 
      << setw(12) << "strange" << endl;
  ofs << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
      << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5]
      << setw(12) << entryKind[6] << setw(12) << entryKind[7] << endl;
  double bkg_tot = 0.;
  for( int i = 1; i <= 7; i++ ){
    bkg_tot += entryKind[i];
  }
  cout << "signal = " << entryKind[0] << ", background = " << bkg_tot << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + bkg_tot );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;
  ofs2 << signi << " " << 100./signi << "%" << endl;

}
