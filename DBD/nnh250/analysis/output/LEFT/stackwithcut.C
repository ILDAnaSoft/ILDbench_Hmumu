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
    else if( str_proc.BeginsWith("4f_") == true ) tmpinfo.prockind = 4;
    else if( str_proc.BeginsWith("aa_") == true ) tmpinfo.prockind = 5;//aa_2f
    else if( str_proc.BeginsWith("ae_") == true || str_proc.BeginsWith("ea_") == true ) tmpinfo.prockind = 6;//1f_3f
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
  TH1F *h[12];
  h[ 0] = new TH1F("all","ILD simulation",nbin,nlow,nhigh);
  h[ 1] = new TH1F("nnh_mumu","nnh_mumu",nbin,nlow,nhigh);
  h[ 2] = new TH1F("other_h_mumu","other_h_mumu",nbin,nlow,nhigh);
  h[ 3] = new TH1F("other_h_decay","other_h_decay",nbin,nlow,nhigh);
  h[ 4] = new TH1F("2f","2f",nbin,nlow,nhigh);
  h[ 5] = new TH1F("4f_nunumumu","4f_nunumumu",nbin,nlow,nhigh);
  h[ 6] = new TH1F("4f_nunutaumu","4f_nunuyaumu",nbin,nlow,nhigh);
  h[ 7] = new TH1F("4f_nunutautau","4f_nunutautau",nbin,nlow,nhigh);
  h[ 8] = new TH1F("4f_other","4f_other",nbin,nlow,nhigh);
  h[ 9] = new TH1F("aa_2f","aa_2f",nbin,nlow,nhigh);
  h[10] = new TH1F("1f_3f","1f_3f",nbin,nlow,nhigh);
  h[11] = new TH1F("other","other",nbin,nlow,nhigh);

  int type;
  nt->SetBranchAddress("type", &type);
  int nn, nm, nt_m;
  nt->SetBranchAddress("nn", &nn);
  nt->SetBranchAddress("nm", &nm);
  nt->SetBranchAddress("nt_m", &nt_m);
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
    if( i % 500000 == 0 && i > 0 ) cout << i << " entries processed." << endl;

    if( tfsel.EvalInstance(i) == 0 ) continue;

    // obtain event info
    int procid = (isiddouble ? (int)idd : idi);
    procinfo &info = g_procinfo[procid];

    // fill stats
    if( info.prockind == 1
	&& ( type == 12 || type == 14 || type == 16 ) ){//nnh, h->mumu
      counter[0]++;
      entryKind[0] += weight;
    }
    else if( info.prockind == 1
	     && !( type == 12 || type == 14 || type == 16 ) ){//qqh/llh, h->mumu
      counter[1]++;
      entryKind[1] += weight;
    }
    else if( info.prockind == 2 ){//ffh, h->other
      counter[2]++;
      entryKind[2] += weight;
    }
    else if( info.prockind == 3 ){//2f
      counter[3]++;
      entryKind[3] += weight;
    }
    else if( info.prockind == 4
	     && ( nn == 2 && nm == 2 ) ){//4f, nunumumu
      counter[4]++;
      entryKind[4] += weight;
    }
    else if( info.prockind == 4
	     && ( nn == 2 && nm == 1 && nt_m == 1 ) ){//4f, nunutaumu
      counter[5]++;
      entryKind[5] += weight;
    }
    else if( info.prockind == 4
	     && ( nn == 2 && nt_m == 2 ) ){//4f, nunutautau
      counter[6]++;
      entryKind[6] += weight;
    }
    else if( info.prockind == 4
	     && !( nn == 2 && nm == 2 )
	     && !( nn == 2 && nm == 1 && nt_m == 1 )
	     && !( nn == 2 && nt_m == 2 ) ){//4f, other
      counter[7]++;
      entryKind[7] += weight;
    }
    else if( info.prockind == 5 ){//aa_2f
      counter[8]++;
      entryKind[8] += weight;
    }
    else if( info.prockind == 6 ){//1f_3f
      counter[9]++;
      entryKind[9] += weight;
    }
    else entryKind[10] += weight;//other

    h[ 0]->Fill(tfexp.EvalInstance(i),weight); //all
    if( info.prockind == 1
	&& ( type == 12 || type == 14 || type == 16 ) ) h[ 1]->Fill(tfexp.EvalInstance(i),weight);//nnh, h->mumu
    else if( info.prockind == 1
	     && !( type == 12 || type == 14 || type == 16 ) ) h[ 2]->Fill(tfexp.EvalInstance(i),weight);//qqh/llh, h->mumu
    else if( info.prockind == 2 ) h[ 3]->Fill(tfexp.EvalInstance(i),weight);//ffh, h->other
    else if( info.prockind == 3 ) h[ 4]->Fill(tfexp.EvalInstance(i),weight);//2f
    else if( info.prockind == 4
	     && ( nn == 2 && nm == 2 ) ) h[ 5]->Fill(tfexp.EvalInstance(i),weight);//4f, nunumumu
    else if( info.prockind == 4
	     && ( nn == 2 && nm == 1 && nt_m == 1 ) ) h[ 6]->Fill(tfexp.EvalInstance(i),weight);//4f, nunutaumu
    else if( info.prockind == 4
	     && ( nn == 2 && nt_m == 2 ) ) h[ 7]->Fill(tfexp.EvalInstance(i),weight);//4f, nunutautau
    else if( info.prockind == 4
	     && !( nn == 2 && nm == 2 )
	     && !( nn == 2 && nm == 1 && nt_m == 1 )
	     && !( nn == 2 && nt_m == 2 ) ) h[ 8]->Fill(tfexp.EvalInstance(i),weight);//4f, other
    else if( info.prockind == 5 ) h[ 9]->Fill(tfexp.EvalInstance(i),weight);//aa_2f
    else if( info.prockind == 6 ) h[10]->Fill(tfexp.EvalInstance(i),weight);//1f_3f
    else h[11]->Fill(tfexp.EvalInstance(i),weight); //other

  }

  // write histograms and stats
  TLegend *leg = new TLegend(0.7,0.5,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  h[0]->SetLineColor(kBlack);h[0]->SetMarkerColor(kBlack);h[0]->Draw();//ALL
  leg->AddEntry(h[0], "All", "lp");
  if( counter[0] != 0 ){//qqh, h->mumu
    h[1]->SetLineColor(kBlue);h[1]->SetMarkerColor(kBlue);h[1]->SetLineWidth(3);
    h[1]->Draw("same");
    leg->AddEntry(h[1], "#nu#nuh, h->#mu#mu", "lp");
  }
  if( counter[1] != 0 ){//nnh/llh, h->mumu
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
  if( counter[4] != 0 ){//4f, nunumumu
    h[5]->SetLineColor(kTeal+4);h[5]->SetMarkerColor(kTeal+4);
    h[5]->Draw("same");
    leg->AddEntry(h[5], "4f #nu#nu#mu#mu", "lp");
  }
  if( counter[5] != 0 ){//4f, nunutaumu
    h[6]->SetLineColor(kTeal+4);h[6]->SetMarkerColor(kTeal+4);h[6]->SetLineStyle(2);
    h[6]->Draw("same");
    leg->AddEntry(h[6], "4f #nu#nu#tau#mu #tau->#mu", "lp");
  }
  if( counter[6] != 0 ){//4f, nunutautau
    h[7]->SetLineColor(kTeal+4);h[7]->SetMarkerColor(kTeal+4);h[7]->SetLineStyle(3);
    h[7]->Draw("same");
    leg->AddEntry(h[7], "4f #nu#nu#tau#tau #tau->#mu", "lp");
  }
  if( counter[7] != 0 ){//4f, other
    h[8]->SetLineColor(kTeal+4);h[8]->SetMarkerColor(kTeal+4);h[8]->SetLineStyle(4);
    h[8]->Draw("same");
    leg->AddEntry(h[8], "4f other", "lp");
  }
  if( counter[8] != 0 ){//aa_2f
    h[9]->SetLineColor(kRed);h[9]->SetMarkerColor(kRed);h[9]->SetLineStyle(2);
    h[9]->Draw("same");
    leg->AddEntry(h[9], "#gamma#gamma->2f", "lp");
  }
  if( counter[9] != 0 ){//1f_3f
    h[10]->SetLineColor(kViolet);h[10]->SetMarkerColor(kViolet);
    h[10]->Draw("same");
    leg->AddEntry(h[10], "e#gamma->3f", "lp");
  }
  leg->Draw();

  //show statistics
  cout << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
       << setw(12) << "2f" << setw(12) << "4f 2nu2mu" << setw(12) << "4f 2nutaumu" << setw(12) << "4f 2nu2tau"
       << setw(12) << "4f other" << setw(12) << "aa_2f" << setw(12) << "1f_3f" << endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
       << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5]
       << setw(12) << entryKind[6] << setw(12) << entryKind[7] << setw(12) << entryKind[8]
       << setw(12) << entryKind[9] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
       << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5]
       << setw(12) << counter[6] << setw(12) << counter[7] << setw(12) << counter[8]
       << setw(12) << counter[9] << endl;
  ofs << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
      << setw(12) << "2f" << setw(12) << "4f 2nu2mu" << setw(12) << "4f 2nutaumu" << setw(12) << "4f 2nu2tau"
      << setw(12) << "4f other" << setw(12) << "aa_2f" << setw(12) << "1f_3f" << endl;
  ofs << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
      << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5]
      << setw(12) << entryKind[6] << setw(12) << entryKind[7] << setw(12) << entryKind[8]
      << setw(12) << entryKind[9] << endl;
  ofs << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
      << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5]
      << setw(12) << counter[6] << setw(12) << counter[7] << setw(12) << counter[8]
      << setw(12) << counter[9] << endl;

  double bkg_tot = 0.;
  for( int i = 1; i <= 10; i++ ){
    bkg_tot += entryKind[i];
  }
  cout << "signal = " << entryKind[0] << ", background = " << bkg_tot << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + bkg_tot );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;
  ofs2 << signi << " " << 100./signi << "%" << endl;

}
