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

void MakeAllWithWeight( TTree *nt, int nbin, double nlow, double nhigh,
			const char *varexp, const char *selection )
{
  if( selection == NULL ) selection = "1"; // select all events

  // Make TreeFormulae
  TTreeFormula tfexp( "tfexp", varexp, nt );
  TTreeFormula tfsel( "tfsel", selection, nt );

  //Make histograms
  TH1F *h[13];
  h[0] = new TH1F("all","ILD simulation",nbin,nlow,nhigh);
  h[1] = new TH1F("nnh_mumu","nnh_mumu",nbin,nlow,nhigh);
  h[2] = new TH1F("other_h_mumu","other_h_mumu",nbin,nlow,nhigh);
  h[3] = new TH1F("other_h_decay","other_h_decay",nbin,nlow,nhigh);
  h[4] = new TH1F("2f","2f",nbin,nlow,nhigh);
  h[5] = new TH1F("4f_2nu2mu","4f_2nu2mu",nbin,nlow,nhigh);
  h[6] = new TH1F("4f_2nu2tau_mu","4f_2nu2tau_mu",nbin,nlow,nhigh);
  h[7] = new TH1F("4f_2nu1mu1tau_mu","4f_2nu1mu1tau_mu",nbin,nlow,nhigh);
  h[8] = new TH1F("4f_other","4f_other",nbin,nlow,nhigh);
  h[9] = new TH1F("aa_4f_2nu2mu","aa_4f_2nu2mu",nbin,nlow,nhigh);
  h[10] = new TH1F("aa_4f_2nu2tau_mu","aa_4f_2nu2tau_mu",nbin,nlow,nhigh);
  h[11] = new TH1F("aa_4f_2nu1mu1tau_mu","aa_4f_2nu1mu1tau_mu",nbin,nlow,nhigh);
  h[12] = new TH1F("aa_4f_other","aa_4f_other",nbin,nlow,nhigh);

  //int higgsdecay1;
  //nt->SetBranchAddress("higgsdecay1", &higgsdecay1);
  int type;
  nt->SetBranchAddress("type", &type);
  int nn, nm, nt_m;
  nt->SetBranchAddress("nn", &nn);
  nt->SetBranchAddress("nm", &nm);
  nt->SetBranchAddress("nt_m", &nt_m);
  double weight;
  nt->SetBranchAddress("weight", &weight);
  int processid;
  nt->SetBranchAddress("processid", &processid);

  // stats                                                                
  double entryKind[15];
  memset(entryKind,0,sizeof(entryKind));
  int counter[15];
  memset(counter,0,sizeof(counter));

  for( int i = 0; i < nt->GetEntries(); i++ ){
    nt->GetEntry(i);
    if( i % 50000 == 0 && i > 0 ) cout << i << " entries processed." << endl;
    if( tfsel.EvalInstance(i) == 0 ) continue;

    // fill histograms and stats
    h[0]->Fill(tfexp.EvalInstance(i),weight); //all
    if( ( processid >= 108161 && processid <= 108164 ) &&
	( type == 12 || type == 14 || type == 16 ) ){//nnh, h->mumu
      counter[0]++;
      entryKind[0] += weight;
      h[1]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 108161 && processid <= 108164 ) &&
	     !( type == 12 || type == 14 || type == 16 ) ){//qqh/llh, h->mumu
      counter[1]++;
      entryKind[1] += weight;
      h[2]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 106515 && processid <= 106526 ) ){//ffh, h->other
      counter[2]++;
      entryKind[2] += weight;
      h[3]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 250101 && processid <= 250108 ) ){//2f
      counter[3]++;
      entryKind[3] += weight;
      h[4]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 250014 && processid <= 250060 ) 
	     && ( nn == 2 && nm == 2 ) ){//4f, 2mu2nu
      counter[4]++;
      entryKind[4] += weight;
      h[5]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 250014 && processid <= 250060 ) 
	     && ( nn == 2 && nt_m == 2 ) ){//4f, 2tau2nu, tau to mu
      counter[5]++;
      entryKind[5] += weight;
      h[6]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 250014 && processid <= 250060 ) 
	     && ( nn == 2 && nm == 1 && nt_m == 1 ) ){//4f, 1tau1mu2nu, tau to mu
      counter[6]++;
      entryKind[6] += weight;
      h[7]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 250014 && processid <= 250060 ) 
	     && !( nn == 2 && nm == 2 )
	     && !( nn == 2 && nm == 1 && nt_m == 1 )
	     && !( nn == 2 && nt_m == 2 ) ){//4f, all other
      counter[7]++;
      entryKind[7] += weight;
      h[8]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 37385 && processid <= 37464 ) 
	     && ( nn == 2 && nm == 2 ) ){//aa_4f, 2mu2nu
      counter[8]++;
      entryKind[8] += weight;
      h[9]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 37385 && processid <= 37464 ) 
	     && ( nn == 2 && nt_m == 2 ) ){//aa_4f, 2tau2nu, tau to mu
      counter[9]++;
      entryKind[9] += weight;
      h[10]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 37385 && processid <= 37464 ) 
	     && ( nn == 2 && nm == 1 && nt_m == 1 ) ){//aa_4f, 1mu1tau2nu, tau to mu
      counter[10]++;
      entryKind[10] += weight;
      h[11]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 37385 && processid <= 37464 )
	     && !( nn == 2 && nm == 2 ) 
	     && !( nn == 2 && nm == 1 && nt_m == 1 )
	     && !( nn == 2 && nt_m == 2 ) ){//aa_4f, other
      counter[11]++;
      entryKind[11] += weight;
      h[12]->Fill(tfexp.EvalInstance(i),weight);
    }
  }

  // write histograms and stats
  TLegend *leg = new TLegend(0.7,0.5,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  h[0]->SetLineColor(kBlack);h[0]->Draw("HIST");//ALL
  leg->AddEntry(h[0], "All", "l");
  if( counter[0] != 0 ){//nnh, h->mumu
    h[1]->SetLineColor(kBlue);h[1]->SetMarkerColor(kBlue);h[1]->SetLineWidth(3);
    h[1]->Draw("same HIST");
    leg->AddEntry(h[1], "#nu#nuh, h->#mu#mu", "l");
  }
  if( counter[1] != 0 ){//qqh/llh, h->mumu
    h[2]->SetLineColor(kBlue);h[2]->SetMarkerColor(kBlue);h[2]->SetLineStyle(2);
    h[2]->Draw("same HIST");
    leg->AddEntry(h[2], "qqh/llh, h->#mu#mu", "l");
  }
  if( counter[2] != 0 ){//ffh, h->others
    h[3]->SetLineColor(kGray);h[3]->SetMarkerColor(kGray);
    h[3]->Draw("same HIST");
    leg->AddEntry(h[3], "ffh, h->other", "l");
  }
  if( counter[3] != 0 ){//2f
    h[4]->SetLineColor(kRed);h[4]->SetMarkerColor(kRed);
    h[4]->Draw("same HIST");
    leg->AddEntry(h[4], "2f", "l");
  }
  if( counter[4] != 0 ){//4f, 2nu2mu
    h[5]->SetLineColor(kTeal+4);h[5]->SetMarkerColor(kTeal+4);
    h[5]->Draw("same HIST");
    leg->AddEntry(h[5], "4f(2#nu2#mu)", "l");
  }
  if( counter[5] != 0 ){//4f, 2nu2tau, tau to mu
    h[6]->SetLineColor(kTeal+4);h[6]->SetMarkerColor(kTeal+4);h[6]->SetLineStyle(2);
    h[6]->Draw("same HIST");
    leg->AddEntry(h[6], "4f(2#nu2#tau->#mu)", "l");
  }
  if( counter[6] != 0 ){//4f, 2nu1mu1tau, tau to mu
    h[7]->SetLineColor(kTeal+4);h[7]->SetMarkerColor(kTeal+4);h[7]->SetLineStyle(3);
    h[7]->Draw("same HIST");
    leg->AddEntry(h[7], "4f(2#nu1#mu1#tau->#mu)", "l");
  }
  if( counter[7] != 0 ){//4f, other
    h[8]->SetLineColor(kTeal+4);h[8]->SetMarkerColor(kTeal+4);h[8]->SetLineStyle(4);
    h[8]->Draw("same HIST");
    leg->AddEntry(h[8], "4f(other)", "l");
  }
  if( counter[8] != 0 ){//aa_4f, 2mu2nu
    h[9]->SetLineColor(kViolet);h[9]->SetMarkerColor(kViolet);
    h[9]->Draw("same HIST");
    leg->AddEntry(h[9], "#gamma#gamma->4f(2#nu2#mu)", "l");
  }
  if( counter[9] != 0 ){//aa_4f, 2tau2nu, tau to mu
    h[10]->SetLineColor(kViolet);h[10]->SetMarkerColor(kViolet);h[10]->SetLineStyle(2);
    h[10]->Draw("same HIST");
    leg->AddEntry(h[10], "#gamma#gamma->4f(2#nu2#tau->#mu)", "l");
  }
  if( counter[10] != 0 ){//aa_4f, 2nu1mu1tau, tau to mu
    h[11]->SetLineColor(kViolet);h[11]->SetMarkerColor(kViolet);h[11]->SetLineStyle(3);
    h[11]->Draw("same HIST");
    leg->AddEntry(h[11], "#gamma#gamma->4f(2#nu1#mu1#tau->#mu)", "l");
  }
  if( counter[11] != 0 ){//aa_4f, other
    h[12]->SetLineColor(kViolet);h[12]->SetMarkerColor(kViolet);h[12]->SetLineStyle(4);
    h[12]->Draw("same HIST");
    leg->AddEntry(h[12], "#gamma#gamma->4f(other)", "l");
  }
  //leg->Draw();

  //normalized to 1
  double norm_sig = h[1]->Integral();
  h[1]->Scale(1./norm_sig);
  h[1]->Draw("HIST");
  for( int i = 2; i <= 12; i++ )
  if( counter[i-1] != 0 ){
    double norm = h[i]->Integral();
    h[i]->Scale(1./norm);
    h[i]->Draw("same HIST");
  }
  leg->Draw();

  //show statistics
  cout << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other" << setw(12) << "2f"
       << setw(12) << "4f 2n2m" << setw(12) << "4f 2n2t"
       << setw(12) << "4f 2nmt" << setw(12) << "4f other"
       << setw(12) << "aa_4f 2n2m" << setw(12) << "aa_4f 2n2t"
       << setw(12) << "aa_4f 2nmt" << setw(12) << "aa_4f other" << endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
       << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5]
       << setw(12) << entryKind[6] << setw(12) << entryKind[7] << setw(12) << entryKind[8]
       << setw(12) << entryKind[9] << setw(12) << entryKind[10] << setw(12) << entryKind[11] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
       << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5]
       << setw(12) << counter[6] << setw(12) << counter[7] << setw(12) << counter[8]
       << setw(12) << counter[9] << setw(12) << counter[10] << setw(12) << counter[11] << endl;

  double bkg_tot = 0.;
  for( int i = 1; i <= 11; i++ ){
    bkg_tot += entryKind[i];
  }
  cout << "signal = " << entryKind[0] << ", background = " << bkg_tot << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + bkg_tot );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;
}
