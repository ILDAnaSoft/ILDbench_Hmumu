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

void MakeAllWithWeight( TString file1, TString file2, TTree *nt,
			int nbin, double nlow, double nhigh,
			const char *varexp, const char *selection )
{
  if( selection == NULL ) selection = "1"; // select all events

  // Make TreeFormulae
  TTreeFormula tfexp( "tfexp", varexp, nt );
  TTreeFormula tfsel( "tfsel", selection, nt );

  //Make histograms
  TH1F *h[7];
  h[0] = new TH1F("all","all",nbin,nlow,nhigh);
  h[1] = new TH1F("nnh_mumu","nnh_mumu",nbin,nlow,nhigh);
  h[2] = new TH1F("other_h_mumu","other_h_mumu",nbin,nlow,nhigh);
  h[3] = new TH1F("other_h_decay","other_h_decay",nbin,nlow,nhigh);
  h[4] = new TH1F("2f","2f",nbin,nlow,nhigh);
  h[5] = new TH1F("4f","4f",nbin,nlow,nhigh);
  h[6] = new TH1F("aa_4f","aa_4f",nbin,nlow,nhigh);

  //int higgsdecay1;
  //nt->SetBranchAddress("higgsdecay1", &higgsdecay1);
  int type;
  nt->SetBranchAddress("type", &type);
  double weight;
  nt->SetBranchAddress("weight", &weight);
  int processid;
  nt->SetBranchAddress("processid", &processid);

  // stats                                                                
  double entryKind[10];
  memset(entryKind,0,sizeof(entryKind));
  int counter[10];
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
    else if( ( processid >= 250014 && processid <= 250060 ) ){//4f
      counter[4]++;
      entryKind[4] += weight;
      h[5]->Fill(tfexp.EvalInstance(i),weight);
    }
    else if( ( processid >= 37385 && processid <= 37464 ) ){//aa_4f
      counter[5]++;
      entryKind[5] += weight;
      h[6]->Fill(tfexp.EvalInstance(i),weight);
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
  if( counter[4] != 0 ){//4f
    h[5]->SetLineColor(kTeal+4);h[5]->SetMarkerColor(kTeal+4);
    h[5]->Draw("same HIST");
    leg->AddEntry(h[5], "4f", "l");
  }
  if( counter[5] != 0 ){//aa_4f
    h[6]->SetLineColor(kTeal+4);h[6]->SetMarkerColor(kTeal+4);h[6]->SetLineStyle(2);
    h[6]->Draw("same HIST");
    leg->AddEntry(h[6], "#gamma#gamma->4f", "l");
  }
  leg->Draw();

  //show statistics
  cout << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
       << setw(12) << "2f" << setw(12) << "4f" << setw(12) << "aa_4f" <<  endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
       << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
       << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5] << endl;

  ofstream ofs(file1);
  ofs << setw(12) << "nnh hmumu" << setw(12) << "bkg hmumu" << setw(12) << "ffh other"
      << setw(12) << "2f" << setw(12) << "4f" << setw(12) << "aa_4f" <<  endl;
  ofs << setw(12) << entryKind[0] << setw(12) << entryKind[1] << setw(12) << entryKind[2]
      << setw(12) << entryKind[3] << setw(12) << entryKind[4] << setw(12) << entryKind[5] << endl;
  ofs << setw(12) << counter[0] << setw(12) << counter[1] << setw(12) << counter[2]
      << setw(12) << counter[3] << setw(12) << counter[4] << setw(12) << counter[5] << endl;

  double bkg_tot = 0.;
  for( int i = 1; i <= 5; i++ ){
    bkg_tot += entryKind[i];
  }
  cout << "signal = " << entryKind[0] << ", background = " << bkg_tot << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + bkg_tot );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;
  ofstream ofs2(file2);
  ofs2 << signi << " " << 100./signi << "%" << endl;

}
