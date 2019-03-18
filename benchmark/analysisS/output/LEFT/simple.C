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
  TH1F *h[3];
  h[0] = new TH1F("all","ILD simulation",nbin,nlow,nhigh);
  h[1] = new TH1F("nnh_mumu","ILD simulation",nbin,nlow,nhigh);
  h[2] = new TH1F("bkg","bkg",nbin,nlow,nhigh);

  //int higgsdecay1;
  //nt->SetBranchAddress("higgsdecay1", &higgsdecay1);
  int type;
  nt->SetBranchAddress("type", &type);
  double weight;
  nt->SetBranchAddress("weight", &weight);
  int processid;
  nt->SetBranchAddress("processid", &processid);

  // stats                                                                
  double entryKind[5];
  memset(entryKind,0,sizeof(entryKind));
  int counter[5];
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
    else{//bkg
      counter[1]++;
      entryKind[1] += weight;
      h[2]->Fill(tfexp.EvalInstance(i),weight);
    }
  }

  // write histograms and stats
  TLegend *leg = new TLegend(0.7,0.7,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  //h[0]->SetLineColor(kBlack);h[0]->Draw("HIST");//ALL
  //leg->AddEntry(h[0], "All", "l");
  if( counter[0] != 0 ){//nnh, h->mumu
    h[1]->SetLineColor(kBlue);h[1]->SetMarkerColor(kBlue);h[1]->SetLineWidth(3);
    h[1]->Draw("same HIST");
    leg->AddEntry(h[1], "signal", "l");
  }
  if( counter[1] != 0 ){//bkg
    h[2]->SetLineColor(kRed);h[2]->SetMarkerColor(kRed);
    h[2]->Draw("same HIST");
    leg->AddEntry(h[2], "background", "l");
  }
  leg->Draw();

  //show statistics
  cout << setw(12) << "nnh hmumu" << setw(12) << "bkg" << endl;
  cout << setw(12) << entryKind[0] << setw(12) << entryKind[1] << endl;
  cout << setw(12) << counter[0] << setw(12) << counter[1] << endl;

  double bkg_tot = 0.;
  for( int i = 1; i <= 1; i++ ){
    bkg_tot += entryKind[i];
  }
  cout << "signal = " << entryKind[0] << ", background = " << bkg_tot << endl;
  double signi = entryKind[0] / TMath::Sqrt( entryKind[0] + bkg_tot );
  cout << "significance = " << signi << endl;
  cout << "precision = " << 100./signi << "%" << endl;

}
