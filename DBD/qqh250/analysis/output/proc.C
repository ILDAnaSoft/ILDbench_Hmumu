#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "TNtupleD.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

using namespace std;

void process()
{
  //gROOT->Reset();

  //TFile *file = new TFile("alldata.root","update");
  TFile *file = new TFile("alldata_select.root","update");

  TH1D *hev = new TH1D("hev","hev",300000,0,300000);

  TTree *tree = (TTree*)file->Get("dataTree");
  //char process[100];


  int processid;
  tree->SetBranchAddress("processid", &processid);
  for( int i = 0; i < tree->GetEntries(); i++ ){
    tree->GetEntry(i);
    //cerr << processid << endl;
    hev->Fill(processid-1);
    if(i%50000==0) cerr<<i<<" "<<processid<<endl;
  }



  file->Write();
  //hev->Draw();
  //hev->Print("proc.root");
  //tree->Write();
  //file->Close();
}
