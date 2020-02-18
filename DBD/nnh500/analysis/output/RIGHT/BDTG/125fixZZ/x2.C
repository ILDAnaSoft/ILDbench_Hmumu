#include "TFile.h"
#include "TNtupleD.h"
void x2(){
  TFile *file = TFile::Open( "alldata_after.root", "update" );
  TTree *nt = (TTree*)file->Get("datatest");
  double weight, weight2;
  TBranch *br = nt->Branch( "weight2", &weight2 );
  nt->SetBranchAddress( "weight", &weight );
  for( int i = 0; i < nt->GetEntries(); i++ ){
    nt->GetEntry(i);
    weight2 = 2.*weight;
    br->Fill();
  }
  file->Write();
}
