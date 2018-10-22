#include <iostream>
#include <stdlib.h>
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1.h"

const double epol_const = 0.1;
const double ppol_const = 0.65;
const double target_lumi = 1600.;

void AddWeight( const char *filename, TH1 *hev ){
  TFile *file = TFile::Open( filename, "update" );
  TTree *nt = (TTree*)file->Get("dataTree");
  double weight;
  int processid;
  float xsec, epol, ppol;
  TBranch *br = nt->Branch( "weight", &weight );
  nt->SetBranchAddress( "epol", &epol );
  nt->SetBranchAddress( "ppol", &ppol );
  nt->SetBranchAddress( "xsec", &xsec );
  nt->SetBranchAddress( "processid", &processid );

  for( int i = 0; i < nt->GetEntries(); i++ ){
    nt->GetEntry(i);
    if( i % 50000 == 0 ){ cout << i << " events processed." << endl; }
    weight = target_lumi * xsec / double( hev->GetBinContent( processid ) );
    weight *= ( epol == 0 ? 1.0 : epol == 1 ? ( 1.0 - epol_const ) : epol_const );
    weight *= ( ppol == 0 ? 1.0 : ppol == 1 ? ( 1.0 - ppol_const ) : ppol_const );
    br->Fill();
  }
  file->Write();
}
