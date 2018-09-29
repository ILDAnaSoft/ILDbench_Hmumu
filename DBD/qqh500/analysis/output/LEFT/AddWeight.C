#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "TString.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1.h"

const int MAX = 1000;
const double epol_const = 0.9;
const double ppol_const = 0.35;
const double target_lumi = 1600.;
int iprocid[MAX];
TString proc[MAX];
TString pol_c[MAX];
double epol[MAX];
double ppol[MAX];
double xsec[MAX];


void MakePolInfo( const char *infofile ){
  for( int i = 0; i < MAX; i++ ){
    iprocid[i] = 0;
    epol[i] = 0;
    ppol[i] = 0;
    xsec[i] = 0;
  }
  
  int i = 0;
  ifstream in( infofile );

  while( !in.eof() ){
    in >> iprocid[i] >> proc[i] >> pol_c[i] >> xsec[i];
    if( pol_c[i] == "eL.pL" ){
      epol[i] = -1;
      ppol[i] = -1;
    }
    else if( pol_c[i] == "eL.pW" || pol_c[i] == "eL.pB" ){
      epol[i] = -1;
      ppol[i] = 0;
    }
    else if( pol_c[i] == "eL.pR" ){
      epol[i] = -1;
      ppol[i] = 1;
    }
    else if( pol_c[i] == "eW.pL" || pol_c[i] == "eB.pL" ){
      epol[i] = 0;
      ppol[i] = -1;
    }
    else if( pol_c[i] == "eW.pR" || pol_c[i] == "eB.pR" ){
      epol[i] = 0;
      ppol[i] = 1;
    }
    else if( pol_c[i] == "eR.pL" ){
      epol[i] = 1;
      ppol[i] = -1;
    }
    else if( pol_c[i] == "eR.pW" || pol_c[i] == "eR.pB" ){
      epol[i] = 1;
      ppol[i] = 0;
    }
    else if( pol_c[i] == "eR.pR" ){
      epol[i] = 1;
      ppol[i] = 1;
    }
    else{
      epol[i] = 0;
      ppol[i] = 0;
    }
    i++;
  }
}

void AddWeight( const char *filename, TH1 *hev ){
  TFile *file = TFile::Open( filename, "update" );
  TTree *nt = (TTree*)file->Get("dataTree");
  double weight;
  int procid;
  TBranch *br = nt->Branch( "weight", &weight );
  nt->SetBranchAddress( "processid", &procid );

  for( int j = 0; j < nt->GetEntries(); j++ ){
    nt->GetEntry(j);
    if( j % 50000 == 0 ) cerr << j << " events processed." << endl;

    for( int k = 0; k < MAX; k++ ){
      if( procid == iprocid[k] ){
	weight = target_lumi * xsec[k] / double( hev->GetBinContent( procid ) );
	weight *= ( epol[k] == 0 ? 1.0 : epol[k] == 1 ? ( 1.0 - epol_const ) : epol_const );
	weight *= ( ppol[k] == 0 ? 1.0 : ppol[k] == 1 ? ( 1.0 - ppol_const ) : ppol_const );
      }
    }

    br->Fill();
  }
  file->Write();
}
