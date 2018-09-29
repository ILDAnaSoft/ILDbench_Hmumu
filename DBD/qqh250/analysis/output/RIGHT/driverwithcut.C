#include <vector>
#include <fstream>
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

#include "stackwithcut.C"

int nstart= 0;
int nend = -1; // -1 means run till end;

vector<TCut> cuts;

// setup cuts to be applied sequentially
void setupCuts() {
  /*
  //precuts
  cuts.push_back( TCut( "n_muminus==1&&n_muplus==1" ) );
  cuts.push_back( TCut( "muminus_Chi2Ndf>0.5&&muminus_Chi2Ndf<1.5&&muplus_Chi2Ndf>0.5&&muplus_Chi2Ndf<1.5" ) );
  cuts.push_back( TCut( "abs(muminus_d0)<0.02&&abs(muplus_d0)<0.02&&abs(muminus_d0-muplus_d0)<0.02" ) );
  cuts.push_back( TCut( "abs(muminus_z0)<0.5&&abs(muplus_z0)<0.5&&abs(muminus_z0-muplus_z0)<0.5" ) );
  cuts.push_back( TCut( "sigma_mumu_mass<0.5" ) );
  cuts.push_back( TCut( "mumu_mass>100&&mumu_mass<130" ) );
  cuts.push_back( TCut( "mumu_costh<-0.45" ) );
  cuts.push_back( TCut( "n_forveto==0" ) );
  cuts.push_back( TCut( "nJets==2" ) );
  cuts.push_back( TCut( "jet1_charged>=4&&jet2_charged>=4" ) );
  cuts.push_back( TCut( "mass_jj>60&&mass_jj<135" ) );
  cuts.push_back( TCut( "recoilmass>70&&recoilmass<130" ) );
  */

  cuts.push_back( TCut( "n_muminus==1&&n_muplus==1" ) );
  cuts.push_back( TCut( "muminus_Chi2Ndf>0.5&&muminus_Chi2Ndf<1.5&&muplus_Chi2Ndf>0.5&&muplus_Chi2Ndf<1.5" ) );
  cuts.push_back( TCut( "abs(muminus_d0)<0.01&&abs(muplus_d0)<0.01&&abs(muminus_d0-muplus_d0)<0.01" ) );
  cuts.push_back( TCut( "abs(muminus_z0)<0.5&&abs(muplus_z0)<0.5&&abs(muminus_z0-muplus_z0)<0.5" ) );
  cuts.push_back( TCut( "sigma_mumu_mass<0.5" ) );
  cuts.push_back( TCut( "mumu_mass>100&&mumu_mass<130" ) );
  cuts.push_back( TCut( "mumu_costh<-0.45" ) );
  cuts.push_back( TCut( "n_forveto==0" ) );
  cuts.push_back( TCut( "nJets==2" ) );
  cuts.push_back( TCut( "mass_jj>60&&mass_jj<135" ) );
  cuts.push_back( TCut( "jet1_charged>=4&&jet2_charged>=4" ) );

  /*
  //optimization
  cuts.push_back( TCut( "mumu_mass>124&&mumu_mass<126" ) );
  cuts.push_back( TCut( "esum<265" ) );
  cuts.push_back( TCut( "Pt_all>60" ) );
  cuts.push_back( TCut( "principalthrust<0.88" ) );
  cuts.push_back( TCut( "muminus_charge_costh>-0.8&&muplus_charge_costh>-0.8" ) );
  //cuts.push_back( TCut( "abs(costh_thrustaxis)<0.96" ) );
  //cuts.push_back( TCut( "abs(costh_missmom)<0.98" ) );
  //cuts.push_back( TCut( "mumu_costh<0.56" ) );
  */
  
  //parameter search
  //cuts.push_back( TCut( "mumu_mass>124&&mumu_mass<126&&esum<265&&Pt_all>60&&principalthrust<0.88&&muplus_charge_costh>-0.8&&muminus_charge_costh>-0.8" ) );
  /*
  cuts.push_back( TCut( "esum<395" ) );
  cuts.push_back( TCut( "esum<390" ) );
  cuts.push_back( TCut( "esum<385" ) );
  cuts.push_back( TCut( "esum<380" ) );
  cuts.push_back( TCut( "esum<375" ) );
  cuts.push_back( TCut( "esum<370" ) );
  cuts.push_back( TCut( "esum<365" ) );
  cuts.push_back( TCut( "esum<360" ) );
  cuts.push_back( TCut( "esum<355" ) );
  cuts.push_back( TCut( "esum<350" ) );
  cuts.push_back( TCut( "esum<345" ) );
  cuts.push_back( TCut( "esum<340" ) );
  cuts.push_back( TCut( "esum<335" ) );
  cuts.push_back( TCut( "esum<330" ) );
  cuts.push_back( TCut( "esum<325" ) );
  cuts.push_back( TCut( "esum<320" ) );
  cuts.push_back( TCut( "esum<315" ) );
  cuts.push_back( TCut( "esum<310" ) );
  cuts.push_back( TCut( "esum<305" ) );
  cuts.push_back( TCut( "esum<300" ) );
  cuts.push_back( TCut( "esum<295" ) );
  cuts.push_back( TCut( "esum<290" ) );
  cuts.push_back( TCut( "esum<285" ) );
  cuts.push_back( TCut( "esum<280" ) );
  cuts.push_back( TCut( "esum<275" ) );
  cuts.push_back( TCut( "esum<270" ) );
  cuts.push_back( TCut( "esum<265" ) );
  cuts.push_back( TCut( "esum<260" ) );
  cuts.push_back( TCut( "esum<255" ) );
  cuts.push_back( TCut( "esum<250" ) );
  */
}


void setup() {
  //MakeProcinfo("procid500.txt");
  MakeProcinfo("procid250.txt");
}

void driver(
	    const char* output, const char* output2,
	    const char* filename, const char* ntpname,
	    int nbin, double histMin, double histMax,
	    const char* var ) {
	setup();
	setupCuts();

	ofstream ofile(output);
	ofstream ofile2(output2);

	TCanvas* c1 = new TCanvas("c1","",800,450);

  TFile* file = new TFile(filename);
	if (file == 0) {
		cout << "Could not open file " << filename << endl;
		return;
	}

	TTree* ntp = (TTree*) file->Get(ntpname);
	if (ntp == 0) {
		cout << "Could not get TTree " << ntpname << endl;
		return;
	}
	/*
	TH1* hev = (TH1*) file->Get("hev");
	if (hev == 0) {
		cout << "Could not get hev" << endl;
		return;
	}
	*/
	if (nend==-1) nend = cuts.size();

	//MakeAllWithWeight( ofile, ntp, hev, nbin, histMin, histMax, var, "1" );
	MakeAllWithWeight( ofile, ofile2, ntp, nbin, histMin, histMax, var, "1" );
	for (int i=nstart; i<nend; ++i) {
		TCut cut;
		for (int j=0; j<=i; ++j) {
			cut += cuts[j];
		}
		//MakeAllWithWeight( ofile, ntp, hev, nbin, histMin, histMax, var, cut );
		MakeAllWithWeight( ofile, ofile2, ntp, nbin, histMin, histMax, var, cut );

		// TODO make histogram
	}

	ofile.close();
}
