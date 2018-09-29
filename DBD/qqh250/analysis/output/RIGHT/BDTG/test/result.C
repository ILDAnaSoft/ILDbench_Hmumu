#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"

#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

using namespace std;

void result( const char *fname )
{
   // --- Create the Reader object
  vector<string> varsFinal;

  varsFinal.push_back("E_jj");
  varsFinal.push_back("costh_jj");
  //varsFinal.push_back("costh_Z");
  varsFinal.push_back("mumu_E");
  varsFinal.push_back("mumu_costh");
  varsFinal.push_back("mumu_costh_tobeam");
  varsFinal.push_back("sum_charge_costh");
  varsFinal.push_back("subleadingmu_E");
  varsFinal.push_back("leadingmu_costh");
  varsFinal.push_back("subleadingmu_costh");

  TMVA::Reader *readerFinal = new TMVA::Reader(varsFinal, "!Color:!Silent" );
  //TMVA::Reader *readerFinal2 = new TMVA::Reader(varsFinal2, "!Color:!Silent" );

   // --- Book the MVA methods
  readerFinal->BookMVA("mymethod","weights/TMVAClassification_BDTG.weights.xml");
  //readerFinal2->BookMVA("mymethod","../final2/weights/TMVAClassification_BDT.weights.xml");

   TFile* input = TFile::Open( fname,"update" );

   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   TTree* tree3 = (TTree*)input->Get("datatest");

   double BDTGoutput;
   //double lfinal2;
   TBranch *br = tree3->Branch("BDTGoutput",&BDTGoutput);
   //TBranch *br2 = tree3->Branch("lfinal2",&lfinal2);
   
   std::cout << "--- Processing: " << tree3->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   Int_t nEvent = tree3->GetEntries();
   for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
			tree3->GetEntry(ievt);

			vector<float> f;
			f.resize(varsFinal.size());
			for(unsigned int nvar=0;nvar<varsFinal.size();nvar++){
				TTreeFormula f1("n",varsFinal[nvar].c_str(),tree3);
				f[nvar] = f1.EvalInstance(ievt);
			}
			BDTGoutput = readerFinal->EvaluateMVA(f,"mymethod");
			/*
			f.resize(varsFinal2.size());
			for(unsigned int nvar=0;nvar<varsFinal2.size();nvar++){
				TTreeFormula f1("n",varsFinal2[nvar].c_str(),tree3);
				f[nvar] = f1.EvalInstance(ievt);
			}
			lfinal2 = readerFinal2->EvaluateMVA(f,"mymethod");
			*/
      br->Fill();
      //br2->Fill();
   }
   tree3->Write();
   input->Close();

   std::cout << "--- Created root file containing the MVA output" << std::endl;
  
   delete readerFinal;
    
   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}
