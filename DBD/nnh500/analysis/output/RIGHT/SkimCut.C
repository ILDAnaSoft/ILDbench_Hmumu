/**

 SkimCut.C

 Author: Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>

 Copies a TTree from a root file, extracting only the events that passes a cut.
 This will reduce the file size to help with the analysis.

 Usage:
   SkimCut(input, output, ntpname, varlist)
     input       = input file name
     output      = output file name
     ntpname     = name of TTree object
     cut         = cut definition

 Note that SkimCut.C can be called directly from the command line.
 For example:

   root -b -q SkimCut.C+\(\"input.root\",\"output.root\",\"dataTree\",\"var1>0\"\)

 The extra backslashes "\" are needed for CINT to process the arguments.
 It is recommended to compile the macro for faster processing.

*/

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>

void SkimCut(const char* input, const char* output, const char* ntpname, const char* cut) {

  cout << "Input file  : " << input   << endl;
  cout << "Output file : " << output  << endl;
  cout << "TTree name  : " << ntpname << endl;
  cout << "Cut         : " << cut     << endl;

  TFile* fi = new TFile(input);
  if (fi == 0) {
    cerr << "Could not open input file:" << input << endl;
    cerr << "Aborting..." << endl;
    return;
  }

  TTree* ntp = (TTree*) fi->Get(ntpname);
  if (ntp == 0) {
    cerr << "Could not find TTree object: " << ntpname << endl;
    cerr << "Aborting..." << endl;
    return;
  }

  TFile* fo = new TFile(output,"RECREATE");
  if (fo == 0) {
    cerr << "Could not open output file:" << output << endl;
    cerr << "Aborting..." << endl;
    return;
  }

  cout << "Cloning TTree..." << endl;

  TTree* clone = ntp->CopyTree(cut);
  clone->Write();

  fo->Close();

  cout << "Done." << endl;
}
