/**

 SkimVar.C

 Author: Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>

 Copies a TTree from a root file, extracting only the specified variables.
 This will reduce the file size to help with the analysis.

 Usage:
   SkimVar(input, output, ntpname, varlist)
     input       = input file name
     output      = output file name
     ntpname     = name of TTree object
     varlist     = list of variables to keep, delimited by a semicolon
                   For example: "var1:var2:var3"

 Note that SkimVar.C can be called directly from the command line.
 For example:

   root -b -q SkimVar.C+\(\"input.root\",\"output.root\",\"dataTree\",\"var1:var2:var3\"\)

 The extra backslashes "\" are needed for CINT to process the arguments.
 It is recommended to compile the macro for faster processing.

*/

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>

void SkimVar(const char* input, const char* output, const char* ntpname, const char* varlist) {

  TObjArray* arr = TString(varlist).Tokenize(":");

  cout << "Input file  : " << input   << endl;
  cout << "Output file : " << output  << endl;
  cout << "TTree name  : " << ntpname << endl;
  cout << "Variables   : " << varlist << endl;

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

  ntp->SetBranchStatus("*",0);

  for (Int_t i=0; i<arr->GetEntries(); ++i) {
    TObjString* p = (TObjString*) arr->At(i);
    cout << " enabling branch: " << p->String() << endl;
    ntp->SetBranchStatus(p->String(),1);
  }

  TFile* fo = new TFile(output,"RECREATE");
  if (fo == 0) {
    cerr << "Could not open output file:" << output << endl;
    cerr << "Aborting..." << endl;
    return;
  }

  cout << "Cloning TTree..." << endl;

  TTree* clone = ntp->CloneTree();
  clone->Write();

  fo->Close();

  cout << "Done." << endl;
}
