void treeseparate(){
  //TFile *file = TFile::Open("alldata_selectvar.root","update");
  TFile *file = TFile::Open("alldata_precuts.root","update");
  //TFile *file = TFile::Open("alldata_basiccuts.root","update");
  TTree *t1 = (TTree*)file->Get("dataTree");
  cerr << "copying to ttrain..." << endl;
  TTree *ttrain = t1->CopyTree("Entry$%2==0");
  ttrain->SetName("datatrain");
  cerr << "copying to ttest..." << endl;
  TTree *ttest = t1->CopyTree("!(Entry$%2==0)");
  ttest->SetName("datatest");
  file->Write();
  cerr << "Done." << endl;
}
