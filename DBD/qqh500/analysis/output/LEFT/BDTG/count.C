void count(){
  //TFile *file = TFile::Open("alldata_selectvar.root");
  TFile *file = TFile::Open("alldata_precuts.root");
  //TFile *file = TFile::Open("alldata_basiccuts.root");
  int processid;
  double weight;
  int type;

  TTree *nttrain = (TTree*)file->Get("datatrain");
  nttrain->SetBranchAddress( "processid", &processid );
  nttrain->SetBranchAddress( "weight", &weight );
  nttrain->SetBranchAddress( "type", &type );

  int Nsigtrain = 0;
  int Nbkgtrain = 0;
  double sigtrainweight = 0.;
  double bkgtrainweight = 0.;

  for( int i = 0; i < nttrain->GetEntries(); i++ ){
    nttrain->GetEntry(i);
    if( processid >= 108161 && processid <= 108164 &&
	( type >= 1 && type <= 5 ) ){
      sigtrainweight += weight;
      Nsigtrain++;
    }
    else{
      bkgtrainweight += weight;
      Nbkgtrain++;
    }
  }
  cout << "----------For training sample with precuts----------" << endl;
  cout << "Nsig = " << Nsigtrain << ", Nsig(weight) = " << sigtrainweight << endl;
  cout << "Nbkg = " << Nbkgtrain << ", Nbkg(weight) = " << bkgtrainweight << endl;
  cout << endl;

  TTree *nt = (TTree*)file->Get("datatest");
  nt->SetBranchAddress( "processid", &processid );
  nt->SetBranchAddress( "weight", &weight );
  nt->SetBranchAddress( "type", &type );

  int Nsig = 0;
  int Nbkg = 0;
  double sigweight = 0.;
  double bkgweight = 0.;

  for( int i = 0; i < nt->GetEntries(); i++ ){
    nt->GetEntry(i);
    if( processid >= 108161 && processid <= 108164 &&
	( type >= 1 && type <= 5 ) ){
      sigweight += weight;
      Nsig++;
    }
    else{
      bkgweight += weight;
      Nbkg++;
    }
  }
  cout << "----------For test sample with precuts----------" << endl;
  cout << "Nsig = " << Nsig << ", Nsig(weight) = " << sigweight << endl;
  cout << "Nbkg = " << Nbkg << ", Nbkg(weight) = " << bkgweight << endl;
  cout << endl;
}
