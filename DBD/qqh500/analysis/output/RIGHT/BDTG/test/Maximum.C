void Maximum(){
  ifstream in("Analysis_result.dat");

  int i = 0;
  const int cases = 720;
  double significance[cases], Nsig[cases], Nbkg[cases], point[cases], shrinkage[cases];
  int cut[cases], depth[cases], tree[cases], node[cases];
  double Msignificance = 0, MNsig = 0, MNbkg = 0, Mpoint = 0, Mshrinkage = 0;
  int Mcut = 0, Mdepth = 0, Mtree = 0, Mnode = 0;
  while( in >> significance[i] >> Nsig[i] >> Nbkg[i] >> point[i] 
	    >> cut[i] >> shrinkage[i] >> depth[i] >> tree[i] >> node[i] ){
    if( significance[i] > Msignificance ){
      Msignificance = significance[i];
      MNsig = Nsig[i]; MNbkg = Nbkg[i]; Mpoint = point[i]; Mshrinkage = shrinkage[i];
      Mcut = cut[i]; Mdepth = depth[i]; Mtree = tree[i]; Mnode = node[i];
    }
    i++;
  }

  cout << "significance = " << Msignificance << ", Nsig = " << MNsig << ", Nbkg = " << MNbkg << ", at point = " << Mpoint << endl;
  cout << "cut = " << Mcut << ", shrinkage = " << Mshrinkage << ", depth = " << Mdepth << ", tree = " << Mtree << ", node = " << Mnode << endl;
}
