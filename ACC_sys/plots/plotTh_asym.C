void plotTh_asym(){
   ifstream infile;
   infile.open("acceptance_smearRMS_shiftTH.txt");

     if(!infile.is_open()){
        printf("Can't open file accfunction file\n");
        exit(0);
     }

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    Double_t smear[7]={0},asym_e[7]={0},asym_e_rms[7]={0},th_e[7]={0},th_e_rms[7]={0};

    while(tmp.ReadLine(infile)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          smear[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          asym_e[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          asym_e_rms[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          th_e[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          th_e_rms[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile.close();

    TGraph *g1 = new TGraph(7);
    for(int ii=0; ii<7; ii++){
	Double_t tmp_asym = asym_e_rms[ii]/asym_e[ii];
	Double_t tmp_th = th_e_rms[ii]/th_e[ii];

	g1->SetPoint(ii,smear[ii],tmp_asym/tmp_th);
    } 

    g1->SetMarkerStyle(8);	
    g1->Draw("AP");
    g1->SetTitle(";smear (%);(dA/A)/(dTh/Th)");
}
