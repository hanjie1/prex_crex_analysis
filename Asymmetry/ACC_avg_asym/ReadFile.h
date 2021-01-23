void ReadACC(TString filename, Double_t angle[], Double_t accp[]){
    ifstream file;
    file.open(filename);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,",");
          angle[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,",");
          accp[nn-1]=atof(content.Data());

          //cout<<x[kin][nn-1]<<"  "<<xavg[kin][nn-1]<<"  "<<Q2[kin][nn-1]<<"  "<<Yield[kin][nn-1]<<"  "<<Y_err[kin][nn-1]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return;
}

void ReadACC_Kent(TString filename, Double_t angle[], Double_t accp[]){
    ifstream file;
    file.open(filename);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,",");
          angle[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,",");
          tmp.Tokenize(content,from,",");
          tmp.Tokenize(content,from,",");
          accp[nn-1]=atof(content.Data());

          //cout<<x[kin][nn-1]<<"  "<<xavg[kin][nn-1]<<"  "<<Q2[kin][nn-1]<<"  "<<Yield[kin][nn-1]<<"  "<<Y_err[kin][nn-1]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return;
}

void ReadAsym(TString filename, Double_t angle[], Double_t asym[], Double_t xs[]){
    ifstream file;
    file.open(filename);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          if(nn==1){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from,",");
          angle[nn-2]=atof(content.Data());
          tmp.Tokenize(content,from,",");
          tmp.Tokenize(content,from,",   ");
          xs[nn-2]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          asym[nn-2]=atof(content.Data());

          //cout<<content<<endl;
          //cout<<angle[nn-2]<<"  "<<asym[nn-2]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return;
}
