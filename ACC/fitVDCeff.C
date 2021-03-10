void fitVDCeff(){
     Double_t rate[10]={45,110,200,510,70,355,242,91,53,15};
     Double_t eff[10]={0.96307,0.94469,0.932779,0.653148,0.960239,0.835,0.9236,0.9545,0.9794,0.972835};

     TGraph *geff = new TGraph(10,rate,eff);
     //geff->Fit("pol3");
     geff->SetMarkerStyle(8);
     geff->Draw("AP");
  
}
