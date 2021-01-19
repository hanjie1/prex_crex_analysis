#include "SetCut.h"

void plotCol(){

     TChain *T = new TChain("T");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_5.root");

     int nbins=100;
     TH2F *hcol = new TH2F("hcol","y vs x distribution at collimator with no cut",100,-0.15,0.15,100,-0.1,0.2);
     TH2F *hcol_dp = new TH2F("hcol_dp","y vs x distribution at collimator with dp cut",100,-0.15,0.15,100,-0.1,0.2);
     TH2F *hcol_acc = new TH2F("hcol_acc","y vs x distribution at collimator with acc cut",100,-0.15,0.15,100,-0.1,0.2);
     TH2F *hcol_all = new TH2F("hcol_all","y vs x distribution at collimator with all cut",100,-0.15,0.15,100,-0.1,0.2);

     TCut ACC = XCUT+colCut+isPb;

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl;

     // dp cut
     TCut DP = Form("(%f-p_ztarg)<2.2",p_peak); 
     
     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     c1->Divide(2,2);
     c1->cd(1);
     T->Draw("y_col_tr:x_col_tr>>hcol",isPb*"rate");

     c1->cd(2);
     T->Draw("y_col_tr:x_col_tr>>hcol_dp",(isPb+DP)*"rate");

     c1->cd(3);
     T->Draw("y_col_tr:x_col_tr>>hcol_acc",(isPb+colCut)*"rate");

     c1->cd(4);
     T->Draw("y_col_tr:x_col_tr>>hcol_all",(isPb+(!colCut))*"rate");
}
