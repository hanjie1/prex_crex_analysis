#include "SetCut.h"

void thetaEp_relation(){

     TChain *T = new TChain("T");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_5.root");

     int nbins=100;
     TH1F *hth_e1 = new TH1F("hth_e1","vertex theta with energy cut 1",nbins,3,8);
     TH1F *hth_e2 = new TH1F("hth_e2","vertex theta with energy cut 2",nbins,3,8);
     TH1F *hth_noe = new TH1F("hth_noe","vertex theta without energy cut",nbins,3,8);

     TH1F *hth_noacc_e1 = new TH1F("hth_noacc_e1","vertex theta with energy cut 1 no acc cut",nbins,3,8);
     TH1F *hth_noacc_e2 = new TH1F("hth_noacc_e2","vertex theta with energy cut 2 no acc cut",nbins,3,8);
     TH1F *hth_noacc_noe = new TH1F("hth_noacc_noe","vertex theta without energy cut no acc cut",nbins,3,8);

     TCut ACC = XCUT+colCut+isPb;

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl;

     // dp cut
     TCut DP1 = Form("(%f-p_ztarg)<2.2",p_peak); 
     TCut DP2 = Form("(%f-p_ztarg)<1.2",p_peak); 

     T->Draw("ev.Th>>hth_e1",(ACC+DP1)*"rate");
     T->Draw("ev.Th>>hth_e2",(ACC+DP2)*"rate");
     T->Draw("ev.Th>>hth_noe",ACC*"rate");

     T->Draw("ev.Th>>hth_noacc_e1",(DP1)*"rate");
     T->Draw("ev.Th>>hth_noacc_e2",(DP2)*"rate");
     T->Draw("ev.Th>>hth_noacc_noe","rate");

     Double_t peak1 = hth_noe->GetMaximum();
     Double_t peak2 = hth_e1->GetMaximum();
     Double_t peak3 = hth_e2->GetMaximum();
     hth_e1->Scale(peak1/peak2);
     hth_e2->Scale(peak1/peak3);

     peak1 = hth_noacc_noe->GetMaximum();
     peak2 = hth_noacc_e1->GetMaximum();
     peak3 = hth_noacc_e2->GetMaximum();
     hth_noacc_e1->Scale(peak1/peak2);
     hth_noacc_e2->Scale(peak1/peak3);

     TH1F *hth_ratio = (TH1F *)hth_e1->Clone("hth_ratio");
     hth_ratio->Divide(hth_noe);
    
     TH1F *hth_noacc_ratio = (TH1F *)hth_noacc_e1->Clone("hth_noacc_ratio");
     hth_noacc_ratio->Divide(hth_noacc_noe);

     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     c1->Divide(2,2);
     c1->cd(1);
     hth_noe->SetLineColor(4);
     hth_noe->SetLineWidth(2);
     hth_e1->SetLineColor(2);
     hth_e1->SetLineWidth(2);
     //hth_e2->SetLineColor(1);
     //hth_e2->SetLineWidth(2);

     hth_noe->Draw("HIST");
     hth_e1->Draw("HIST same");
     //hth_e2->Draw("HIST same");

    TLegend *leg = new TLegend(0.65,0.7,0.8,0.85);
    leg->AddEntry(hth_noe,"no dp cut"); 
    leg->AddEntry(hth_e1,"dp<2.2 MeV"); 
    //leg->AddEntry(hth_e2,"dp<1.2 MeV"); 
    leg->Draw();

     c1->cd(2);
     
     hth_noacc_noe->SetLineColor(4);
     hth_noacc_noe->SetLineWidth(2);
     hth_noacc_e1->SetLineColor(2);
     hth_noacc_e1->SetLineWidth(2);
     //hth_noacc_e2->SetLineColor(1);
     //hth_noacc_e2->SetLineWidth(2);

     hth_noacc_noe->Draw("HIST");
     hth_noacc_e1->Draw("HIST same");
     //hth_noacc_e2->Draw("HIST same");

    TLegend *leg1 = new TLegend(0.65,0.7,0.8,0.85);
    leg1->AddEntry(hth_noe,"no dp cut and no acc cut"); 
    leg1->AddEntry(hth_e1,"dp<2.2 MeV"); 
 //   leg1->AddEntry(hth_e2,"dp<1.2 MeV"); 
    leg1->Draw();

    c1->cd(3);
    hth_ratio->SetMarkerStyle(8);
    hth_ratio->Draw(); 

    c1->cd(4);
    hth_noacc_ratio->SetMarkerStyle(8);
    hth_noacc_ratio->Draw(); 
}
