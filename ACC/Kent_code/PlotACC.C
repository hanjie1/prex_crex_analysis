#include "LoadACC.h"

void PlotACC(){
     double accp_angle[200]={0}, accp[200]={0}, accp_err[200]={0};
     int status = LoadACC("accfunction_norm.csv",accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     double accp_angle_k[100]={0}, accp_k[100]={0},accp_err_k[100]={0};
     status = LoadACC_kent("accfunc_central_average.csv",accp_angle_k, accp_k);

     double accp_angle_1[1000]={0}, accp_1[1000]={0}, accp_err_1[1000]={0};
     status = LoadACC("accfunction_norm_1000.csv",accp_angle_1, accp_1, accp_err_1);
    
     TGraph *gACC = new TGraph(1000,accp_angle_1,accp_1);
 
     TH1F *hACC = new TH1F("hACC","acceptance function",200,3,8);
     TH1F *hACC_k = new TH1F("hACC_k","kent's acceptance function",100,3,8);

     int nbin_th = 100000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius
     for(int ii=0; ii<nbin_th; ii++){
        double thisAngle = 3.0 + dtheta*ii;
        double thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,200);  // acceptance to the corresponding angle
        double thisACC_k = FindACC(thisAngle,accp_angle_k,accp_k,accp_err_k,100);  // acceptance to the corresponding angle
        if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}
        if(thisACC_k<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}
        
        hACC->Fill(thisAngle,thisACC/500.);
        hACC_k->Fill(thisAngle,thisACC_k/1000.);
     }
   
     Double_t mean = hACC->GetMean();
     Double_t RMS = hACC->GetRMS();
     Double_t mean_bin = hACC->FindBin(mean);
     Double_t peak = hACC->GetBinContent(mean_bin);
     cout<<"mean, rms"<<mean<<"    "<<RMS<<endl;

     mean = hACC_k->GetMean();
     RMS = hACC_k->GetRMS();
     mean_bin = hACC_k->FindBin(mean);
     peak = hACC_k->GetBinContent(mean_bin);
     cout<<"mean, rms"<<mean<<"    "<<RMS<<endl;

     TCanvas *c1 = new TCanvas("c1","c1",1000,1000); 
     c1->Divide(2,1);
     c1->cd(1);
     hACC->SetLineColor(2);
     hACC_k->SetLineColor(4);
     hACC->SetLineWidth(2);
     hACC_k->SetLineWidth(2);
     gACC->SetMarkerStyle(8);
     gACC->SetMarkerColor(5);

     hACC->Draw("HIST"); 
     hACC_k->Draw("HIST sames"); 

     
   gPad->Update();
   TPaveStats* statA1 = (TPaveStats*)hACC->FindObject("stats");
   TPaveStats* statA2 = (TPaveStats*)hACC_k->FindObject("stats");
   statA1->SetY1NDC(0.90);
   statA1->SetY2NDC(0.75);
   statA2->SetY1NDC(0.75);
   statA2->SetY2NDC(0.60);
   statA1->SetTextColor(1);
   statA2->SetTextColor(6);
   gPad->Modified();

   c1->cd(2);
   gACC->Draw("AP");

}
