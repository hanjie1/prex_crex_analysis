#include "LoadAsym.h"
#include "LoadACC.h"
#include "SetCut.h"

void CheckVertexAsym(){
     double ppb = 1e9;
     double amu_c2 = 0.931494028; // GeV

     TChain *TL = new TChain("T");
     TL->Add("../Rootfiles/PREXLHRSnom*");

     TChain *TR = new TChain("T");
     TR->Add("../Rootfiles/PREXRHRSnom*");

     Double_t p_peakL = 952.312;

     // dp cut for the accepted events
     TString dpcutL = Form("(p_ztarg_tr-50*th_ztarg_tr)>(%f-2.2)",p_peakL);
     TCut DPL = Form("%s",dpcutL.Data());

     // beam E cut for the accepted events
     TString ecutL = Form("(ev.beamp*1000.-%f)>-5.0",p_peakL);
     TCut ECUTL = Form("%s",ecutL.Data());

     TCut ACCL = DPL+isPb+colCutL+XCUT+UPCutL+DownCutL;

     Double_t p_peakR = 952.312;

     // dp cut for the accepted events
     TString dpcutR = Form("(p_ztarg_tr-50*th_ztarg_tr)>(%f-2.2)",p_peakR);
     TCut DPR = Form("%s",dpcutR.Data());

     // beam E cut for the accepted events
     TString ecutR = Form("(ev.beamp*1000.-%f)>-5.0",p_peakR);
     TCut ECUTR = Form("%s",ecutR.Data());

     TCut ACCR = DPR+isPb+colCutR+XCUT+UPCutR+DownCutR;

     // load acceptance table
     const int nbin_acc=200;
     double accp_angle[nbin_acc]={0}, accp[nbin_acc]={0}, accp_err[nbin_acc]={0};
     TString acc_shape = "original";  // option: original, smear, box, gaus, shift

     int status = LoadACC("accfunction.csv",accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     bool check_ACC = true; // check mean and rms of acceptance function
     TH1F *hACC_check = new TH1F("hACC_check","acceptance function distribution",nbin_acc,3,8);
     if(check_ACC){
        int check_nbin_th = 100000;
        double dtheta = (8.0-3.0)/(check_nbin_th*1.0);  // delta theta in radius
        for(int ii=0; ii<check_nbin_th; ii++){
          double thisAngle = 3.0 + dtheta*ii;
          double thisACC = 0;
	  if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,nbin_acc);  // original acceptance function

          if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

          hACC_check->Fill(thisAngle,thisACC/(1.0*100000/(1.0*nbin_acc)));
        }

        Double_t mean = hACC_check->GetMean();
        Double_t RMS = hACC_check->GetRMS();
        Double_t mean_bin = hACC_check->FindBin(mean);
        Double_t peak = hACC_check->GetBinContent(mean_bin);
        cout<<acc_shape<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
     }

      
     // load asymmetry table
     LoadTable("horpb.dat", 0);

     int nbins=250;
     // LHRS
     TH1F *hasymL_origin = new TH1F("hasymL_origin","LHRS vertex asymmetry in g4hrs with ACC cut",nbins,200,1000);
     TH1F *hthL_origin = new TH1F("hthL_origin","LHRS vertex theta in g4hrs with ACC cut",150,3,8);
     TH1F *hq2L_origin = new TH1F("hq2L_origin","LHRS vertex Q2 in g4hrs with ACC cut",150,0,0.015);

     // RHRS
     TH1F *hasymR_origin = new TH1F("hasymR_origin","RHRS vertex asymmetry in g4hrs with ACC cut",nbins,200,1000);
     TH1F *hthR_origin = new TH1F("hthR_origin","RHRS vertex theta in g4hrs with ACC cut",150,3,8);
     TH1F *hq2R_origin = new TH1F("hq2R_origin","RHRS vertex Q2 in g4hrs with ACC cut",150,0,0.015);

     TL->Draw("ev.A>>hasymL_origin",ACCL*"rate");   // vertex asymmetry calculated in g4hrs with acceptance cuts applied
     TL->Draw("ev.Th>>hthL_origin",ACCL*"rate");   // vertex theta calculated in g4hrs with acceptance cuts applied
     TL->Draw("ev.Q2>>hq2L_origin",ACCL*"rate");   // vertex Q2 calculated in g4hrs with acceptance cuts applied

     TR->Draw("ev.A>>hasymR_origin",ACCR*"rate");   // vertex asymmetry calculated in g4hrs with acceptance cuts applied
     TR->Draw("ev.Th>>hthR_origin",ACCR*"rate");   // vertex theta calculated in g4hrs with acceptance cuts applied
     TR->Draw("ev.Q2>>hq2R_origin",ACCR*"rate");   // vertex Q2 calculated in g4hrs with acceptance cuts applied

     Double_t totL = hasymL_origin->Integral();
     hasymL_origin->Scale(1./totL);
     hthL_origin->Scale(1./totL);
     hq2L_origin->Scale(1./totL);

     Double_t totR = hasymR_origin->Integral();
     hasymR_origin->Scale(1./totR);
     hthR_origin->Scale(1./totR);
     hq2R_origin->Scale(1./totR);

     TH1F *hasym_origin = (TH1F *)hasymL_origin->Clone("hasym_origin"); 
     hasym_origin->Add(hasymR_origin);

     TH1F *hth_origin = (TH1F *)hthL_origin->Clone("hth_origin"); 
     hth_origin->Add(hthR_origin);

     TH1F *hq2_origin = (TH1F *)hq2L_origin->Clone("hq2_origin"); 
     hq2_origin->Add(hq2R_origin);

     // total
     TH1F *hasym_accL = new TH1F("hasym_accL","LHRS vertex asymmetry calculated from acceptance function with dp cut",nbins,200,1000);
     TH1F *hasym_accR = new TH1F("hasym_accR","RHRS vertex asymmetry calculated from acceptance function with dp cut",nbins,200,1000);

     double thisBeamE_L,thisTh_L,thisRate_L,p_ztarg_L,xfp_L;
     int thisA_L;
     TL->SetBranchAddress("ev.Th",&thisTh_L);
     TL->SetBranchAddress("ev.beamp",&thisBeamE_L);
     TL->SetBranchAddress("ev.nuclA",&thisA_L);
     TL->SetBranchAddress("rate",&thisRate_L);
     TL->SetBranchAddress("p_ztarg",&p_ztarg_L);
     TL->SetBranchAddress("x_fp_tr",&xfp_L);

     Double_t nentries = TL->GetEntries();
     for(int ii=0; ii<nentries; ii++){
        TL->GetEntry(ii);

	double thisAsym = Interpolate(thisBeamE_L,thisTh_L,0,1);  // vertex asymmetry (should be the same as ev.A)
	double thisACC = 0;
	if(acc_shape == "original" ) thisACC = FindACC(thisTh_L,accp_angle,accp,accp_err,nbin_acc);  // original acceptance function

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	if( thisA_L==208 && (thisBeamE_L*1000.-p_peakL)>-5) hasym_accL->Fill(thisAsym*ppb,thisRate_L*thisACC);
     }

     double thisBeamE_R,thisTh_R,thisRate_R,p_ztarg_R,xfp_R;
     int thisA_R;
     TR->SetBranchAddress("ev.Th",&thisTh_R);
     TR->SetBranchAddress("ev.beamp",&thisBeamE_R);
     TR->SetBranchAddress("ev.nuclA",&thisA_R);
     TR->SetBranchAddress("rate",&thisRate_R);
     TR->SetBranchAddress("p_ztarg",&p_ztarg_R);
     TR->SetBranchAddress("x_fp_tr",&xfp_R);

     nentries = TR->GetEntries();
     for(int ii=0; ii<nentries; ii++){
        TR->GetEntry(ii);

	double thisAsym = Interpolate(thisBeamE_R,thisTh_R,0,1);  // vertex asymmetry (should be the same as ev.A)
	double thisACC = 0;
	if(acc_shape == "original" ) thisACC = FindACC(thisTh_R,accp_angle,accp,accp_err,nbin_acc);  // original acceptance function

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	if( thisA_R==208 && (thisBeamE_R*1000.-p_peakR)>-5) hasym_accR->Fill(thisAsym*ppb,thisRate_R*thisACC);
     }

     Double_t totL_acc = hasym_accL->Integral();
     hasym_accL->Scale(1./totL_acc);

     Double_t totR_acc = hasym_accR->Integral();
     hasym_accR->Scale(1./totR_acc);

     TH1F *hasym_acc1 = (TH1F *)hasym_accL->Clone("hasym_acc1"); 
     hasym_acc1->Add(hasym_accR);

     TH1F *hasym_e_acc = new TH1F("hasym_e_acc","elemental asymmetry with acceptance function applied",nbins,200,1000);
     TH1F *hth_e_acc = new TH1F("hth_e_acc","elemental theta with acceptance function applied",150,3,8);
     TH1F *hq2_e_acc = new TH1F("hq2_e_acc","elemental Q2 with acceptance function applied",150,0,0.015);

     double setBeamE = 0.9534; 
     int nbin_th = 10000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius

     for(int ii=0; ii<nbin_th; ii++){
	double thisAngle = 3.0 + dtheta*ii;
	double thisAsym = Interpolate(setBeamE,thisAngle,0,1);  // elemental asymmetry
	double thisXS = Interpolate(setBeamE,thisAngle,0,0);  // elemental cross section

	double deg_to_rad = TMath::Pi()/180.;
	double thisEp = setBeamE/(1.+setBeamE/(208.*amu_c2)*(1.-TMath::Cos(thisAngle*deg_to_rad)));
	double thisQ2 = 4.*setBeamE*thisEp*pow(TMath::Sin(thisAngle*deg_to_rad/2.),2);

	double thisACC=0;
	if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,nbin_acc);  // original acceptance function

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	double omega = TMath::Sin(thisAngle*deg_to_rad)*dtheta*deg_to_rad;
	hasym_e_acc->Fill(thisAsym*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
	hth_e_acc->Fill(thisAngle,thisXS*omega*thisACC);   // elemetal angle with the acceptance function applied
	hq2_e_acc->Fill(thisQ2,thisXS*omega*thisACC);   // elemetal Q2 with the acceptance function applied
     } 

     cout<<"================================="<<endl;
     cout<<"orignal vertex asymmetry mean:     "<<setw(8)<<hasym_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_origin->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean:   "<<setw(8)<<hasym_acc1->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc1->GetRMS()<<endl;
     cout<<"elemental asymmetry mean:          "<<setw(8)<<hasym_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex angle mean:   "<<setw(8)<<hth_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_origin->GetRMS()<<endl;
     cout<<"elemental angle mean:        "<<setw(8)<<hth_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex Q2 mean:   "<<setw(8)<<hq2_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_origin->GetRMS()<<endl;
     cout<<"elemental Q2 mean:        "<<setw(8)<<hq2_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;

     Double_t inte1 = hasym_origin->Integral();
     Double_t inte2 = hasym_e_acc->Integral();
     Double_t inte3 = hasym_acc1->Integral();
     hasym_e_acc->Scale(inte1/inte2);
     //hasym_e_acc->Scale(inte1/inte2);

     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     hasym_origin->SetLineColor(1);
     hasym_acc1->SetLineColor(6);
     hasym_e_acc->SetLineColor(2);

     hasym_origin->SetLineWidth(2);
     hasym_acc1->SetLineWidth(2);
     hasym_e_acc->SetLineWidth(2);

     if(inte1>inte3){
       hasym_origin->Draw("HIST");    
       hasym_acc1->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_origin->SetTitle("Both arm asymmetry with acceptance; asym(ppb)");
     }
     else{
       hasym_acc1->Draw("HIST");    
       hasym_origin->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_acc1->SetTitle("Both arm asymmetry with acceptance; asym(ppb)");
     }

     TLegend *leg = new TLegend(0.15,0.65,0.35,0.8);
     leg->AddEntry(hasym_origin,"origin","L");
     //leg->AddEntry(hasym_acc1,Form("reproduce (dp_inc/p=%.3f)",delta_p_percent),"L");
     leg->AddEntry(hasym_acc1,Form("reproduce (beamE cut)"),"L");
     leg->AddEntry(hasym_e_acc,"elemental","L");
     leg->Draw();

   gPad->Update();
   TPaveStats* statA1 = (TPaveStats*)hasym_origin->FindObject("stats");
   TPaveStats* statA2 = (TPaveStats*)hasym_acc1->FindObject("stats");
   TPaveStats* statA3 = (TPaveStats*)hasym_e_acc->FindObject("stats");
   statA1->SetY1NDC(0.90);
   statA1->SetY2NDC(0.75);
   statA2->SetY1NDC(0.75);
   statA2->SetY2NDC(0.60);
   statA3->SetY1NDC(0.60);
   statA3->SetY2NDC(0.45);
   statA1->SetTextColor(1);
   statA2->SetTextColor(6);
   statA3->SetTextColor(2);
   gPad->Modified();

     TLatex tex1;
     tex1.SetTextAlign(11);
     tex1.SetTextSize(0.02);
     tex1.DrawLatexNDC(0.65,0.3,Form("%s acceptance",acc_shape.Data()));
    
     inte1 = hth_origin->Integral();
     inte2 = hth_e_acc->Integral();
     hth_e_acc->Scale(inte1/inte2);

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hth_origin->SetLineColor(1);
     hth_e_acc->SetLineColor(2);

     hth_origin->SetLineWidth(2);
     hth_e_acc->SetLineWidth(2);

     Double_t hpeak1 = hth_origin->GetMaximum();
     Double_t hpeak2 = hth_e_acc->GetMaximum();
  
     if(hpeak2<hpeak1){
       hth_origin->Draw("HIST");    
       hth_e_acc->Draw("HIST sames");
       hth_origin->SetTitle("Both arm angle with acceptance; angle");
     }
     else{
       hth_e_acc->Draw("HIST");
       hth_origin->Draw("HIST sames");    
       hth_e_acc->SetTitle("Both arm angle with acceptance; angle");
     }


     TLegend *leg1 = new TLegend(0.5,0.65,0.7,0.8);
     leg1->AddEntry(hth_origin,"origin","L");
     leg1->AddEntry(hth_e_acc,"elemental","L");
     leg1->Draw();

   gPad->Update();
   TPaveStats* statB1 = (TPaveStats*)hth_origin->FindObject("stats");
   TPaveStats* statB2 = (TPaveStats*)hth_e_acc->FindObject("stats");
   statB1->SetY1NDC(0.90);
   statB1->SetY2NDC(0.75);
   statB2->SetY1NDC(0.75);
   statB2->SetY2NDC(0.60);
   statB1->SetTextColor(1);
   statB2->SetTextColor(6);
   gPad->Modified();
    
     c2->cd(2);
     inte1 = hq2_origin->Integral();
     inte2 = hq2_e_acc->Integral();
     hq2_e_acc->Scale(inte1/inte2);

     hq2_origin->SetLineColor(1);
     hq2_e_acc->SetLineColor(2);

     hq2_origin->SetLineWidth(2);
     hq2_e_acc->SetLineWidth(2);

     hpeak1 = hq2_origin->GetMaximum();
     hpeak2 = hq2_e_acc->GetMaximum();

     if(hpeak1>hpeak2){
       hq2_origin->Draw("HIST");    
       hq2_e_acc->Draw("HIST sames");
       hq2_origin->SetTitle("Both arm Q2 with acceptance; Q2");
     }
     else{
       hq2_e_acc->Draw("HIST");
       hq2_origin->Draw("HIST sames");    
       hq2_e_acc->SetTitle("Both arm Q2 with acceptance; Q2");
     }


     TLegend *leg2 = new TLegend(0.15,0.65,0.35,0.8);
     leg2->AddEntry(hq2_origin,"origin","L");
     leg2->AddEntry(hq2_e_acc,"elemental","L");
     leg2->Draw();

   gPad->Update();
   TPaveStats* statC1 = (TPaveStats*)hq2_origin->FindObject("stats");
   TPaveStats* statC2 = (TPaveStats*)hq2_e_acc->FindObject("stats");
   statC1->SetY1NDC(0.90);
   statC1->SetY2NDC(0.75);
   statC2->SetY1NDC(0.75);
   statC2->SetY2NDC(0.60);
   statC1->SetTextColor(1);
   statC2->SetTextColor(6);
   gPad->Modified();
    
     TLatex tex2;
     tex2.SetTextAlign(11);
     tex2.SetTextSize(0.03);
     tex2.DrawLatexNDC(0.5,0.5,Form("%s acceptance",acc_shape.Data()));
  
     if(check_ACC){
       TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
       hACC_check->Draw();
     }

}
