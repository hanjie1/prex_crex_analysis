#include "SetCut.h"
#include "LoadACC.h"

void Accfunc(){
     TChain *TL = new TChain("T");
     TL->Add("../Rootfiles/PREXLHRSnom*");

     TChain *TR = new TChain("T");
     TR->Add("../Rootfiles/PREXRHRSnom*");

     // elastic peak position
     TH1F *hpztargL = new TH1F("hpztargL","LHRS p_ztarg distribution with cut",200,945,954);
     TL->Draw("p_ztarg>>hpztargL",(colCutL+isPb)*"rate");
     hpztargL->Draw("HIST");

     Int_t pbin = hpztargL->GetMaximumBin();
     Double_t p_peakL = hpztargL->GetBinCenter(pbin);
     cout<<"Left p peak:  "<<p_peakL<<endl; 

     // dp cut for the accepted events
     TString dpcutL = Form("(p_ztarg_tr-50*th_ztarg_tr)>(%f-2.2)",p_peakL);
     TCut DPL = Form("%s",dpcutL.Data());

     // beam E cut for the accepted events
     TString ecutL = Form("(ev.beamp*1000.-%f)>-5.0",p_peakL);
     TCut ECUTL = Form("%s",ecutL.Data());

     TCut ACCL = DPL+isPb+colCutL+XCUT+UPCutL+DownCutL;
     TCut INCL = isPb+ECUTL;

     // elastic peak position
     TH1F *hpztargR = new TH1F("hpztargR","RHRS p_ztarg distribution with cut",200,945,954);
     TR->Draw("p_ztarg>>hpztargR",(colCutR+isPb)*"rate");
     hpztargR->Draw("HIST");

     Int_t pbin_R = hpztargR->GetMaximumBin();
     Double_t p_peakR = hpztargR->GetBinCenter(pbin_R);
     cout<<"Right p peak:  "<<p_peakR<<endl; 

     // dp cut for the accepted events
     TString dpcutR = Form("(p_ztarg_tr-50*th_ztarg_tr)>(%f-2.2)",p_peakR);
     TCut DPR = Form("%s",dpcutR.Data());

     // beam E cut for the accepted events
     TString ecutR = Form("(ev.beamp*1000.-%f)>-5.0",p_peakR);
     TCut ECUTR = Form("%s",ecutR.Data());

     TCut ACCR = DPR+isPb+colCutR+XCUT+UPCutR+DownCutR;
     TCut INCR = isPb+ECUTR;

     const int nbin = 200;

     TH1F *hincL_v = new TH1F("hincL_v","LHRS incident events angle distribution",nbin,3,8); 
     TH1F *haccL_v = new TH1F("haccL_v","LHRS accepted events angle distribution",nbin,3,8); 

     hincL_v->Sumw2();
     haccL_v->Sumw2();

     TL->Draw("ev.Th>>hincL_v",INCL*"rate");
     TL->Draw("ev.Th>>haccL_v",ACCL*"rate");

     Double_t tot_incL = hincL_v->Integral();
     hincL_v->Scale(1./tot_incL);
     haccL_v->Scale(1./tot_incL);

     TH1F *hincR_v = new TH1F("hincR_v","RHRS incident events angle distribution",nbin,3,8); 
     TH1F *haccR_v = new TH1F("haccR_v","RHRS accepted events angle distribution",nbin,3,8); 

     hincR_v->Sumw2();
     haccR_v->Sumw2();

     TR->Draw("ev.Th>>hincR_v",INCR*"rate");
     TR->Draw("ev.Th>>haccR_v",ACCR*"rate");

     Double_t tot_incR = hincR_v->Integral();
     hincR_v->Scale(1./tot_incR);
     haccR_v->Scale(1./tot_incR);

     TH1F *hinc_v = (TH1F *)hincL_v->Clone("hinc_v"); 
     hinc_v->Add(hincR_v);

     TH1F *hacc_v = (TH1F *)haccL_v->Clone("hacc_v"); 
     hacc_v->Add(haccR_v);


     TGraphErrors *gACC_v = new TGraphErrors();
     TGraphErrors *gACC_v_norm = new TGraphErrors();

     int nn = 0;
     Double_t del_th = (8.0-3.0)/(1.0*nbin);
     Double_t inte=0;
     Double_t th_bin[nbin]={0}, acc[nbin]={0}, acc_err[nbin]={0};

     ofstream outfile;
     outfile.open("accfunction.csv");
     outfile<<"vertex angle,acceptance,stat_err"<<endl;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc_v->GetBinContent(ii);
	th_bin[ii-1] = hinc_v->GetBinCenter(ii);
	if(ninc==0) {
	   acc[ii-1]=0;
	   acc_err[ii-1]=0;
	   continue;
	}

	Double_t nacc = hacc_v->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;
 	acc[ii-1] = tmp_ratio;

	inte += tmp_ratio*del_th;

        Double_t inc_sumw2_sqrt = hinc_v->GetBinError(ii); // sqrt(sumw2)
	Double_t tmp_err2 = tmp_ratio * (1.-tmp_ratio) * pow(inc_sumw2_sqrt,2)/pow(ninc,2);
	Double_t tmp_err = sqrt(tmp_err2);

 	acc_err[ii-1] = tmp_err;

	gACC_v->SetPoint(nn,th_bin[ii-1],tmp_ratio);	
	gACC_v->SetPointError(nn,0,tmp_err);	

	outfile<<th_bin[ii-1]<<","<<tmp_ratio<<","<<tmp_err<<endl;
 	nn++;
     } 
     outfile.close();

     ofstream outfile1;
     outfile1.open("accfunction_norm.csv");
     outfile1<<"vertex angle,acceptance,stat_err"<<endl;
     Double_t acc_norm[nbin]={0},acc_norm_err[nbin]={0};
     for(int ii=0; ii<nbin; ii++){
	acc_norm[ii] = acc[ii]/inte;
	acc_norm_err[ii] = acc_err[ii]/inte;

	gACC_v_norm->SetPoint(ii,th_bin[ii],acc_norm[ii]);
	gACC_v_norm->SetPointError(ii,0,acc_norm_err[ii]);
	outfile1<<th_bin[ii]<<","<<acc_norm[ii]<<","<<acc_norm_err[ii]<<endl;
     }

     outfile1.close();

     ofstream outfile2;
     outfile2.open("accfunction_norm_1000.csv");
     outfile2<<"vertex angle,acceptance,stat_err"<<endl;

     double dtheta = (8.0-3.0)/1000.;  // delta theta in radius
     for(int ii=0; ii<1000; ii++){
	double tmp_th = 3.0+ii*dtheta;
	double tmp_ACC = FindACC(tmp_th,th_bin,acc_norm,acc_norm_err,nbin);
	outfile2<<tmp_th<<","<<tmp_ACC<<","<<0<<endl;
     }

     outfile2.close();
     


     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(4);
     gACC_v->Draw("AP");
     gACC_v->SetTitle("general acceptance function;ev.Th (deg)");

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hinc_v->Draw();
     hinc_v->SetLineColor(4);
     hinc_v->SetTitle("general incident angle distribution;ev.Th (deg);");

     c2->cd(2);
     hacc_v->Draw();
     hacc_v->SetLineColor(4);
     hacc_v->SetTitle(" general accepted angle distribution;ev.Th (deg);");

     TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
     gACC_v_norm->SetMarkerStyle(8);
     gACC_v_norm->SetMarkerColor(4);
     gACC_v_norm->Draw("AP");
     gACC_v_norm->SetTitle("general acceptance function (normalized);ev.Th (deg)");


}
