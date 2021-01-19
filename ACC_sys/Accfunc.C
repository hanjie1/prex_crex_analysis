#include "SetCut.h"

void Accfunc(){
     TChain *T = new TChain("T");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_5.root");

     // elastic peak position
     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl; 

     // dp cut for the accepted events
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

     // beam E cut for the accepted events
     TString ecut = Form("abs(ev.beamp*1000.-%f)<2.2",p_peak);
     TCut ECUT = Form("%s",ecut.Data());

     // dp cut for the incident events
     double delta_p_percent = 2./100.;
     double delta_p = p_peak*delta_p_percent;
     TString dpcut_inc = Form("(%f-p_ztarg)<%f",p_peak,delta_p);
     TCut DP_INC = Form("%s",dpcut_inc.Data());

     TCut ACC = DP+isPb+colCut+XCUT+UPCut+DownCut;
     TCut INC = isPb+ECUT;

     const int nbin = 100;

     TH1F *hinc_v = new TH1F("hinc_v","incident events angle distribution",nbin,3,8); 
     TH1F *hacc_v = new TH1F("hacc_v","accepted events angle distribution",nbin,3,8); 

     hinc_v->Sumw2();
     hacc_v->Sumw2();

     T->Draw("ev.Th>>hinc_v",INC*"rate");
     T->Draw("ev.Th>>hacc_v",ACC*"rate");

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

     for(int ii=0; ii<nbin; ii++){
	Double_t acc_norm = acc[ii]/inte;
	Double_t acc_norm_err = acc_err[ii]/inte;

	gACC_v_norm->SetPoint(ii,th_bin[ii],acc_norm);
	gACC_v_norm->SetPointError(ii,0,acc_norm_err);
	outfile1<<th_bin[ii]<<","<<acc_norm<<","<<acc_norm_err<<endl;
     }

     outfile1.close();

     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     c1->Draw();
     TPad *p1 = new TPad("p1","p1",0.0,0.0,0.8,1.);
     p1->SetRightMargin(0.1);
     p1->Draw();
     p1->cd();

     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(4);
     gACC_v->Draw("AP");
     gACC_v->SetTitle("LHRS acceptance function;ev.Th (deg)");

     c1->cd(0);
     TPad *p2 = new TPad("p2","p2",0.8,0.0,1.,1.);
     p2->SetLeftMargin(0.);
     p2->Draw();
     p2->cd();
     TLatex tex;
     tex.SetTextAlign(11);
     tex.SetTextSize(0.05);
     tex.DrawLatexNDC(0.0,0.8,"thrown-in cut:");
     tex.DrawLatexNDC(0.1,0.75,ispb);
     tex.DrawLatexNDC(0.1,0.7,ecut);
     tex.DrawLatexNDC(0.0,0.65,"accepted cut:");
     tex.DrawLatexNDC(0.1,0.6,ispb);
     tex.DrawLatexNDC(0.1,0.55,colcut);
     tex.DrawLatexNDC(0.1,0.5,dpcut);
     tex.DrawLatexNDC(0.1,0.45,xcut);
     tex.DrawLatexNDC(0.1,0.4,"upcut");
     tex.DrawLatexNDC(0.1,0.35,"downcut");

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hinc_v->Draw();
     hinc_v->SetLineColor(4);
     hinc_v->SetTitle("LHRS incident angle distribution;ev.Th (deg);");

     c2->cd(2);
     hacc_v->Draw();
     hacc_v->SetLineColor(4);
     hacc_v->SetTitle(" LHRS accepted angle distribution;ev.Th (deg);");

     TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
     c3->Draw();
     TPad *p3 = new TPad("p3","p3",0.0,0.0,0.8,1.);
     p3->SetRightMargin(0.1);
     p3->Draw();
     p3->cd();

     gACC_v_norm->SetMarkerStyle(8);
     gACC_v_norm->SetMarkerColor(4);
     gACC_v_norm->Draw("AP");
     gACC_v_norm->SetTitle("LHRS acceptance function (normalized);ev.Th (deg)");

     c3->cd(0);
     TPad *p4 = new TPad("p4","p4",0.8,0.0,1.,1.);
     p4->SetLeftMargin(0.);
     p4->Draw();
     p4->cd();
     TLatex tex1;
     tex1.SetTextAlign(11);
     tex1.SetTextSize(0.05);
     tex1.DrawLatexNDC(0.0,0.8,"thrown-in cut:");
     tex1.DrawLatexNDC(0.1,0.75,ispb);
     tex1.DrawLatexNDC(0.1,0.7,ecut);
     tex1.DrawLatexNDC(0.0,0.65,"accepted cut:");
     tex1.DrawLatexNDC(0.1,0.6,ispb);
     tex1.DrawLatexNDC(0.1,0.55,colcut);
     tex1.DrawLatexNDC(0.1,0.5,dpcut);
     tex1.DrawLatexNDC(0.1,0.45,xcut);
     tex1.DrawLatexNDC(0.1,0.4,"upcut");
     tex1.DrawLatexNDC(0.1,0.35,"downcut");


}
