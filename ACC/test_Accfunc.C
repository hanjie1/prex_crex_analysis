#include "CollimatorL.C"

void test_Accfunc(){
     TChain *T = new TChain("T");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Sandwich/DPbD_1_xs_A.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Sandwich/DPbD_2_xs_A.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Sandwich/DPbD_3_xs_A.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Sandwich/DPbD_4_xs_A.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Sandwich/DPbD_5_xs_A.root");

     /**  cuts **/
     double tol=0.0;
     double xmin[11] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
     double xmax[11] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
     double ymin[11] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
     double ymax[11] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };
     
     TCut xcut = "x_fp_tr!=-333.";//&&x_vdc_tr>0.0";
     TCut colCut = "CollimatorL(x_col_tr,y_col_tr)";
     TCut radCut1 = Form("x_zup1 > %f && x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
     TCut radCut2 = Form("x_zup2 > %f && x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
     TCut downCut1 = Form("x_zdown1 > (%f+%f) && x_zdown1 < (%f-%f) && y_zdown1 > (%f+1.5*%f) && y_zdown1 < (%f-1.5*%f)",xmin[2],tol,xmax[2],tol,ymin[2],tol,ymax[2],tol);
     TCut downCut2 = Form("x_zdown2 > (%f+%f) && x_zdown2 < (%f-%f) && y_zdown2 > (%f+1.5*%f) && y_zdown2 < (%f-1.5*%f)",xmin[3],tol,xmax[3],tol,ymin[3],tol,ymax[3],tol);
     TCut downCut3 = Form("x_zdown3 > (%f+%f) && x_zdown3 < (%f-%f) && y_zdown3 > (%f+1.5*%f) && y_zdown3 < (%f-1.5*%f)",xmin[4],tol,xmax[4],tol,ymin[4],tol,ymax[4],tol);
     TCut downCut4 = Form("x_zdown4 > (%f+%f) && x_zdown4 < (%f-%f) && y_zdown4 > (%f+1.5*%f) && y_zdown4 < (%f-1.5*%f)",xmin[5],tol,xmax[5],tol,ymin[5],tol,ymax[5],tol);
     TCut downCut5 = Form("x_zdown5 > (%f+%f) && x_zdown5 < (%f-%f) && y_zdown5 > (%f+1.5*%f) && y_zdown5 < (%f-1.5*%f)",xmin[6],tol,xmax[6],tol,ymin[6],tol,ymax[6],tol);
     TCut downCut6 = Form("x_zdown6 > (%f+%f) && x_zdown6 < (%f-%f) && y_zdown6 > (%f+1.5*%f) && y_zdown6 < (%f-1.5*%f)",xmin[7],tol,xmax[7],tol,ymin[7],tol,ymax[7],tol);
     TCut downCut7 = Form("x_zdown7 > (%f+%f) && x_zdown7 < (%f-%f) && y_zdown7 > (%f+1.5*%f) && y_zdown7 < (%f-1.5*%f)",xmin[8],tol,xmax[8],tol,ymin[8],tol,ymax[8],tol);
     TCut downCut8 = Form("x_zdown8 > (%f+%f) && x_zdown8 < (%f-%f) && y_zdown8 > (%f+1.5*%f) && y_zdown8 < (%f-1.5*%f)",xmin[9],tol,xmax[9],tol,ymin[9],tol,ymax[9],tol);
     TCut downCut9 = Form("x_zdown9 > (%f+%f) && x_zdown9 < (%f-%f) && y_zdown9 > (%f+1.5*%f) && y_zdown9 < (%f-1.5*%f)",xmin[10],tol,xmax[10],tol,ymin[10],tol,ymax[10],tol);

     TCut COLCUT = colCut+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9+xcut;

     // lead cut
     TCut isPb = "ev.nuclA==208";

     int nbin = 100;

     TH1F *hinc = new TH1F("hinc","incident events angle distribution",nbin,2,8); 
     TH1F *hacc = new TH1F("hacc","accepted events angle distribution",nbin,2,8); 

     T->Draw("th_ztarg*180./TMath::Pi()>>hinc",isPb);
     T->Draw("th_ztarg*180./TMath::Pi()>>hacc",COLCUT+isPb);

     TGraphErrors *gACC = new TGraphErrors();

     int nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc->GetBinContent(ii);
	if(ninc==0) continue;

	Double_t th_bin = hinc->GetBinCenter(ii);
	Double_t nacc = hacc->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;
	Double_t tmp_err = tmp_ratio * sqrt(1./nacc - 1./ninc);
		
	gACC->SetPoint(nn,th_bin,tmp_ratio);	
	gACC->SetPointError(nn,0,0);	
 	nn++;
     } 

     TH1F *hinc_v = new TH1F("hinc_v","incident events angle distribution*rate",nbin,2,8); 
     TH1F *hacc_v = new TH1F("hacc_v","accepted events angle distribution*rate",nbin,2,8); 

     T->Draw("th_ztarg*180./TMath::Pi()>>hinc_v",isPb*"rate");
     T->Draw("th_ztarg*180./TMath::Pi()>>hacc_v",(COLCUT+isPb)*"rate");


     TGraphErrors *gACC_v = new TGraphErrors();

     nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc_v->GetBinContent(ii);
	if(ninc==0) continue;

	Double_t th_bin = hinc_v->GetBinCenter(ii);
	Double_t nacc = hacc_v->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;
	Double_t tmp_err = tmp_ratio * sqrt(1./nacc - 1./ninc);
	
	if(tmp_ratio==0) tmp_err=0;	

	gACC_v->SetPoint(nn,th_bin,tmp_ratio);	
	gACC_v->SetPointError(nn,0,0);	

 	nn++;
     } 

     TCanvas *c1 = new TCanvas("c1","c1");
     TMultiGraph *mg = new TMultiGraph();
     mg->Add(gACC);
     mg->Add(gACC_v);

     gACC->SetMarkerStyle(8);
     gACC->SetMarkerColor(4);
     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(2);
     mg->Draw("AP");
     mg->SetTitle("acceptance function;;");

     TLegend *leg = new TLegend(0.15,0.65,0.4,0.8);
     leg->AddEntry(gACC,"without rate","P");
     leg->AddEntry(gACC_v,"with rate","P");
     leg->Draw();

     TLatex tex;
     tex.SetTextAlign(12);
     tex.SetTextSize(0.035);
     tex.DrawLatexNDC(0.65,0.8,"thrown-in cut: isPb");
     tex.DrawLatexNDC(0.65,0.7,"accepted cut: isPb+COLCUT");


}
