#include "CollimatorL.C"

void Accfunc(){
     TChain *T = new TChain("T");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("Rootfiles/Zero1_SandwichLHRS_PREX_0.0_5.root");

     /**  cuts **/
     double tol=0.0;
     double xmin[11] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
     double xmax[11] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
     double ymin[11] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
     double ymax[11] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };
 
     TString xcut = "x_fp_tr!=-333.";
     TString colcut = "CollimatorL(x_col_tr,y_col_tr)";
    
     TCut XCUT = Form("%s",xcut.Data());//&&x_vdc_tr>0.0";
     TCut colCut = Form("%s",colcut.Data());
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


     // lead cut
     TString ispb = "ev.nuclA==208";
     TCut isPb = Form("%s",ispb.Data());

     //TCanvas *c0 = new TCanvas("c0","c0",1000,1000);
     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl; 

     // dp cut
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

     // dp cut for the incident events
     double delta_p = p_peak*0.02;
     TString dpcut_inc = Form("(%f-p_ztarg)<%f",p_peak,delta_p);
     TCut DP_INC = Form("%s",dpcut_inc.Data());

     // beam E cut
     TString ecut = Form("ev.beamp>%f",p_peak/1000.);
     TCut ECUT = Form("%s",ecut.Data());

     TCut ACC = isPb+colCut+XCUT+DP;
     TCut INC = isPb+ECUT;

     int nbin = 100;

    // TCanvas *c00 = new TCanvas("c00","c00",1000,1000);
     TH1F *hinc = new TH1F("hinc","incident events angle distribution",nbin,3,8); 
     TH1F *hacc = new TH1F("hacc","accepted events angle distribution",nbin,3,8); 

     hinc->Sumw2();
     hacc->Sumw2();

     T->Draw("ev.th*180./TMath::Pi()>>hinc",INC*"rate");
     T->Draw("ev.th*180./TMath::Pi()>>hacc",ACC*"rate");

     TGraphErrors *gACC = new TGraphErrors();

     int nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc->GetBinContent(ii);
	if(ninc==0) continue;

	Double_t th_bin = hinc->GetBinCenter(ii);
	Double_t nacc = hacc->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;

        Double_t inc_sumw2_sqrt = hinc->GetBinError(ii); // sqrt(sumw2)
	Double_t tmp_err2 = tmp_ratio * (1.-tmp_ratio) * pow(inc_sumw2_sqrt,2)/pow(ninc,2);
	Double_t tmp_err = sqrt(tmp_err2);
		
	gACC->SetPoint(nn,th_bin,tmp_ratio);	
	gACC->SetPointError(nn,0,tmp_err);	
 	nn++;
     } 

     TH1F *hinc_v = new TH1F("hinc_v","incident events angle distribution",nbin,3,8); 
     TH1F *hacc_v = new TH1F("hacc_v","accepted events angle distribution",nbin,3,8); 

     hinc_v->Sumw2();
     hacc_v->Sumw2();

     //T->Draw("atan(sqrt(ev.tpx*ev.tpx+ev.tpy*ev.tpy)/ev.tpz)*180./TMath::Pi()>>hinc_v",INC*"rate");
     //T->Draw("atan(sqrt(ev.tpx*ev.tpx+ev.tpy*ev.tpy)/ev.tpz)*180./TMath::Pi()>>hacc_v",ACC*"rate");

     T->Draw("ev.Th>>hinc_v",INC*"rate");
     T->Draw("ev.Th>>hacc_v",ACC*"rate");

     TGraphErrors *gACC_v = new TGraphErrors();

     ofstream outfile;
     outfile.open("accfunction.csv");
     outfile<<"vertex angle,acceptance,stat_err"<<endl;

     nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc_v->GetBinContent(ii);
	Double_t th_bin = hinc_v->GetBinCenter(ii);
	if(ninc==0) {
	   outfile<<th_bin<<","<<0<<","<<0<<endl;
	   continue;
	}

	Double_t nacc = hacc_v->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;

        Double_t inc_sumw2_sqrt = hinc_v->GetBinError(ii); // sqrt(sumw2)
	Double_t tmp_err2 = tmp_ratio * (1.-tmp_ratio) * pow(inc_sumw2_sqrt,2)/pow(ninc,2);
	Double_t tmp_err = sqrt(tmp_err2);

	gACC_v->SetPoint(nn,th_bin,tmp_ratio);	
	gACC_v->SetPointError(nn,0,tmp_err);	

	outfile<<th_bin<<","<<tmp_ratio<<","<<tmp_err<<endl;
 	nn++;
     } 
     outfile.close();

     TH1F *hinc_v_norate = new TH1F("hinc_v_norate","incident events angle distribution without rate applied",nbin,3,8); 
     TH1F *hacc_v_norate = new TH1F("hacc_v_norate","accepted events angle distribution without rate applied",nbin,3,8); 

     hinc_v_norate->Sumw2();
     hacc_v_norate->Sumw2();

     //T->Draw("atan(sqrt(ev.tpx*ev.tpx+ev.tpy*ev.tpy)/ev.tpz)*180./TMath::Pi()>>hinc_v_norate",INC);
     //T->Draw("atan(sqrt(ev.tpx*ev.tpx+ev.tpy*ev.tpy)/ev.tpz)*180./TMath::Pi()>>hacc_v_norate",ACC);

     T->Draw("ev.Th>>hinc_v_norate",INC);
     T->Draw("ev.Th>>hacc_v_norate",ACC);

     TGraphErrors *gACC_v_norate = new TGraphErrors();

     nn = 0;
     for(int ii=1; ii<nbin+1; ii++){
	Double_t ninc = hinc_v_norate->GetBinContent(ii);
	if(ninc==0) continue;

	Double_t th_bin = hinc_v_norate->GetBinCenter(ii);
	Double_t nacc = hacc_v_norate->GetBinContent(ii);
	Double_t tmp_ratio = nacc/ninc;

        Double_t inc_sumw2_sqrt = hinc_v_norate->GetBinError(ii); // sqrt(sumw2)
	Double_t tmp_err2 = tmp_ratio * (1.-tmp_ratio) * pow(inc_sumw2_sqrt,2)/pow(ninc,2);
	Double_t tmp_err = sqrt(tmp_err2);

	gACC_v_norate->SetPoint(nn,th_bin,tmp_ratio);	
	gACC_v_norate->SetPointError(nn,0,tmp_err);	

 	nn++;
     } 


     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     c1->Draw();
     TPad *p1 = new TPad("p1","p1",0.0,0.0,0.8,1.);
     p1->SetRightMargin(0.1);
     p1->Draw();
     p1->cd();
     TMultiGraph *mg = new TMultiGraph();
     mg->Add(gACC);
     mg->Add(gACC_v);

     gACC->SetMarkerStyle(8);
     gACC->SetMarkerColor(4);
     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(2);
     mg->Draw("AP");
     mg->SetTitle("acceptance function;angle;");

     TLegend *leg = new TLegend(0.15,0.65,0.4,0.8);
     leg->AddEntry(gACC,"ev.th","P");
     //leg->AddEntry(gACC_v,"atan(#frac{#sqrt{ev.tpx^{2}+ev.tpy^{2}}}{ev.tpz})","P");
     leg->AddEntry(gACC_v,"ev.Th","P");
     leg->Draw();

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
     //tex.DrawLatexNDC(0.1,0.7,dpcut_inc);
     tex.DrawLatexNDC(0.1,0.7,ecut);
     tex.DrawLatexNDC(0.0,0.65,"accepted cut:");
     tex.DrawLatexNDC(0.1,0.6,ispb);
     tex.DrawLatexNDC(0.1,0.55,colcut);
     tex.DrawLatexNDC(0.1,0.5,dpcut);
     tex.DrawLatexNDC(0.1,0.45,xcut);

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hinc_v->Draw();
     hinc->Draw("same");
     hinc_v->SetLineColor(2);

     c2->cd(2);
     hacc->Draw();
     hacc_v->Draw("same");
     hacc_v->SetLineColor(2);

     TLegend *leg1 = new TLegend(0.15,0.65,0.4,0.8);
     leg1->AddEntry(hacc,"ev.th","L");
     //leg1->AddEntry(hacc_v,"atan(#frac{#sqrt{ev.tpx^{2}+ev.tpy^{2}}}{ev.tpz})","L");
     leg1->AddEntry(hacc_v,"ev.Th","L");
     leg1->Draw();

     TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
     TMultiGraph *mg1 = new TMultiGraph();
     mg1->Add(gACC_v);
     mg1->Add(gACC_v_norate);

     gACC_v_norate->SetMarkerStyle(8);
     gACC_v_norate->SetMarkerColor(4);
     gACC_v->SetMarkerStyle(8);
     gACC_v->SetMarkerColor(2);
     mg1->Draw("AP");
     mg1->SetTitle("acceptance function;;");

     TLegend *leg2 = new TLegend(0.15,0.65,0.4,0.8);
     leg2->AddEntry(gACC_v,"with rate","P");
     leg2->AddEntry(gACC_v_norate,"without rate","P");
     leg2->Draw();
/*
     TCanvas *c4 = new TCanvas("c4","c4",1500,1500);
     c4->Divide(2,1);
     c4->cd(1);
     //TH2F *hth_rate = new TH2F("hth_rate","rate vs. scattering angle",100,2,8,100,0,9e6); 
     T->Draw("rate:atan(sqrt(ev.tpx*ev.tpx+ev.tpy*ev.tpy)/ev.tpz)*180./TMath::Pi()",ACC,"COLZ");
   
     c4->cd(2);
     //TH1F *hth_ck = new TH1F("hth_ck","check scattering angle");
     T->Draw("(1.-208*0.931494028*(1./ev.ep-1./ev.beamp))-(ev.tpz/sqrt(pow(ev.tpx,2)+pow(ev.tpy,2)+pow(ev.tpz,2)))",isPb+"ev.beamp>0.6");
*/    


}
