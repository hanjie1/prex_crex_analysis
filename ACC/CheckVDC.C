#include "CollimatorL.C"

void CheckVDC(){
     TChain *T = new TChain("T");
     //T->Add("Rootfiles/PREXLHRSnom*");
     T->Add("/lustre19/expphy/volatile/halla/parity/ryanrich/Ebeam953_4by6/TargNominal/Zero1_SandwichLHRS_PREX*");

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

     TH1F *hvdc_x = new TH1F("hvdc_x","vdc x distribution",100,-0.8,0.8);
     T->Draw("x_vdc_tr>>hvdc_x",ACC*"rate/10.","hist");

     TGraph *geff = new TGraph();
     TGraph *gnew_x = new TGraph();
     Double_t eff[100]={0};
     for(int ii=1; ii<101; ii++){
	Double_t bin_rate = hvdc_x->GetBinContent(ii)/1000.;
   	Double_t bin_x = hvdc_x->GetBinCenter(ii);

        eff[ii-1] = 0.985157-0.000566364*bin_rate+3.87401e-06*pow(bin_rate,2)-1.38876e-08*pow(bin_rate,3)+1.16954e-11*pow(bin_rate,4);
        geff->SetPoint(ii-1,bin_x,eff[ii-1]);
        gnew_x->SetPoint(ii-1,bin_x,eff[ii-1]*bin_rate);
     }     

     cout<<"orignal x mean                   "<<hvdc_x->GetMean()<<endl;
     cout<<"after acceptance applied x mean  "<<gnew_x->GetMean()<<endl;

     TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
     c1->Divide(3,1);
     c1->cd(1);
     hvdc_x->Draw("HIST");
     hvdc_x->SetTitle("rate vs. x;x_vdc_tr;rate(Hz);");

     TLatex tex;
     tex.SetTextAlign(11);
     tex.SetTextSize(0.03);
     tex.DrawLatexNDC(0.6,0.6,Form("x mean=%.4f",hvdc_x->GetMean()));
    
     c1->cd(2);
     geff->SetMarkerStyle(8);
     geff->SetMarkerColor(4);
     geff->Draw("AP"); 
     geff->SetTitle("eff vs. x;x_vdc_tr;eff;");
     c1->cd(3);
     gnew_x->SetMarkerStyle(8);
     gnew_x->SetMarkerColor(4);
     gnew_x->Draw("AP"); 
     gnew_x->SetTitle("rate*eff vs. x;x_vdc_tr;rate*eff;");

     TLatex tex1;
     tex1.SetTextAlign(11);
     tex1.SetTextSize(0.03);
     tex1.DrawLatexNDC(0.6,0.6,Form("x mean=%.4f",gnew_x->GetMean()));
    
     T->Draw(">>electron",ACC);
     TEventList *electron;
     gDirectory->GetObject("electron",electron);
     T->SetEventList(electron);

     T->SetBranchStatus("*",0);
     T->SetBranchStatus("rate",1);
     T->SetBranchStatus("ev.Q2",1);
     T->SetBranchStatus("ev.Th",1);
     T->SetBranchStatus("ev.A",1);
     T->SetBranchStatus("x_vdc_tr",1);

     Double_t rate,aQ2, aTh, x_vdc_tr,aAsy;
     T->SetBranchAddress("rate",&rate);
     T->SetBranchAddress("ev.Q2",&aQ2);
     T->SetBranchAddress("ev.Th",&aTh);
     T->SetBranchAddress("ev.A",&aAsy);
     T->SetBranchAddress("x_vdc_tr",&x_vdc_tr);

     TH1F *hQ2 = new TH1F("hQ2","ev.Q2",100,0,0.015); 
     TH1F *hAsy = new TH1F("hAsy","ev.A",200,200,1000); 
     TH1F *hTH= new TH1F("hTH","ev.Th",100,3,8); 

     TH1F *hQ2_new = new TH1F("hQ2_new","ev.Q2 with efficiency applied",100,0,0.015); 
     TH1F *hTH_new= new TH1F("hTH_new","ev.Th with efficiency applied",100,3,8); 
     TH1F *hAsy_new = new TH1F("hAsy_new","ev.A with efficiency applied",200,200,1000); 
     Int_t nentries=electron->GetN();
     for(int ii=0;ii<nentries;ii++){
         T->GetEntry(electron->GetEntry(ii));

	 Double_t delta_x=0.8*2/100.0;
	 int nbin=0;
   	 for(int jj=0; jj<100; jj++){
	    Double_t xlow=-0.8+jj*delta_x;
 	    Double_t xhigh=xlow+delta_x;
	    if(x_vdc_tr<xhigh && x_vdc_tr>=xlow){ nbin=jj; break;} 
	 } 
         if(nbin>0){
	   hQ2->Fill(aQ2,rate); 
	   hTH->Fill(aTh,rate); 
	   hAsy->Fill(aAsy,rate); 

	   hQ2_new->Fill(aQ2,rate*eff[nbin]); 
	   hTH_new->Fill(aTh,rate*eff[nbin]); 
	   hAsy_new->Fill(aAsy,rate*eff[nbin]); 
	 }
     }

     TCanvas*c2 = new TCanvas("c2","c2",1200,1000);
     c2->Divide(3,1);
     c2->cd(1);
     hQ2->Draw();
     hQ2_new->SetLineColor(2);
     hQ2_new->Draw("same");
     TLatex tex2;
     tex2.SetTextAlign(11);
     tex2.SetTextSize(0.03);
     tex2.DrawLatexNDC(0.6,0.6,Form("#color[4]{origin mean=%.4f}",hQ2->GetMean()));
     tex2.DrawLatexNDC(0.6,0.5,Form("#color[2]{*eff mean=%.4f}",hQ2_new->GetMean()));

     c2->cd(2);
     hTH->Draw();
     hTH_new->SetLineColor(2);
     hTH_new->Draw("same");

     TLatex tex3;
     tex3.SetTextAlign(11);
     tex3.SetTextSize(0.03);
     tex3.DrawLatexNDC(0.6,0.6,Form("#color[4]{origin mean=%.4f}",hTH->GetMean()));
     tex3.DrawLatexNDC(0.6,0.5,Form("#color[2]{*eff mean=%.4f}",hTH_new->GetMean()));

     c2->cd(3);
     hAsy->Draw();
     hAsy_new->SetLineColor(2);
     hAsy_new->Draw("same");

     TLatex tex4;
     tex4.SetTextAlign(11);
     tex4.SetTextSize(0.03);
     tex4.DrawLatexNDC(0.6,0.6,Form("#color[4]{origin mean=%.4f}",hAsy->GetMean()));
     tex4.DrawLatexNDC(0.6,0.5,Form("#color[2]{*eff mean=%.4f}",hAsy_new->GetMean()));


     cout<<"Q2 original mean  "<<hQ2->GetMean()<<endl;
     cout<<"Q2 new mean       "<<hQ2_new->GetMean()<<endl;
     cout<<"TH original mean  "<<hTH->GetMean()<<endl;
     cout<<"TH new mean       "<<hTH_new->GetMean()<<endl;
     cout<<"Asy original mean  "<<hAsy->GetMean()<<endl;
     cout<<"Asy new mean       "<<hAsy_new->GetMean()<<endl;
}
