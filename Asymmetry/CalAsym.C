#include "CollimatorL.C"
#include "LoadAsym.h"

void CalAsym(){

     Double_t ppb = 1e9;

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

     //TCut COLCUT = colCut+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9+xcut;

     // lead cut
     TString ispb = "ev.nuclA==208";
     TCut isPb = Form("%s",ispb.Data());

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl;

     // dp cut
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

//     TCut ACC = DP+isPb+colCut+XCUT+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9;
     TCut ACC = DP+isPb+colCut+XCUT;

     int nbins = 200;
     TH1F *hA_v = new TH1F("hA_v","vertex asymmetry",nbins,200,1000);
     TH1F *hA_ap = new TH1F("hA_ap","apparent asymmetry",nbins,200,1000);
     TH1F *hA_e = new TH1F("hA_e","elemental asymmetry",nbins,200,1000);

     TH1F *hQ2_v = new TH1F("hQ2_v","vertex Q2",nbins,0,0.015);
     TH1F *hQ2_ap = new TH1F("hQ2_ap","apparent Q2",nbins,0,0.015);
     TH1F *hQ2_e = new TH1F("hQ2_e","elemental Q2",nbins,0,0.015);

     // vertex asymmetry
     T->Draw("ev.A>>hA_v",ACC*"rate");
     T->Draw("ev.Q2>>hQ2_v",ACC*"rate");

     // apparent asymmetry
     LoadTable("Tables/horpb.dat",0);   // not streched 
     LoadTable("Tables/horpb1.dat",1);   // streched 

     double beamE = 0.9534; //GeV  set in g4hrs
     double tpx, tpy, tpz, th_ztarg,rate,thisV,thisTh;
     T->SetBranchAddress("ev.tpx",&tpx);     
     T->SetBranchAddress("ev.tpy",&tpy);     
     T->SetBranchAddress("ev.tpz",&tpz);     
     T->SetBranchAddress("th_ztarg",&th_ztarg);     
     T->SetBranchAddress("rate",&rate);     
     T->SetBranchAddress("ev.V",&thisV);     
     T->SetBranchAddress("ev.Th",&thisTh);     
	
     T->Draw(">>goodMC",ACC);
     TEventList *goodMC;
     gDirectory->GetObject("goodMC",goodMC);
     T->SetEventList(goodMC);

     Double_t nMC = goodMC->GetN();
     for(int ii=0; ii<nMC; ii++){
	T->GetEntry(goodMC->GetEntry(ii));
	// apparent asymmetry 
 	double theta_ap = th_ztarg*180./TMath::Pi();
	double asym_ap = Interpolate(beamE,theta_ap,0,1);
	double asym_st_ap = Interpolate(beamE,theta_ap,1,1);  // stretched asymmetry

	double ep_ap = beamE/(1.+beamE/(208*0.931494028)*(1.-cos(th_ztarg)));
	double q2_ap = 4.*beamE*ep_ap*pow(sin(th_ztarg/2.),2);

	// elemental asymmetry 
	//double theta_e = atan(sqrt(tpx*tpx+tpy*tpy)/tpz)*180./TMath::Pi();
	double theta_e = thisTh;
	double asym_e = Interpolate(beamE,theta_e,0,1);
	double xs_e = Interpolate(beamE,theta_e,0,0);
	double asym_st_e = Interpolate(beamE,theta_e,1,1);  // stretched asymmetry

	double ep_e = beamE/(1.+beamE/(208*0.931494028)*(1.-cos(th_ztarg)));
	double q2_e = 4.*beamE*ep_e*pow(sin(theta_e*TMath::Pi()/180./2.),2);

	hA_ap->Fill(asym_ap*ppb,rate);
	hQ2_ap->Fill(q2_ap,rate);
	hA_e->Fill(asym_e*ppb,xs_e*thisV);
	hQ2_e->Fill(q2_e,xs_e*thisV);
     }

     Double_t inte1 = hA_v->Integral();
     Double_t inte2 = hA_e->Integral();
     hA_e->Scale(inte1/inte2);

     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     hA_e->SetLineWidth(2);
     hA_ap->SetLineWidth(2);
     hA_v->SetLineWidth(2);
     hA_e->SetLineColor(2);
     hA_ap->SetLineColor(1);
     hA_v->SetLineColor(4);

//     hA_ap->Draw();
     hA_e->Draw();
     hA_v->Draw("sames");
     hA_ap->SetTitle("asymmetry distribution");

     TLegend *leg1 = new TLegend(0.15,0.65,0.3,0.8);
     leg1->AddEntry(hA_e,"elemental");
     leg1->AddEntry(hA_v,"vertex");
     //leg1->AddEntry(hA_ap,"apparent");
     leg1->Draw();

   gPad->Update();
   TPaveStats* statA1 = (TPaveStats*)hA_e->FindObject("stats");
   TPaveStats* statA2 = (TPaveStats*)hA_v->FindObject("stats");
   statA1->SetY1NDC(0.90);
   statA1->SetY2NDC(0.75);
   statA2->SetY1NDC(0.75);
   statA2->SetY2NDC(0.60);
   statA1->SetTextColor(1);
   statA2->SetTextColor(6);
   gPad->Modified();

     inte1 = hQ2_v->Integral();
     inte2 = hQ2_e->Integral();
     hQ2_e->Scale(inte1/inte2);

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     hQ2_e->SetLineWidth(2);
     hQ2_ap->SetLineWidth(2);
     hQ2_v->SetLineWidth(2);
     hQ2_e->SetLineColor(2);
     hQ2_ap->SetLineColor(1);
     hQ2_v->SetLineColor(4);

   //  hQ2_ap->Draw();
     hQ2_e->Draw();
     hQ2_v->Draw("sames");
     hQ2_ap->SetTitle("Q2 distribution");

     TLegend *leg2 = new TLegend(0.6,0.65,0.75,0.8);
     leg2->AddEntry(hQ2_e,"elemental");
     leg2->AddEntry(hQ2_v,"vertex");
     //leg2->AddEntry(hQ2_ap,"apparent");
     leg2->Draw();

   gPad->Update();
   TPaveStats* statA3 = (TPaveStats*)hQ2_e->FindObject("stats");
   TPaveStats* statA4 = (TPaveStats*)hQ2_v->FindObject("stats");
   statA3->SetY1NDC(0.90);
   statA3->SetY2NDC(0.75);
   statA4->SetY1NDC(0.75);
   statA4->SetY2NDC(0.60);
   statA3->SetTextColor(1);
   statA4->SetTextColor(6);
   gPad->Modified();

     Double_t mean_v = hA_v->GetMean();
     Double_t mean_v_neff = hA_v->GetEffectiveEntries(); 
     Double_t mean_v_std= hA_v->GetMeanError(); 
     Double_t mean_v_err = mean_v_std * sqrt(mean_v_neff/(mean_v_neff-1.));

     Double_t mean_e = hA_e->GetMean();
     Double_t mean_e_neff = hA_e->GetEffectiveEntries(); 
     Double_t mean_e_std= hA_e->GetMeanError(); 
     Double_t mean_e_err = mean_e_std * sqrt(mean_e_neff/(mean_e_neff-1.));

     Double_t RC = mean_e/mean_v;
     Double_t RC_err = RC*sqrt(pow(mean_e_err/mean_e,2)+pow(mean_v_err/mean_v,2)); 

     cout<<"Asym:  "<<RC<<"  "<<RC_err<<endl;

     Double_t q2mean_v = hQ2_v->GetMean();
     Double_t q2mean_v_neff = hQ2_v->GetEffectiveEntries(); 
     Double_t q2mean_v_std= hQ2_v->GetMeanError(); 
     Double_t q2mean_v_err = q2mean_v_std * sqrt(q2mean_v_neff/(q2mean_v_neff-1.));

     Double_t q2mean_e = hQ2_e->GetMean();
     Double_t q2mean_e_neff = hQ2_e->GetEffectiveEntries(); 
     Double_t q2mean_e_std= hQ2_e->GetMeanError(); 
     Double_t q2mean_e_err = q2mean_e_std * sqrt(q2mean_e_neff/(q2mean_e_neff-1.));

     Double_t q2_ratio = q2mean_e/q2mean_v;
     Double_t q2_ratio_err = q2_ratio*sqrt(pow(q2mean_e_err/q2mean_e,2)+pow(q2mean_v_err/q2mean_v,2)); 

     cout<<"Q2:  "<<q2_ratio<<"  "<<q2_ratio_err<<endl;



}
