#include "CollimatorL.C"

void CheckMC(){
   
   bool islogy = false;

   TChain *T_MC = new TChain("T");

   TString target;
   TString tg;
   cout<<"Which target? (Pb, C12, Ca40) "<<endl;
   cin>>target;
   int run_number=0;
   if(target == "Pb") { run_number = 2055; tg = "pb208"; }
   if(target == "C12") { run_number = 2302; tg = "c12"; }
   if(target == "Ca40") { run_number = 2297; tg = "ca40"; }

   TString MCfile = Form("g4hrs_rootfiles/prex_%s.root",tg.Data()); 
   T_MC->Add(MCfile);

/*** MC cut ***/
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
 
 
   TCut collcut = colCut+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9+xcut;
   TCut phtg_tail_cut = "ph_tg_tr<-0.014 && ph_tg_tr>-0.02";
   TCut Q2_tail_cut = "ev.Q2<0.0045";

   TCut MCcut = collcut + Q2_tail_cut;

     gStyle->SetStatY(0.9);
     gStyle->SetStatX(0.31);
     gStyle->SetStatW(0.15);
     gStyle->SetStatH(0.2);

//======== focal plane ===========//

   Double_t trx_min=0.,trx_max=0.,try_min=0,try_max=0;
   Double_t trth_min=0.,trth_max=0.,trph_min=0,trph_max=0;

   if(islogy){
	trx_min=-1; trx_max=1; try_min=-0.3; try_max=0.3;
	trth_min=-0.5; trth_max=0.5; trph_min=-0.5; trph_max=0.5;
   }
   else{
	trx_min=-0.3; trx_max=0.1; try_min=-0.1; try_max=0.04;
	trth_min=-0.2; trth_max=0.05; trph_min=-0.08; trph_max=0.08;
   }

   TH1F *htrx_MC = new TH1F("htrx_MC","x_fp_tr distribution",100,trx_min,trx_max);
   TH1F *htry_MC = new TH1F("htry_MC","y_fp_tr distribution",100,try_min,try_max);
   TH1F *htrth_MC = new TH1F("htrth_MC","th_fp_tr distribution",100,trth_min,trth_max);
   TH1F *htrph_MC = new TH1F("htrph_MC","ph_fp_tr distribution",100,trph_min,trph_max);
   TH1F *htrx_MC_tail = new TH1F("htrx_MC_tail","x_fp_tr distribution",100,trx_min,trx_max);
   TH1F *htry_MC_tail = new TH1F("htry_MC_tail","y_fp_tr distribution",100,try_min,try_max);
   TH1F *htrth_MC_tail = new TH1F("htrth_MC_tail","th_fp_tr distribution",100,trth_min,trth_max);
   TH1F *htrph_MC_tail = new TH1F("htrph_MC_tail","ph_fp_tr distribution",100,trph_min,trph_max);

   TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
   c1->Divide(2,2);
   c1->cd(1);
   T_MC->Draw("x_fp_tr>>htrx_MC",collcut*"rate","");
   T_MC->Draw("x_fp_tr>>htrx_MC_tail",MCcut*"rate","");

   htrx_MC->SetLineColor(4);
   htrx_MC->SetLineWidth(2);
   htrx_MC_tail->SetLineColor(2);
   htrx_MC_tail->SetLineWidth(2);

   htrx_MC->Draw("HIST");
   htrx_MC_tail->Draw("HIST same");
   htrx_MC->SetTitle(Form("%s MC focal plane x; x;",tg.Data()));

   if(islogy) gPad->SetLogy();

   c1->cd(2);
   T_MC->Draw("y_fp_tr>>htry_MC",collcut*"rate","");
   T_MC->Draw("y_fp_tr>>htry_MC_tail",MCcut*"rate","");

   htry_MC->SetLineColor(4);
   htry_MC->SetLineWidth(2);
   htry_MC_tail->SetLineColor(2);
   htry_MC_tail->SetLineWidth(2);

   htry_MC->Draw("HIST");
   htry_MC_tail->Draw("HIST same");
   htry_MC->SetTitle(Form("%s MC focal plane y; y;",tg.Data()));
   if(islogy) gPad->SetLogy();

   c1->cd(3);
   T_MC->Draw("th_fp_tr>>htrth_MC",collcut*"rate","");
   T_MC->Draw("th_fp_tr>>htrth_MC_tail",MCcut*"rate","");

   htrth_MC->SetLineColor(4);
   htrth_MC->SetLineWidth(2);
   htrth_MC_tail->SetLineColor(2);
   htrth_MC_tail->SetLineWidth(2);

   htrth_MC->Draw("HIST");
   htrth_MC_tail->Draw("HIST same");
   htrth_MC->SetTitle(Form("%s MC focal plane theta; theta;",tg.Data()));
   if(islogy) gPad->SetLogy();

   c1->cd(4);
   T_MC->Draw("ph_fp_tr>>htrph_MC",collcut*"rate","");
   T_MC->Draw("ph_fp_tr>>htrph_MC_tail",MCcut*"rate","");

   htrph_MC->SetLineColor(4);
   htrph_MC->SetLineWidth(2);
   htrph_MC->SetLineColor(2);
   htrph_MC->SetLineWidth(2);

   htrph_MC->Draw("HIST");
   htrph_MC_tail->Draw("HIST same");
   htrph_MC->SetTitle(Form("%s MC focal plane phi; phi;",tg.Data()));
   if(islogy) gPad->SetLogy();


//======== target ===========//

   Double_t tgth_min=0.0,tgth_max=0,tgph_min=0,tgph_max=0;
   if(islogy){
	tgth_min = -0.5; tgth_max=0.5; tgph_min=-0.5; tgph_max=0.5;
   }
   else{
	tgth_min = -0.08; tgth_max=0.08; tgph_min=-0.04; tgph_max=0.04;
   }

   TH1F *htgth_MC = new TH1F("htgth_MC","th_tg_tr distgibution",100,tgth_min,tgth_max);
   TH1F *htgph_MC = new TH1F("htgph_MC","ph_tg_tr distgibution",100,tgph_min,tgph_max);
   TH1F *htgth_MC_tail = new TH1F("htgth_MC_tail","th_tg_tr distgibution",100,tgth_min,tgth_max);
   TH1F *htgph_MC_tail = new TH1F("htgph_MC_tail","ph_tg_tr distgibution",100,tgph_min,tgph_max);

   TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
   c2->Divide(2,2);
   c2->cd(1);
   T_MC->Draw("th_tg_tr>>htgth_MC",collcut*"rate","");
   T_MC->Draw("th_tg_tr>>htgth_MC_tail",MCcut*"rate","");

   htgth_MC->SetLineColor(4);
   htgth_MC->SetLineWidth(2);
   htgth_MC_tail->SetLineColor(2);
   htgth_MC_tail->SetLineWidth(2);

   htgth_MC->Draw("HIST");
   htgth_MC_tail->Draw("HIST same");
   htgth_MC->SetTitle(Form("%s MC target theta; theta;",tg.Data()));
   if(islogy) gPad->SetLogy();

   c2->cd(2);
   T_MC->Draw("ph_tg_tr>>htgph_MC",collcut*"rate","");
   T_MC->Draw("ph_tg_tr>>htgph_MC_tail",MCcut*"rate","");

   htgph_MC->SetLineColor(4);
   htgph_MC->SetLineWidth(2);
   htgph_MC_tail->SetLineColor(2);
   htgph_MC_tail->SetLineWidth(2);

   htgph_MC->Draw("HIST");
   htgph_MC_tail->Draw("HIST same");
   htgph_MC->SetTitle(Form("%s MC target phi; phi;",tg.Data()));
   if(islogy) gPad->SetLogy();

   Double_t tgdp_min=0, tgdp_max=0;
   if(islogy) { tgdp_min = -1; tgdp_max = 1; }
   else { tgdp_min = -0.008; tgdp_max = 0.; }

   TH1F *htgdp_MC = new TH1F("htgdp_MC","dp_fp_tr distgibution",100,tgdp_min,tgdp_max);
   TH1F *htgdp_MC_tail = new TH1F("htgdp_MC_tail","dp_dp_tr distgibution",100,tgdp_min,tgdp_max);

   c2->cd(3);
   T_MC->Draw("(p_fp_tr-950.0)/950.0>>htgdp_MC",collcut*"rate","");
   T_MC->Draw("(p_fp_tr-950.0)/950.0>>htgdp_MC_tail",MCcut*"rate","");
   
   htgdp_MC->SetLineColor(4);
   htgdp_MC->SetLineWidth(2);
   htgdp_MC_tail->SetLineColor(2);
   htgdp_MC_tail->SetLineWidth(2);

   htgdp_MC->Draw("HIST");
   htgdp_MC_tail->Draw("HIST same");
   htgdp_MC->SetTitle(Form("%s MC focal plane dp; dp;",tg.Data()));
   if(islogy) gPad->SetLogy();

   Double_t q2_min=0, q2_max=0;
   if(islogy) { q2_min = 0; q2_max = 0.1; }
   else { q2_min = 0; q2_max = 0.015; }

   TH1F *hq2_MC = new TH1F("hq2_MC","MC Q2 distgibution",100,q2_min,q2_max);
   TH1F *hq2_MC_tail = new TH1F("hq2_MC_tail","MC Q2 distgibution",100,q2_min,q2_max);

   TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
   T_MC->Draw("ev.Q2>>hq2_MC",collcut*"rate","");
   T_MC->Draw("ev.Q2>>hq2_MC_tail",MCcut*"rate","");
   
   hq2_MC->SetLineColor(4);
   hq2_MC->SetLineWidth(2);
   hq2_MC_tail->SetLineColor(2);
   hq2_MC_tail->SetLineWidth(2);

   hq2_MC->Draw("HIST");
   hq2_MC_tail->Draw("HIST same");
   hq2_MC->SetTitle(Form("%s MC ev.Q2; Q2;",tg.Data()));
   if(islogy) gPad->SetLogy();


   TString out1, out2, out3;
   if(islogy){
     out1 = Form("outfiles/%s_checkMC_log.pdf[",tg.Data());
     out2 = Form("outfiles/%s_checkMC_log.pdf",tg.Data());
     out3 = Form("outfiles/%s_checkMC_log.pdf]",tg.Data());
   }
   else{
     out1 = Form("outfiles/%s_checkMC.pdf[",tg.Data());
     out2 = Form("outfiles/%s_checkMC.pdf",tg.Data());
     out3 = Form("outfiles/%s_checkMC.pdf]",tg.Data());
   }



   c1->Print(out1.Data());
   c1->Print(out2.Data());
   c2->Print(out2.Data());
   c3->Print(out2.Data());
   c3->Print(out3.Data());

}
