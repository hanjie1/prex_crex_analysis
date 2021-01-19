#include "CollimatorR.C"

void CompareFP(){
   
   bool islogy = false;

   TChain *T_data = new TChain("T");
   TChain *T_MC = new TChain("T");

   TString target;
   TString tg;
   cout<<"Which target? (Pb, C12, Ca40) "<<endl;
   cin>>target;
   int run_number=0;
   if(target == "Pb") { run_number = 21188; tg = "pb208"; }
   if(target == "C12") { run_number = 21424; tg = "c12"; }
   if(target == "Ca40") { run_number = 21419; tg = "ca40"; }


   TString datafile = Form("/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexRHRS_%d_-1_file0.root",run_number);
   TString MCfile = Form("../g4hrs_rootfiles/prex_%s_RHRS.root",tg.Data()); 

   T_data->Add(datafile);
   T_MC->Add(MCfile);

/*** data cut ***/
   TCut TRK = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";

/*** MC cut ***/
double tol=0.0;
double xmin[11] = { 0.0, 0.0, 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[11] = {0.1, 0.1, 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[11] = {-0.05, -0.05, -0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[11] = {0.05, 0.05,0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

TCut xcut = "x_fp_tr!=-333.";//&&x_vdc_tr>0.0";
TCut colCut = "CollimatorR(x_col_tr,y_col_tr)";
TCut radCut1 = Form("-x_zup1 > %f && -x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
TCut radCut2 = Form("-x_zup2 > %f && -x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
TCut downCut1 = Form("-x_zdown1 > (%f+%f) && -x_zdown1 < (%f-%f) && y_zdown1 > (%f+1.5*%f) && y_zdown1 < (%f-1.5*%f)",xmin[2],tol,xmax[2],tol,ymin[2],tol,ymax[2],tol);
TCut downCut2 = Form("-x_zdown2 > (%f+%f) && -x_zdown2 < (%f-%f) && y_zdown2 > (%f+1.5*%f) && y_zdown2 < (%f-1.5*%f)",xmin[3],tol,xmax[3],tol,ymin[3],tol,ymax[3],tol);
TCut downCut3 = Form("-x_zdown3 > (%f+%f) && -x_zdown3 < (%f-%f) && y_zdown3 > (%f+1.5*%f) && y_zdown3 < (%f-1.5*%f)",xmin[4],tol,xmax[4],tol,ymin[4],tol,ymax[4],tol);
TCut downCut4 = Form("-x_zdown4 > (%f+%f) && -x_zdown4 < (%f-%f) && y_zdown4 > (%f+1.5*%f) && y_zdown4 < (%f-1.5*%f)",xmin[5],tol,xmax[5],tol,ymin[5],tol,ymax[5],tol);
TCut downCut5 = Form("-x_zdown5 > (%f+%f) && -x_zdown5 < (%f-%f) && y_zdown5 > (%f+1.5*%f) && y_zdown5 < (%f-1.5*%f)",xmin[6],tol,xmax[6],tol,ymin[6],tol,ymax[6],tol);
TCut downCut6 = Form("-x_zdown6 > (%f+%f) && -x_zdown6 < (%f-%f) && y_zdown6 > (%f+1.5*%f) && y_zdown6 < (%f-1.5*%f)",xmin[7],tol,xmax[7],tol,ymin[7],tol,ymax[7],tol);
TCut downCut7 = Form("-x_zdown7 > (%f+%f) && -x_zdown7 < (%f-%f) && y_zdown7 > (%f+1.5*%f) && y_zdown7 < (%f-1.5*%f)",xmin[8],tol,xmax[8],tol,ymin[8],tol,ymax[8],tol);
TCut downCut8 = Form("-x_zdown8 > (%f+%f) && -x_zdown8 < (%f-%f) && y_zdown8 > (%f+1.5*%f) && y_zdown8 < (%f-1.5*%f)",xmin[9],tol,xmax[9],tol,ymin[9],tol,ymax[9],tol);
TCut downCut9 = Form("-x_zdown9 > (%f+%f) && -x_zdown9 < (%f-%f) && y_zdown9 > (%f+1.5*%f) && y_zdown9 < (%f-1.5*%f)",xmin[10],tol,xmax[10],tol,ymin[10],tol,ymax[10],tol);
 
   TCut MCcut = colCut+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9+xcut;

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

   TH1F *htrx_vdc = new TH1F("htrx_vdc","R.tr.x distribution",100,trx_min,trx_max);
   TH1F *htrx_fp = new TH1F("htrx_fp","R.tr.x at fp distribution",100,trx_min,trx_max);
   TH1F *htrx_MC_vdc = new TH1F("htrx_MC_vdc","x_vdc_tr distribution",100,trx_min,trx_max);
   TH1F *htrx_MC_fp = new TH1F("htrx_MC_fp","x_fp_tr distribution",100,trx_min,trx_max);

   TH1F *htry_vdc = new TH1F("htry_vdc","R.tr.y distribution",100,try_min,try_max);
   TH1F *htry_fp = new TH1F("htry_fp","R.tr.y at fp distribution",100,try_min,try_max);
   TH1F *htry_MC_vdc = new TH1F("htry_MC_vdc","y_vdc_tr distribution",100,try_min,try_max);
   TH1F *htry_MC_fp = new TH1F("htry_MC_fp","y_fp_tr distribution",100,try_min,try_max);

   TH1F *htrth= new TH1F("htrth","R.tr.th distribution",100,trth_min,trth_max);
   TH1F *htrth_MC_vdc = new TH1F("htrth_MC_vdc","th_vdc_tr distribution",100,trth_min,trth_max);
   TH1F *htrth_MC_fp = new TH1F("htrth_MC_fp","th_fp_tr distribution",100,trth_min,trth_max);

   TH1F *htrph = new TH1F("htrph","R.tr.ph distribution",100,trph_min,trph_max);
   TH1F *htrph_MC_vdc = new TH1F("htrph_MC_vdc","ph_vdc_tr distribution",100,trph_min,trph_max);
   TH1F *htrph_MC_fp = new TH1F("htrph_MC_fp","ph_fp_tr distribution",100,trph_min,trph_max);

   TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
   c1->Divide(2,2);
   c1->cd(1);
   T_data->Draw("R.tr.x>>htrx_vdc",TRK);
   T_data->Draw("R.tr.x+0.9*R.tr.th>>htrx_fp",TRK);
   T_MC->Draw("x_vdc_tr>>htrx_MC_vdc",MCcut*"rate","same");
   T_MC->Draw("x_fp_tr>>htrx_MC_fp",MCcut*"rate","same");

   Double_t ntrx_vdc = htrx_vdc->Integral();
   Double_t ntrx_fp = htrx_fp->Integral();
   Double_t ntrx_mc_vdc = htrx_MC_vdc->Integral();
   Double_t ntrx_mc_fp = htrx_MC_fp->Integral();
   htrx_MC_vdc->Scale(ntrx_vdc/ntrx_mc_vdc);
   htrx_MC_fp->Scale(ntrx_vdc/ntrx_mc_fp);
   htrx_fp->Scale(ntrx_vdc/ntrx_fp);

   Double_t h1_max = htrx_vdc->GetBinContent(htrx_vdc->GetMaximumBin());
   Double_t h2_max = htrx_MC_vdc->GetBinContent(htrx_MC_vdc->GetMaximumBin());
  
   if(h1_max>h2_max){ 
      htrx_vdc->Draw("HIST");
      htrx_fp->Draw("HIST same");
      htrx_MC_vdc->Draw("HIST same");
      htrx_MC_fp->Draw("HIST same");
      htrx_vdc->SetTitle(Form("%s focal plane X; x;",target.Data()));
   }
   else{
      htrx_MC_vdc->Draw("HIST");
      htrx_MC_fp->Draw("HIST same");
      htrx_vdc->Draw("HIST same");
      htrx_fp->Draw("HIST same");
      htrx_MC_vdc->SetTitle(Form("RHRS %s focal plane X; x;",target.Data()));
   }
   htrx_vdc->SetLineColor(2);
   htrx_vdc->SetLineWidth(2);
   htrx_fp->SetLineColor(8);
   htrx_fp->SetLineWidth(2);
   htrx_MC_vdc->SetLineColor(4);
   htrx_MC_vdc->SetLineWidth(2);
   htrx_MC_fp->SetLineColor(6);
   htrx_MC_fp->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   TLegend *leg1 = new TLegend(0.15,0.6,0.35,0.8);
   leg1->AddEntry(htrx_vdc,"data vdc","L");
   leg1->AddEntry(htrx_fp,"data fp","L");
   leg1->AddEntry(htrx_MC_vdc,"MC vdc","L");
   leg1->AddEntry(htrx_MC_fp,"MC fp","L");
   leg1->Draw();

   c1->cd(2);
   T_data->Draw("R.tr.y>>htry_vdc",TRK);
   T_data->Draw("R.tr.y+0.9*R.tr.ph>>htry_fp",TRK);
   T_MC->Draw("y_vdc_tr>>htry_MC_vdc",MCcut*"rate","same");
   T_MC->Draw("y_fp_tr>>htry_MC_fp",MCcut*"rate","same");

   Double_t ntry_vdc = htry_vdc->Integral();
   Double_t ntry_fp = htry_fp->Integral();
   Double_t ntry_mc_vdc = htry_MC_vdc->Integral();
   Double_t ntry_mc_fp = htry_MC_fp->Integral();
   htry_MC_vdc->Scale(ntry_vdc/ntry_mc_vdc);
   htry_MC_fp->Scale(ntry_vdc/ntry_mc_fp);
   htry_fp->Scale(ntry_vdc/ntry_fp);

   h1_max = htry_vdc->GetBinContent(htry_vdc->GetMaximumBin());
   h2_max = htry_MC_vdc->GetBinContent(htry_MC_vdc->GetMaximumBin());
  
   if(h1_max>h2_max){ 
      htry_vdc->Draw("HIST");
      htry_fp->Draw("HIST same");
      htry_MC_vdc->Draw("HIST same");
      htry_MC_fp->Draw("HIST same");
      htry_vdc->SetTitle("focal plane Y; y;");
   }
   else{
      htry_MC_fp->Draw("HIST");
      htry_MC_vdc->Draw("HIST same");
      htry_vdc->Draw("HIST same");
      htry_fp->Draw("HIST same");
      htry_MC_fp->SetTitle("focal plane Y; y;");
   }
   htry_vdc->SetLineColor(2);
   htry_vdc->SetLineWidth(2);
   htry_fp->SetLineColor(8);
   htry_fp->SetLineWidth(2);
   htry_MC_vdc->SetLineColor(4);
   htry_MC_vdc->SetLineWidth(2);
   htry_MC_fp->SetLineColor(6);
   htry_MC_fp->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   c1->cd(3);
   T_data->Draw("R.tr.th>>htrth",TRK);
   T_MC->Draw("th_vdc_tr>>htrth_MC_vdc",MCcut*"rate","same");
   T_MC->Draw("th_fp_tr>>htrth_MC_fp",MCcut*"rate","same");

   Double_t ntrth = htrth->Integral();
   Double_t ntrth_mc_vdc = htrth_MC_vdc->Integral();
   Double_t ntrth_mc_fp = htrth_MC_fp->Integral();
   htrth_MC_vdc->Scale(ntrth/ntrth_mc_vdc);
   htrth_MC_fp->Scale(ntrth/ntrth_mc_fp);
   
   h1_max = htrth->GetBinContent(htrth->GetMaximumBin());
   h2_max = htrth_MC_vdc->GetBinContent(htrth_MC_vdc->GetMaximumBin());
  
   if(h1_max>h2_max){ 
      htrth->Draw("HIST");
      htrth_MC_vdc->Draw("HIST same");
      htrth_MC_fp->Draw("HIST same");
      htrth->SetTitle("focal plane Theta; theta;");
   }
   else{
      htrth_MC_vdc->Draw("HIST");
      htrth_MC_fp->Draw("HIST same");
      htrth->Draw("HIST same");
      htrth_MC_vdc->SetTitle("focal plane Theta; theta;");
   }
   htrth->SetLineColor(2);
   htrth->SetLineWidth(2);
   htrth_MC_vdc->SetLineColor(4);
   htrth_MC_vdc->SetLineWidth(2);
   htrth_MC_fp->SetLineColor(6);
   htrth_MC_fp->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   c1->cd(4);
   T_data->Draw("R.tr.ph>>htrph",TRK);
   T_MC->Draw("ph_vdc_tr>>htrph_MC_vdc",MCcut*"rate","same");
   T_MC->Draw("ph_fp_tr>>htrph_MC_fp",MCcut*"rate","same");

   Double_t ntrph = htrph->Integral();
   Double_t ntrph_mc_vdc = htrph_MC_vdc->Integral();
   Double_t ntrph_mc_fp = htrph_MC_fp->Integral();
   htrph_MC_vdc->Scale(ntrph/ntrph_mc_vdc);
   htrph_MC_fp->Scale(ntrph/ntrph_mc_fp);
   
   h1_max = htrph->GetBinContent(htrph->GetMaximumBin());
   h2_max = htrph_MC_vdc->GetBinContent(htrph_MC_vdc->GetMaximumBin());
  
   if(h1_max>h2_max){ 
      htrph->Draw("HIST");
      htrph_MC_vdc->Draw("HIST same");
      htrph_MC_fp->Draw("HIST same");
      htrph->SetTitle("focal plane Phi; phi;");
   }
   else{
      htrph_MC_vdc->Draw("HIST");
      htrph_MC_fp->Draw("HIST same");
      htrph->Draw("HIST same");
      htrph_MC_vdc->SetTitle("focal plane Phi; phi;");
   }
   htrph->SetLineColor(2);
   htrph->SetLineWidth(2);
   htrph_MC_vdc->SetLineColor(4);
   htrph_MC_vdc->SetLineWidth(2);
   htrph_MC_fp->SetLineColor(6);
   htrph_MC_fp->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   c1->Print(Form("outfiles/%s_FP_RHRS.pdf",tg.Data()));
}
