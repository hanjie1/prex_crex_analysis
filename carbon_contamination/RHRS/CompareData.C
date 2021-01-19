#include "CollimatorR.C"

void CompareData(){
   
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


   TString datafile = Form("/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexRHRS_%d_-1.root",run_number);
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

   TH1F *htrx = new TH1F("htrx","R.tr.x distribution",100,trx_min,trx_max);
   TH1F *htrx_MC = new TH1F("htrx_MC","x_fp_tr distribution",100,trx_min,trx_max);
   TH1F *htry = new TH1F("htry","R.tr.y distribution",100,try_min,try_max);
   TH1F *htry_MC = new TH1F("htry_MC","y_fp_tr distribution",100,try_min,try_max);
   TH1F *htrth = new TH1F("htrth","R.tr.th distribution",100,trth_min,trth_max);
   TH1F *htrth_MC = new TH1F("htrth_MC","th_fp_tr distribution",100,trth_min,trth_max);
   TH1F *htrph = new TH1F("htrph","R.tr.ph distribution",100,trph_min,trph_max);
   TH1F *htrph_MC = new TH1F("htrph_MC","ph_fp_tr distribution",100,trph_min,trph_max);

   TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
   c1->Divide(2,2);
   c1->cd(1);
   T_data->Draw("(R.tr.x+0.9*R.tr.th)>>htrx",TRK);
   T_MC->Draw("x_fp_tr>>htrx_MC",MCcut*"rate","same");

   Double_t ntrx = htrx->Integral();
   Double_t ntrx_mc = htrx_MC->Integral();
   htrx_MC->Scale(ntrx/ntrx_mc);

   Double_t h1_max = htrx->GetBinContent(htrx->GetMaximumBin());
   Double_t h2_max = htrx_MC->GetBinContent(htrx_MC->GetMaximumBin());

   if(h1_max>h2_max){
      htrx->Draw("HIST");
      htrx_MC->Draw("HIST same");
      htrx->SetTitle(Form("%s focal plane X; x;",target.Data()));
   }
   else{
      htrx_MC->Draw("HIST");
      htrx->Draw("HIST same");
      htrx_MC->SetTitle(Form("%s focal plane X; x;",target.Data()));
   }
   htrx->SetLineColor(2);
   htrx->SetLineWidth(2);
   htrx_MC->SetLineColor(4);
   htrx_MC->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   TLegend *leg1 = new TLegend(0.15,0.7,0.3,0.8);
   leg1->AddEntry(htrx,"data","L");
   leg1->AddEntry(htrx_MC,"MC","L");
   leg1->Draw();


   c1->cd(2);
   T_data->Draw("(R.tr.y+0.9*R.tr.ph)>>htry",TRK);
   T_MC->Draw("y_fp_tr>>htry_MC",MCcut*"rate","same");

   Double_t ntry = htry->Integral();
   Double_t ntry_mc = htry_MC->Integral();
   htry_MC->Scale(ntry/ntry_mc);

   h1_max = htry->GetBinContent(htry->GetMaximumBin());
   h2_max = htry_MC->GetBinContent(htry_MC->GetMaximumBin());

   if(h1_max>h2_max){
      htry->Draw("HIST");
      htry_MC->Draw("HIST same");
      htry->SetTitle("focal plane Y; y;");
   }
   else{
      htry_MC->Draw("HIST");
      htry->Draw("HIST same");
      htry_MC->SetTitle("focal plane Y; y;");
   }
   htry->SetLineColor(2);
   htry->SetLineWidth(2);
   htry_MC->SetLineColor(4);
   htry_MC->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

   c1->cd(3);
   T_data->Draw("R.tr.th>>htrth",TRK);
   T_MC->Draw("th_fp_tr>>htrth_MC",MCcut*"rate","same");

   Double_t ntrth = htrth->Integral();
   Double_t ntrth_mc = htrth_MC->Integral();
   htrth_MC->Scale(ntrth/ntrth_mc);

   h1_max = htrth->GetBinContent(htrth->GetMaximumBin());
   h2_max = htrth_MC->GetBinContent(htrth_MC->GetMaximumBin());

   if(h1_max>h2_max){
      htrth->Draw("HIST");
      htrth_MC->Draw("HIST same");
      htrth->SetTitle("focal plane Theta; theta;");
   }
   else{
      htrth_MC->Draw("HIST");
      htrth->Draw("HIST same");
      htrth_MC->SetTitle("focal plane Theta; theta;");
   }
   htrth->SetLineColor(2);
   htrth->SetLineWidth(2);
   htrth_MC->SetLineColor(4);
   htrth_MC->SetLineWidth(2);
   if(islogy) gPad->SetLogy();


   c1->cd(4);
   T_data->Draw("R.tr.ph>>htrph",TRK);
   T_MC->Draw("ph_fp_tr>>htrph_MC",MCcut*"rate","same");

   Double_t ntrph = htrph->Integral();
   Double_t ntrph_mc = htrph_MC->Integral();
   htrph_MC->Scale(ntrph/ntrph_mc);

   h1_max = htrph->GetBinContent(htrph->GetMaximumBin());
   h2_max = htrph_MC->GetBinContent(htrph_MC->GetMaximumBin());

   if(h1_max>h2_max){
      htrph->Draw("HIST");
      htrph_MC->Draw("HIST same");
      htrph->SetTitle("focal plane Phi; phi;");
   }
   else{
      htrph_MC->Draw("HIST");
      htrph->Draw("HIST same");
      htrph_MC->SetTitle("focal plane Phi; phi;");
   }
   htrph->SetLineColor(2);
   htrph->SetLineWidth(2);
   htrph_MC->SetLineColor(4);
   htrph_MC->SetLineWidth(2);
   if(islogy) gPad->SetLogy();

//======== target ===========//
   gStyle->SetOptStat(0);

   Double_t tgth_min=0.0,tgth_max=0,tgph_min=0,tgph_max=0;
   if(islogy){
	tgth_min = -0.5; tgth_max=0.5; tgph_min=-0.5; tgph_max=0.5;
   }
   else{
	tgth_min = -0.08; tgth_max=0.08; tgph_min=-0.04; tgph_max=0.04;
   }

   const int nplane = 4;

   TH1F *htgth = new TH1F("htgth","R.gold.th distgibution",100,tgth_min,tgth_max);
   TH1F *htgth_MC[nplane];
   TH1F *htgph = new TH1F("htgph","R.gold.ph distgibution",100,tgph_min,tgph_max);
   TH1F *htgph_MC[nplane];

   TString plane_name[nplane]= {"tg","zup1","zup2","zsieve"};
   int color[nplane]={1,4,6,8};

   for(int ii=0; ii<nplane; ii++){
      htgth_MC[ii] = new TH1F(Form("htgth_MC%d",ii+1),"th MC distgibution",100,tgth_min,tgth_max);
      htgph_MC[ii] = new TH1F(Form("htgph_MC%d",ii+1),"ph MC distgibution",100,tgph_min,tgph_max);
   }

   TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
   c2->Divide(2,2);
   c2->cd(1);
   T_data->Draw("R.gold.th>>htgth",TRK);
   htgth->SetLineColor(2);
   htgth->SetLineWidth(2);
   Double_t ntgth = htgth->Integral();
   for(int ii=0; ii<nplane; ii++){
      T_MC->Draw(Form("th_%s_tr>>htgth_MC%d",plane_name[ii].Data(),ii+1),MCcut*"rate");
      Double_t ntgth_mc = htgth_MC[ii]->Integral();
      htgth_MC[ii]->Scale(ntgth/ntgth_mc);
      htgth_MC[ii]->SetLineColor(color[ii]); 
      htgth_MC[ii]->SetLineWidth(2); 
   }

   htgth_MC[1]->Draw("HIST");
   htgth->Draw("HIST same");
   htgth_MC[0]->Draw("HIST same");
   for(int ii=2; ii<nplane; ii++)
     htgth_MC[ii]->Draw("HIST same");

   htgth_MC[1]->SetTitle(Form("%s target theta; th;",tg.Data()));
   if(islogy) gPad->SetLogy();

   TLegend *leg2 = new TLegend(0.15,0.65,0.3,0.9);
   leg2->AddEntry(htgth,"data","L");
   for(int ii=0; ii<nplane; ii++)
      leg2->AddEntry(htgth_MC[ii],Form("%s",plane_name[ii].Data()),"L");
   leg2->Draw();

   c2->cd(2);
   T_data->Draw("R.gold.ph>>htgph",TRK);
   htgph->SetLineColor(2);
   htgph->SetLineWidth(2);
   Double_t ntgph = htgph->Integral();
   for(int ii=0; ii<nplane; ii++){
      T_MC->Draw(Form("ph_%s_tr>>htgph_MC%d",plane_name[ii].Data(),ii+1),MCcut*"rate");
      Double_t ntgph_mc = htgph_MC[ii]->Integral();
      htgph_MC[ii]->Scale(ntgph/ntgph_mc);
      htgph_MC[ii]->SetLineColor(color[ii]); 
      htgph_MC[ii]->SetLineWidth(2); 
   }

   htgph_MC[nplane-1]->Draw("HIST");
   htgph->Draw("HIST same");
   for(int ii=0; ii<nplane-1; ii++)
     htgph_MC[ii]->Draw("HIST same");

   htgph_MC[nplane-1]->SetTitle(Form("%s target phi; ph;",tg.Data()));
   if(islogy) gPad->SetLogy();

   c2->cd(3);
   Double_t tgdp_min=0, tgdp_max=0;
   if(islogy) { tgdp_min = -1; tgdp_max = 1; }
   else { tgdp_min = -0.01; tgdp_max = 0.; }

   TH1F *htgdp = new TH1F("htgdp","L.gold.dp distgibution",100,tgdp_min,tgdp_max);
   TH1F *htgdp_MC[nplane];

   for(int ii=0; ii<nplane; ii++){
      htgdp_MC[ii] = new TH1F(Form("htgdp_MC%d",ii+1),"dp MC distgibution",100,tgdp_min,tgdp_max);
   }


   T_data->Draw("R.gold.dp>>htgdp",TRK);
   htgdp->SetLineColor(2);
   Double_t ntgdp = htgdp->Integral();

   for(int ii=0; ii<nplane; ii++){
      T_MC->Draw(Form("(p_%s_tr-950.0)/950.0>>htgdp_MC%d",plane_name[ii].Data(),ii+1),MCcut*"rate");
      Double_t ntgdp_mc = htgdp_MC[ii]->Integral();
      htgdp_MC[ii]->Scale(ntgdp/ntgdp_mc);
    
      htgdp_MC[ii]->SetLineColor(color[ii]); 
   }

   htgdp_MC[nplane-1]->Draw("HIST");
   htgdp->Draw("HIST same");
   for(int ii=0; ii<nplane-1; ii++)
     htgdp_MC[ii]->Draw("HIST same");

   htgdp_MC[nplane-1]->SetTitle(Form("%s dp; dp;",tg.Data()));
   if(islogy) gPad->SetLogy();

//=========  Q2 ============ //

   Double_t q2_min=0, q2_max=0;
   if(islogy) { q2_min = 0; q2_max = 0.1; }
   else { q2_min = 0.; q2_max = 0.015; }

   TH1F *hq2 = new TH1F("hq2","EK_R.Q2 distgibution",100,q2_min,q2_max);
   TH1F *hq2_MC[nplane];

   for(int ii=0; ii<nplane; ii++){
      hq2_MC[ii] = new TH1F(Form("hq2_MC%d",ii+1),"Q2 MC distgibution",100,q2_min,q2_max);
   }

   TCanvas *c3 = new TCanvas("c3","c3",1500,1500);
   T_data->Draw("EK_R.Q2>>hq2",TRK);
   hq2->SetLineColor(2);
   hq2->SetLineWidth(2); 
   Double_t nq2 = hq2->Integral();
   Double_t data_mean = hq2->GetMean();
   cout<<"data  "<<data_mean<<endl;

  Double_t ep_MC[nplane]={0},eth_MC[nplane]={0},eph_MC[nplane]={0};
  Double_t beamp=0.0, rate=0.0;

  T_MC->SetBranchAddress("ev.beamp",&beamp);
  T_MC->SetBranchAddress("rate",&rate);

  for(int ii=0; ii<nplane; ii++){
    T_MC->SetBranchAddress(Form("p_%s",plane_name[ii].Data()),&ep_MC[ii]);
    T_MC->SetBranchAddress(Form("ph_%s_tr",plane_name[ii].Data()),&eph_MC[ii]);
    T_MC->SetBranchAddress(Form("th_%s_tr",plane_name[ii].Data()),&eth_MC[ii]);
  }

   T_MC->Draw(">>goodMC",MCcut);
   TEventList *goodMC;
   gDirectory->GetObject("goodMC",goodMC);
   T_MC->SetEventList(goodMC);

   Double_t sep_angle = -4.74/180.*M_PI;
   Double_t nMC = goodMC->GetN();
   for(int ii=0; ii<nMC; ii++){
       T_MC->GetEntry(goodMC->GetEntry(ii));

       for(int jj=0; jj<nplane; jj++){
         Double_t tmp_th=0, tmpq2=0, tmp_tgth=0, tmp_tgph=0, tmp_tgp=0;
         tmp_tgth = eth_MC[jj];
         tmp_tgph = eph_MC[jj];
         tmp_tgp = ep_MC[jj];

         tmp_th = acos( (cos(sep_angle) - tmp_tgph*sin(sep_angle)) / (sqrt(1+pow(tmp_tgth,2)+pow(tmp_tgph,2)))  );
         tmpq2 = 4.*beamp*tmp_tgp/1000.* pow(sin(tmp_th/2.0),2);
         hq2_MC[jj]->Fill(tmpq2,rate);
       }
   }

   for(int ii=0; ii<nplane; ii++){
      Double_t tmp_mean = hq2_MC[ii]->GetMean();
      cout<<ii<<"  "<<tmp_mean<<endl;

      Double_t nq2_mc = hq2_MC[ii]->Integral();
      hq2_MC[ii]->Scale(nq2/nq2_mc);
      hq2_MC[ii]->SetLineColor(color[ii]); 
      hq2_MC[ii]->SetLineWidth(2); 
   }

   hq2_MC[nplane-1]->Draw("HIST");
   hq2->Draw("HIST same");
   for(int ii=0; ii<nplane-1; ii++)
     hq2_MC[ii]->Draw("HIST same");

   hq2_MC[nplane-1]->SetTitle(Form("%s Q2; Q2;",tg.Data()));
   if(islogy) gPad->SetLogy();

   TLegend *leg3 = new TLegend(0.15,0.65,0.3,0.9);
   leg3->AddEntry(hq2,"data","L");
   for(int ii=0; ii<nplane; ii++)
      leg3->AddEntry(hq2_MC[ii],Form("%s",plane_name[ii].Data()),"L");
   leg3->Draw();

   
   c1->Print(Form("outfiles/%s_RHRS.pdf[",tg.Data()));
   c1->Print(Form("outfiles/%s_RHRS.pdf",tg.Data()));
   c2->Print(Form("outfiles/%s_RHRS.pdf",tg.Data()));
   c3->Print(Form("outfiles/%s_RHRS.pdf",tg.Data()));
   c3->Print(Form("outfiles/%s_RHRS.pdf]",tg.Data()));





}
