#include "CollimatorL.C"

void CompareDp(){
   
   bool islogy = false;

   TChain *T_data = new TChain("T");
   TChain *T_MC = new TChain("T");

   TString target;
   TString tg;
   cout<<"Which target? (Pb, C12, Ca40) "<<endl;
   cin>>target;
   int run_number=0;
   if(target == "Pb") { run_number = 2055; tg = "pb208"; }
   if(target == "C12") { run_number = 2302; tg = "c12"; }
   if(target == "Ca40") { run_number = 2297; tg = "ca40"; }

   TString datafile = Form("/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexLHRS_%d_-1.root",run_number);
   TString MCfile = Form("g4hrs_rootfiles/prex_%s.root",tg.Data()); 

   T_data->Add(datafile);
   T_MC->Add(MCfile);

/*** data cut ***/
   TCut TRK = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";

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
 
 
   TCut MCcut = colCut+radCut1+radCut2+downCut1+downCut2+downCut3+downCut4+downCut5+downCut6+downCut7+downCut8+downCut9+xcut;

//======== target ===========//

   Double_t tgth_min=0.0,tgth_max=0,tgph_min=0,tgph_max=0;
   if(islogy){
	tgth_min = -0.5; tgth_max=0.5; tgph_min=-0.5; tgph_max=0.5;
   }
   else{
	tgth_min = -0.08; tgth_max=0.08; tgph_min=-0.04; tgph_max=0.04;
   }


   Double_t tgdp_min=0, tgdp_max=0;
   if(islogy) { tgdp_min = -1; tgdp_max = 1; }
   else { tgdp_min = -0.008; tgdp_max = 0.; }

   TH1F *htgdp = new TH1F("htgdp","L.gold.dp distgibution",100,tgdp_min,tgdp_max);
   TH1F *htgdp_MC[14];

   TString plane_name[14]= {"sen","sm","sex","col","q1en","q1ex","q2en","q2ex","den","dex","q3en","q3ex","vdc","fp"};

   for(int ii=0; ii<14; ii++){
      htgdp_MC[ii] = new TH1F(Form("htgdp_MC%d",ii+1),"dp MC distgibution",100,tgdp_min,tgdp_max);
   }


   T_data->Draw("L.gold.dp>>htgdp",TRK);
   htgdp->SetLineColor(2);
   Double_t ntgdp = htgdp->Integral();

   int color[14]={1,3,4,5,6,7,8,9,30,40,46,38,41};
   for(int ii=0; ii<14; ii++){
      T_MC->Draw(Form("(p_%s_tr-950.0)/950.0>>htgdp_MC%d",plane_name[ii].Data(),ii+1),MCcut*"rate");
      Double_t ntgdp_mc = htgdp_MC[ii]->Integral();
cout<<ntgdp_mc<<endl;
      htgdp_MC[ii]->Scale(ntgdp/ntgdp_mc);
    
      htgdp_MC[ii]->SetLineColor(color[ii]); 
   }

   gStyle->SetOptStat(0);

   TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
   htgdp_MC[0]->Draw("HIST");
   htgdp->Draw("HIST same");
   for(int ii=1; ii<14; ii++)
     htgdp_MC[ii]->Draw("HIST same");

   htgdp_MC[0]->SetTitle(Form("%s dp; dp;",tg.Data()));
   if(islogy) gPad->SetLogy();

   TLegend *leg1 = new TLegend(0.15,0.75,0.25,0.8);
   leg1->AddEntry(htgdp,"data","L");
   leg1->Draw();



}
