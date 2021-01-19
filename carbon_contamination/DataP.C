#include "CollimatorL.C"

void DataP(){
   TChain *T_c = new TChain("T");
   TChain *T_pb = new TChain("T");
   TChain *T_ca = new TChain("T");

   TString datafile1 = "/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexLHRS_2302_-1_file0.root";
   TString datafile2 = "/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexLHRS_2055_-1_file0.root";
   TString datafile3 = "/lustre/expphy/volatile/halla/parity/ryanrich/prexRootFiles/prexLHRS_2297_-1_file0.root";

   T_c->Add(datafile1);
   T_pb->Add(datafile2);
   T_ca->Add(datafile3);

   TString MCfile1 = "g4hrs_rootfiles/prex_c12_wIon.root";
   TString MCfile2 = "g4hrs_rootfiles/prex_pb208_wIon.root";
   TString MCfile3 = "g4hrs_rootfiles/prex_ca40.root";

   TChain *TMC_c = new TChain("T");
   TChain *TMC_pb = new TChain("T");
   TChain *TMC_ca = new TChain("T");

   TMC_c->Add(MCfile1);
   TMC_pb->Add(MCfile2);
   TMC_ca->Add(MCfile3);


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

   // =====  p central  ===== //

   TH1F *hp_c = new TH1F("hp_c","momentum for c12",100,940,952);
   TH1F *hp_pb = new TH1F("hp_pb","momentum for pb208",100,940,952);
   TH1F *hp_ca = new TH1F("hp_ca","momentum for ca40",100,940,952);

   T_c->Draw("HacL_D1_NMR_SIG>>hNMR_c");
   TH1F *hNMR_c = (TH1F *) gDirectory->Get("hNMR_c");
   Double_t NMR_c = hNMR_c->GetMean();
   Double_t pcentral_c = 2.702*(NMR_c)-1.6e-03*pow(NMR_c,3);

   T_pb->Draw("HacL_D1_NMR_SIG>>hNMR_pb");
   TH1F *hNMR_pb = (TH1F *) gDirectory->Get("hNMR_pb");
   Double_t NMR_pb = hNMR_pb->GetMean();
   Double_t pcentral_pb = 2.702*(NMR_pb)-1.6e-03*pow(NMR_pb,3);

   T_ca->Draw("HacL_D1_NMR_SIG>>hNMR_ca");
   TH1F *hNMR_ca = (TH1F *) gDirectory->Get("hNMR_ca");
   Double_t NMR_ca = hNMR_ca->GetMean();
   Double_t pcentral_ca = 2.702*(NMR_ca)-1.6e-03*pow(NMR_ca,3);

   cout<<"p central"<<endl;
   cout<<"c12:    "<<pcentral_c<<endl;
   cout<<"pb208:  "<<pcentral_pb<<endl;
   cout<<"c40:    "<<pcentral_ca<<endl;

   T_c->Draw("(1.+L.gold.dp)*(2.702*(HacL_D1_NMR_SIG)-1.6e-03*pow(HacL_D1_NMR_SIG,3))*1000>>hp_c",TRK);
   T_pb->Draw("(1.+L.gold.dp)*(2.702*(HacL_D1_NMR_SIG)-1.6e-03*pow(HacL_D1_NMR_SIG,3))*1000>>hp_pb",TRK);
   T_ca->Draw("(1.+L.gold.dp)*(2.702*(HacL_D1_NMR_SIG)-1.6e-03*pow(HacL_D1_NMR_SIG,3))*1000>>hp_ca",TRK);

   TH1F *hpMC_c = new TH1F("hpMC_c","MC momentum for c12",100,940,952);
   TH1F *hpMC_pb = new TH1F("hpMC_pb","MC momentum for pb208",100,940,952);
   TH1F *hpMC_ca = new TH1F("hpMC_ca","MC momentum for ca40",100,940,952);

   TMC_c->Draw("p_zsieve>>hpMC_c",MCcut*"rate");
   TMC_pb->Draw("p_zsieve>>hpMC_pb",MCcut*"rate");
   TMC_ca->Draw("p_zsieve>>hpMC_ca",MCcut*"rate");

   Double_t ndata_c = hp_c->Integral();  
   Double_t nMC_c = hpMC_c->Integral();  
   hpMC_c->Scale(ndata_c/nMC_c);

   Double_t ndata_pb = hp_pb->Integral();  
   Double_t nMC_pb = hpMC_pb->Integral();  
   hpMC_pb->Scale(ndata_pb/nMC_pb);

   Double_t ndata_ca = hp_ca->Integral();  
   Double_t nMC_ca = hpMC_ca->Integral();  
   hpMC_ca->Scale(ndata_ca/nMC_ca);

   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
   c1->Divide(2,2);
   c1->cd(1);
   hpMC_c->Draw("HIST");
   hp_c->Draw("HIST same");
   hp_c->SetLineWidth(2);
   hpMC_c->SetLineWidth(2);
   hp_c->SetLineColor(2);
   hpMC_c->SetLineColor(4);
   hpMC_c->SetTitle("LHRS C12 momentum (wIon); p (MeV);");
 
   c1->cd(2);
   hpMC_pb->Draw("HIST");
   hp_pb->Draw("HIST same");
   hp_pb->SetLineWidth(2);
   hpMC_pb->SetLineWidth(2);
   hp_pb->SetLineColor(2);
   hpMC_pb->SetLineColor(4);
   hpMC_pb->SetTitle("LHRS Pb208 momentum (wIon); p (MeV);");

   c1->cd(3);
   hp_ca->Draw("HIST");
   hpMC_ca->Draw("HIST same");
   hp_ca->SetLineWidth(2);
   hpMC_ca->SetLineWidth(2);
   hp_ca->SetLineColor(2);
   hpMC_ca->SetLineColor(4);
   hp_ca->SetTitle("LHRS Ca40 momentum; p (MeV);");
 
   c1->Print("outfiles/momentum_LHRS_wIon.pdf");

}
