#include "CollimatorR.C"

void ComparePbC(){
   TChain *T_pb = new TChain("T");  // pb208
   TChain *T_c = new TChain("T");   // c12

   T_pb->Add("../g4hrs_rootfiles/prex_pb208_RHRS.root");
   T_c->Add("../g4hrs_rootfiles/prex_c12_RHRS.root");

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

   TH1F *hpb_fpx = new TH1F("hpb_fpx","pb208 focal plane x",100,-0.6,0.1);	
   TH1F *hc_fpx = new TH1F("hc_fpx","c12 focal plane x",100,-0.6,0.1);	

   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
   T_pb->Draw("x_fp_tr>>hpb_fpx",MCcut*"rate");
   T_c->Draw("x_fp_tr>>hc_fpx",MCcut*"rate","same");

   int binmax_pb = hpb_fpx->GetMaximumBin(); 
   double ymax_pb = hpb_fpx->GetBinContent(binmax_pb);
   int binmax_c = hc_fpx->GetMaximumBin(); 
   double ymax_c = hc_fpx->GetBinContent(binmax_c);

   cout<<ymax_pb<<"  "<<ymax_c<<"  "<<hpb_fpx->Integral()<<"  "<<hc_fpx->Integral()<<endl;

   hc_fpx->Scale(ymax_pb/ymax_c);

   gStyle->SetOptStat(0);

   hpb_fpx->Draw("HIST");
   hc_fpx->Draw("HIST same");
   hpb_fpx->SetTitle("RHRS focal plane x; x_fp_tr;");
   hpb_fpx->SetLineColor(4);
   hpb_fpx->SetLineWidth(2);
   hc_fpx->SetLineColor(2);
   hc_fpx->SetLineWidth(2);

   TLegend *leg = new TLegend(0.2,0.7,0.35,0.85);
   leg->AddEntry(hpb_fpx,"Pb208","L");
   leg->AddEntry(hc_fpx,"C12","L");
   leg->Draw();
}
