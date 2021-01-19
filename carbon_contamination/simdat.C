#include "CollimatorL.C"

void simdat(){

TChain *T[2];
TH1F *x[2], *y[2], *th[2], *ph[2];
TH2F *acc[2];
TH1F *qsq[2];
TH1F *thtg[2], *phtg[2];


gStyle->SetOptStat(0);

x[0] = new TH1F("x[0]","X distribution at Detector Plane",100,-0.8,0.1);
y[0] = new TH1F("y[0]", "Y distribution at Detector Plane",100,-0.05,0.05);
th[0] = new TH1F("th[0]","Theta distribution at Detector Plane", 100,-0.1,0.03);
ph[0] = new TH1F("ph[0]","Phi distribution at Detector Plane",100,-0.045,0.045);
acc[0] = new TH2F("acc[0]","Theta-Phi Acceptance",100,-0.04,0.04,100,-0.05,0.05);

thtg[0] = new TH1F("thtg[0]","Theta Target",100,-0.06,0.06);
phtg[0] = new TH1F("phtg[0]","Phi Target",100,-0.03,0.03);


qsq[0] = new TH1F("qsq[0]","Q^{2}",100,0,0.012);
qsq[1] = new TH1F("qsq[1]","Q^{2}",100,0,0.012);


x[1] = new TH1F("x[1]","X distribution at Detector Plane",100,-0.8,0.1);
y[1] = new TH1F("y[1]", "Y distribution at Detector Plane",100,-0.05,0.05);
th[1] = new TH1F("th[1]","Theta distribution at Detector Plane", 100,-0.1,0.03);
ph[1] = new TH1F("ph[1]","Phi distribution at Detector Plane",100,-0.045,0.045);
acc[1] = new TH2F("acc[1]","Theta-Phi Acceptance at Target",150,-0.15,0.15,150,-0.15,0.15);

thtg[1] = new TH1F("thtg[1]","Theta Target",100,-0.06,0.06);
phtg[1] = new TH1F("phtg[1]","Phi Target",100,-0.03,0.03);




for(int i = 0; i < 2; i++) { T[i] = new TChain("T");

acc[i]->GetYaxis()->SetTitle("#theta_{tg} (rad)");
acc[i]->GetXaxis()->SetTitle("#phi_{tg} (rad)");
qsq[i]->GetYaxis()->SetTitle("Q^{2} (GeV^{2}/c)");

x[i]->GetXaxis()->SetTitle("x (m)");
y[i]->GetXaxis()->SetTitle("y (m)");
th[i]->GetXaxis()->SetTitle("#theta (rad)");
ph[i]->GetXaxis()->SetTitle("#phi (rad)");

thtg[i]->GetXaxis()->SetTitle("#theta_{tg} (rad)");
phtg[i]->GetXaxis()->SetTitle("#phi_{tg} (rad)");


if( i == 0 ) {  
x[i]->SetLineColor(kBlack); y[i]->SetLineColor(kBlack); th[i]->SetLineColor(kBlack); ph[i]->SetLineColor(kBlack); acc[i]->SetMarkerColor(kBlack); qsq[i]->SetLineColor(kBlack);thtg[i]->SetLineColor(kBlack); phtg[i]->SetLineColor(kBlack);  }
else {  x[i]->SetLineColor(kRed); y[i]->SetLineColor(kRed); th[i]->SetLineColor(kRed); ph[i]->SetLineColor(kRed); acc[i]->SetMarkerColor(kRed); qsq[i]->SetLineColor(kRed); phtg[i]->SetLineColor(kRed); thtg[i]->SetLineColor(kRed); }
 }


T[0]->Add("prexLHRS_2052_50000.root");

double Q1[5] = { 0.225106, 0.224926,  0.223931, 0.224715, 0.224802 };

//for(int i = 1; i < 6; i++) { T[1]->Add(Form("../septmistune/Rootfiles/PREXfiles/septscaleLHRS_PREX_0.0_%i.root",i)); } 

for(int i = 1; i < 6; i++) { T[1]->Add(Form("../HRSRate/rootfiles/RateLHRS_%g_%i.root",Q1[0],i)); }

T[1]->SetAlias("den","sqrt(1 + th_tg_tr**2 + ph_tg_tr**2)");
T[1]->SetAlias("num","cos(0.08476592938) - ph_tg_tr*sin(0.08476592938)");//+1 for R, -1 for L

T[1]->SetAlias("Cos","num/den");

T[1]->SetAlias("var","1-Cos");

T[1]->SetAlias("try","2*0.95*0.95*var");


char ctrig[50],vdccut[200],tgtcut[200], tgtcut1[200], qtz1[200], ccut1[800],ccut[800], qtz[200], cand[5];

sprintf(ctrig,"fEvtHdr.fEvtType==1");
sprintf(vdccut,"L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1");
sprintf(cand,"&&");
sprintf(tgtcut,"L.gold.th[0]>-0.055&&L.gold.th[0]<0.055&&L.gold.ph[0]>-0.018&&L.gold.ph[0]<0.026");
sprintf(qtz,"P.loQadcL>600");
sprintf(qtz1,"P.upQadcL>490");
sprintf(tgtcut1,"L.gold.th[0]>-0.08&&L.gold.th[0]<0.08&&L.gold.ph[0]>-0.05&&L.gold.ph[0]<0.05");


strcpy(ccut1,ctrig);
strcat(ccut1,cand);
strcat(ccut1,vdccut);
strcat(ccut1,cand);
strcat(ccut1,tgtcut1);
strcat(ccut1,cand);
strcat(ccut1,qtz);
strcat(ccut1,cand);
strcat(ccut1,qtz1);

strcpy(ccut,ctrig);
strcat(ccut,cand);
strcat(ccut,vdccut);
strcat(ccut,cand);
strcat(ccut,tgtcut);
//strcat(ccut,cand);
//strcat(ccut,qtz);

double xmin[9] = { 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[9] = {0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[9] = {-0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[9] = {0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };


TCut totalCut,totalCut1;
TCut colCut = "CollimatorL(x_col_tr,y_col_tr)";
TCut radCut1 = Form("x_zup1 > %f && x_zup1 < %f && y_zup1 > %f && y_zup1 < %f",0.0448,0.0707,-0.0258,0.0255);
TCut radCut2 = Form("x_zup2 > %f && x_zup2 < %f && y_zup2 > %f && y_zup2 < %f",0.0482,0.0757,-0.0275,0.0272);
TCut downCut1 = Form("x_zdown1 > %f && x_zdown1 < %f && y_zdown1 > %f && y_zdown1 < %f",xmin[0],xmax[0],ymin[0],ymax[0]);
TCut downCut2 = Form("x_zdown2 > %f && x_zdown2 < %f && y_zdown2 > %f && y_zdown2 < %f",xmin[1],xmax[1],ymin[1],ymax[1]);
TCut downCut3 = Form("x_zdown3 > %f && x_zdown3 < %f && y_zdown3 > %f && y_zdown3 < %f",xmin[2],xmax[2],ymin[2],ymax[2]);
TCut downCut4 = Form("x_zdown4 > %f && x_zdown4 < %f && y_zdown4 > %f && y_zdown4 < %f",xmin[3],xmax[3],ymin[3],ymax[3]);
TCut downCut5 = Form("x_zdown5 > %f && x_zdown5 < %f && y_zdown5 > %f && y_zdown5 < %f",xmin[4],xmax[4],ymin[4],ymax[4]);
TCut downCut6 = Form("x_zdown6 > %f && x_zdown6 < %f && y_zdown6 > %f && y_zdown6 < %f",xmin[5],xmax[5],ymin[5],ymax[5]);
TCut downCut7 = Form("x_zdown7 > %f && x_zdown7 < %f && y_zdown7 > %f && y_zdown7 < %f",xmin[6],xmax[6],ymin[6],ymax[6]);
TCut downCut8 = Form("x_zdown8 > %f && x_zdown8 < %f && y_zdown8 > %f && y_zdown8 < %f",xmin[7],xmax[7],ymin[7],ymax[7]);
TCut downCut9 = Form("x_zdown9 > %f && x_zdown9 < %f && y_zdown9 > %f && y_zdown9 < %f",xmin[8],xmax[8],ymin[8],ymax[8]);
TCut xcut = "x_fp_tr!=-333.";
TCut momCut = "(p_col_tr - 950)/950 > -0.002736842";
TCut ph_cut = "(ph_tg_tr>-0.014&&ph_tg_tr<0.026)"; 

totalCut = ph_cut&&colCut&&xcut&&radCut1&&radCut2&&downCut1&&downCut2&&downCut3&&downCut4&&downCut5&&downCut6&&downCut7&&downCut8;

totalCut1 = totalCut&&momCut;


TCanvas *c = new TCanvas("c","c");
c->Divide(2,2);
c->cd(1);
T[0]->Draw("(L.tr.x[0]+0.9*L.tr.th[0])+0.01 >> x[0]",ccut);
x[0]->Scale(1.0/x[0]->GetBinContent(x[0]->GetMaximumBin()));
T[1]->Draw("x_fp_tr >> x[1]",totalCut*"rate","same");
x[1]->Scale(1.0/x[1]->GetBinContent(x[1]->GetMaximumBin()));
c->cd(2);
T[0]->Draw("(L.tr.y[0]+0.9*L.tr.ph[0])+0.02 >> y[0]",ccut);
y[0]->Scale(1.0/y[0]->GetBinContent(y[0]->GetMaximumBin()));
T[1]->Draw("y_fp_tr >> y[1]",totalCut*"rate","same");
y[1]->Scale(1.0/y[1]->GetBinContent(y[1]->GetMaximumBin()));
c->cd(3);
T[0]->Draw("L.tr.th[0] >> th[0]",ccut);
th[0]->Scale(1.0/th[0]->GetBinContent(th[0]->GetMaximumBin()));
T[1]->Draw("th_fp_tr >> th[1]",totalCut*"rate","same");
th[1]->Scale(1.0/th[1]->GetBinContent(th[1]->GetMaximumBin()));
c->cd(4);
T[0]->Draw("L.tr.ph[0]+0.002 >> ph[0]",ccut);
ph[0]->Scale(1.0/ph[0]->GetBinContent(ph[0]->GetMaximumBin()));
T[1]->Draw("ph_fp_tr >> ph[1]",totalCut*"rate","same");
ph[1]->Scale(1.0/ph[1]->GetBinContent(ph[1]->GetMaximumBin()));
//c->Print("focalPlane.png");



TCanvas *c1 = new TCanvas("c1","c1");
c1->cd();
//T[1]->Draw("th_tg_tr:ph_tg_tr >> acc[1]",totalCut1);
//T[1]->Draw("th_tg_tr:ph_tg_tr >> acc[1]",totalCut1);
T[1]->Draw("x_col_tr:y_col_tr >> acc[1]",totalCut1);
//T[0]->Draw("L.gold.th[0]:L.gold.ph[0] >> acc[0]",ccut1,"same");
//c1->Print("acc2D_n.png");


/*
TCanvas *c2 = new TCanvas("c2","c2");
c2->cd();
T[0]->Draw("EK_L.Q2[0] >> qsq[0]",ccut1);
qsq[0]->Scale(1.0/qsq[0]->GetBinContent(qsq[0]->GetMaximumBin()));
T[1]->Draw("try >> qsq[1]",totalCut1*"rate/3500","same");
qsq[1]->Scale(1.0/qsq[1]->GetBinContent(qsq[1]->GetMaximumBin()));
c2->Print("qsqcomp.png");
*/


TCanvas *c3 = new TCanvas("c3","c3");
c3->Divide(2,1);
c3->cd(1);
T[0]->Draw("L.gold.th[0] >> thtg[0]",ccut1);
thtg[0]->Scale(1.0/thtg[0]->GetBinContent(thtg[0]->GetMaximumBin()));
T[1]->Draw("th_tg_tr >> thtg[1]",totalCut1,"same");
thtg[1]->Scale(1.0/thtg[1]->GetBinContent(thtg[1]->GetMaximumBin()));
c3->cd(2);
T[0]->Draw("L.gold.ph[0] >> phtg[0]",ccut1);
phtg[0]->Scale(1.0/phtg[0]->GetBinContent(phtg[0]->GetMaximumBin()));
T[1]->Draw("ph_tg_tr >> phtg[1]",totalCut1,"same");
phtg[1]->Scale(1.0/phtg[1]->GetBinContent(phtg[1]->GetMaximumBin()));
//c3->Print("onedacc.png");


//T[1]->Draw("ph_tg_tr >> phtg[1]",totalCut1);


} 
