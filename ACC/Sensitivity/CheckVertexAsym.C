#include "LoadAsym.h"
#include "LoadACC.h"
#include "SetCut.h"

void CheckVertexAsym(){
     double ppb = 1e9;
     double amu_c2 = 0.931494028; // GeV

     TChain *T = new TChain("T");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("../Rootfiles/Zero1_SandwichLHRS_PREX_0.0_5.root");

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     hpztarg->Draw("HIST");

     Int_t pbin = hpztarg->GetMaximumBin();
     Double_t p_peak = hpztarg->GetBinCenter(pbin);
     cout<<"p peak:  "<<p_peak<<endl;

     // dp cut
     TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     TCut DP = Form("%s",dpcut.Data());

     // final cuts
     TCut ACC = DP+isPb+colCut+XCUT;

     // load acceptance table
     double accp_angle[100]={0}, accp[100]={0}, accp_err[100]={0};
     int status = LoadACC(accp_angle, accp, accp_err);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     // load asymmetry table
     LoadTable("horpb.dat", 0);

     int nbins=200;
     TH1F *hasym_origin = new TH1F("hasym_origin","vertex asymmetry in g4hrs with ACC cut",nbins,200,1000);
     TH1F *hasym_acc = new TH1F("hasym_acc","vertex asymmetry calculated from acceptance function",nbins,200,1000);
     TH1F *hasym_acc1 = new TH1F("hasym_acc1","vertex asymmetry calculated from acceptance function with dp cut",nbins,200,1000);
     TH1F *hasym_e_acc = new TH1F("hasym_e_acc","elemental asymmetry with acceptance function applied",nbins,200,1000);

     TH1F *hth_origin = new TH1F("hth_origin","vertex theta in g4hrs with ACC cut",100,3,8);
     TH1F *hth_e_acc = new TH1F("hth_e_acc","elemental theta with acceptance function applied",100,3,8);

     TH1F *hq2_origin = new TH1F("hq2_origin","vertex Q2 in g4hrs with ACC cut",100,0,0.015);
     TH1F *hq2_e_acc = new TH1F("hq2_e_acc","elemental Q2 with acceptance function applied",100,0,0.015);

     T->Draw("ev.A>>hasym_origin",ACC*"rate");   // vertex asymmetry calculated in g4hrs with acceptance cuts applied
     T->Draw("ev.Th>>hth_origin",ACC*"rate");   // vertex theta calculated in g4hrs with acceptance cuts applied
     T->Draw("ev.Q2>>hq2_origin",ACC*"rate");   // vertex Q2 calculated in g4hrs with acceptance cuts applied

     double thisBeamE,thisTh,thisRate,p_ztarg,xfp;
     int thisA;
     T->SetBranchAddress("ev.Th",&thisTh);
     T->SetBranchAddress("ev.beamp",&thisBeamE);
     T->SetBranchAddress("ev.nuclA",&thisA);
     T->SetBranchAddress("rate",&thisRate);
     T->SetBranchAddress("p_ztarg",&p_ztarg);
     T->SetBranchAddress("x_fp_tr",&xfp);

     Double_t nentries = T->GetEntries();
     for(int ii=0; ii<nentries; ii++){
        T->GetEntry(ii);

	double thisAsym = Interpolate(thisBeamE,thisTh,0,1);  // vertex asymmetry (should be the same as ev.A)
	double thisACC = FindACC(thisTh,accp_angle,accp,accp_err,100);  // acceptance to the corresponding angle

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}
	if(thisA==208 && thisBeamE>p_peak/1000.0)hasym_acc->Fill(thisAsym*ppb,thisRate*thisACC);  // vertex asymmetry with the acceptance function applied

	//if( (p_peak-p_ztarg)<2.2 ) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
	if( thisA==208 && (p_peak-p_ztarg)<2.2) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
     }

     double setBeamE = 0.9534; 
     int nbin_th = 10000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius

     for(int ii=0; ii<nbin_th; ii++){
	double thisAngle = 3.0 + dtheta*ii;
	double thisAsym = Interpolate(setBeamE,thisAngle,0,1);  // elemental asymmetry
	double thisXS = Interpolate(setBeamE,thisAngle,0,0);  // elemental cross section

	double deg_to_rad = TMath::Pi()/180.;
	double thisEp = setBeamE/(1.+setBeamE/(208.*amu_c2)*(1.-TMath::Cos(thisAngle*deg_to_rad)));
	double thisQ2 = 4.*setBeamE*thisEp*pow(TMath::Sin(thisAngle*deg_to_rad/2.),2);

	double thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,100);  // acceptance to the corresponding angle
	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

//cout<<ii<<"  "<<accp_angle[ii]<<"  "<<thisEp<<"  "<<thisQ2<<"  "<<thisAsym*ppb<<"  "<<thisXS<<"  "<<thisACC<<endl;
//
	double omega = TMath::Sin(thisAngle*deg_to_rad)*dtheta*deg_to_rad;
	hasym_e_acc->Fill(thisAsym*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
	hth_e_acc->Fill(thisAngle,thisXS*omega*thisACC);   // elemetal angle with the acceptance function applied
	hq2_e_acc->Fill(thisQ2,thisXS*omega*thisACC);   // elemetal Q2 with the acceptance function applied
     } 

     cout<<"================================="<<endl;
     cout<<"orignal vertex asymmetry mean:     "<<setw(8)<<hasym_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_origin->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean:              "<<setw(8)<<hasym_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean(with dp cut):   "<<setw(8)<<hasym_acc1->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc1->GetRMS()<<endl;
     cout<<"elemental asymmetry mean:          "<<setw(8)<<hasym_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex angle mean:   "<<setw(8)<<hth_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_origin->GetRMS()<<endl;
     cout<<"elemental angle mean:        "<<setw(8)<<hth_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex Q2 mean:   "<<setw(8)<<hq2_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_origin->GetRMS()<<endl;
     cout<<"elemental Q2 mean:        "<<setw(8)<<hq2_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;

     Double_t inte1 = hasym_origin->GetMaximum();
     Double_t inte2 = hasym_e_acc->GetMaximum();
     Double_t inte3 = hasym_acc1->GetMaximum();
     hasym_e_acc->Scale(inte1/inte2);
     hasym_acc1->Scale(inte1/inte3);

     TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
     hasym_origin->SetLineColor(1);
     hasym_acc->SetLineColor(4);
     hasym_acc1->SetLineColor(6);
     hasym_e_acc->SetLineColor(2);

     hasym_origin->SetLineWidth(2);
     hasym_acc->SetLineWidth(2);
     hasym_acc1->SetLineWidth(2);
     hasym_e_acc->SetLineWidth(2);

     hasym_origin->Draw("HIST");    
     hasym_acc->Draw("HIST same");    
     hasym_acc1->Draw("HIST same");    
     hasym_e_acc->Draw("HIST same");
     hasym_origin->SetTitle("Asymmetry with acceptance; asym(ppb)");

     TLegend *leg = new TLegend(0.15,0.65,0.35,0.8);
     leg->AddEntry(hasym_origin,"origin","L");
     leg->AddEntry(hasym_acc,"reproduce","L");
     leg->AddEntry(hasym_acc1,"reproduce (dp and pb cut)","L");
     leg->AddEntry(hasym_e_acc,"elemental","L");
     leg->Draw();
    
     inte1 = hth_origin->GetMaximum();
     inte2 = hth_e_acc->GetMaximum();
     hth_e_acc->Scale(inte1/inte2);

     TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
     c2->Divide(2,1);
     c2->cd(1);
     hth_origin->SetLineColor(1);
     hth_e_acc->SetLineColor(2);

     hth_origin->SetLineWidth(2);
     hth_e_acc->SetLineWidth(2);

     hth_origin->Draw("HIST");    
     hth_e_acc->Draw("HIST same");
     hth_origin->SetTitle("angle with acceptance; angle");

     TLegend *leg1 = new TLegend(0.65,0.65,0.8,0.8);
     leg1->AddEntry(hth_origin,"origin","L");
     leg1->AddEntry(hth_e_acc,"elemental","L");
     leg1->Draw();
    
     c2->cd(2);
     inte1 = hq2_origin->GetMaximum();
     inte2 = hq2_e_acc->GetMaximum();
     hq2_e_acc->Scale(inte1/inte2);

     hq2_origin->SetLineColor(1);
     hq2_e_acc->SetLineColor(2);

     hq2_origin->SetLineWidth(2);
     hq2_e_acc->SetLineWidth(2);

     hq2_origin->Draw("HIST");    
     hq2_e_acc->Draw("HIST same");
     hq2_origin->SetTitle("Q2 with acceptance; Q2");

     TLegend *leg2 = new TLegend(0.15,0.65,0.35,0.8);
     leg2->AddEntry(hq2_origin,"origin","L");
     leg2->AddEntry(hq2_e_acc,"elemental","L");
     leg2->Draw();
    

}
