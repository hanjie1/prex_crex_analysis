#include "LoadAsym.h"
#include "LoadACC.h"
#include "SetCut.h"

void CheckVertexAsym(){
     double ppb = 1e9;

     TChain *T = new TChain("T");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_1.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_2.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_3.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_4.root");
     T->Add("/w/halla-scifs17exp/parity/disk1/ryanrich/Ebeam953/TargNominal/Zero1_SandwichLHRS_PREX_0.0_5.root");

     TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,953);
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

     TH1F *hasym_origin = new TH1F("hasym_origin","vertex asymmetry in g4hrs with ACC cut",100,50,800);
     TH1F *hasym_acc = new TH1F("hasym_acc","vertex asymmetry calculated from acceptance function",100,50,800);
     TH1F *hasym_acc1 = new TH1F("hasym_acc1","vertex asymmetry calculated from acceptance function with dp cut",100,50,800);
     TH1F *hasym_e_acc = new TH1F("hasym_e_acc","elemental asymmetry with acceptance function applied",100,50,800);

     T->Draw("ev.A>>hasym_origin",ACC*"rate");   // vertex asymmetry calculated in g4hrs with acceptance cuts applied

     double thisBeamE,thisTh,thisRate,p_ztarg,xfp;
     T->SetBranchAddress("ev.Th",&thisTh);
     T->SetBranchAddress("ev.beamp",&thisBeamE);
     T->SetBranchAddress("rate",&thisRate);
     T->SetBranchAddress("p_ztarg",&p_ztarg);
     T->SetBranchAddress("x_fp_tr",&xfp);

     Double_t nentries = T->GetEntries();
     for(int ii=0; ii<nentries; ii++){
        T->GetEntry(ii);

	double thisAsym = Interpolate(thisBeamE,thisTh,0,1);  // vertex asymmetry (should be the same as ev.A)
	double thisACC = FindACC(thisTh,accp_angle,accp,accp_err,100);  // acceptance to the corresponding angle

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}
	hasym_acc->Fill(thisAsym*ppb,thisRate*thisACC);  // vertex asymmetry with the acceptance function applied

	if( (p_peak-p_ztarg)<2.2 ) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
     }

     double setBeamE = 0.953; 
     double dtheta = (accp_angle[1]-accp_angle[0])/180.*TMath::Pi();  // delta theta in radius
     for(int ii=0; ii<100; ii++){
	double thisAsym = Interpolate(setBeamE,accp_angle[ii],0,1);  // elemental asymmetry
	double thisXS = Interpolate(setBeamE,accp_angle[ii],0,0);  // elemental cross section
	double thisACC = 0;

	if(accp[ii]<accp_err[ii]) thisACC=0;
	else thisACC = accp[ii];
cout<<ii<<"  "<<accp_angle[ii]<<"  "<<thisAsym<<"  "<<thisXS<<"  "<<thisACC<<endl;
	double omega = TMath::Sin(accp_angle[ii]/180.*TMath::Pi())*dtheta;
	hasym_e_acc->Fill(thisAsym*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
     } 

     cout<<"================================="<<endl;
     cout<<"orignal vertex asymmetry mean:                "<<setw(8)<<hasym_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_origin->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean:              "<<setw(8)<<hasym_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry with dp cut mean:  "<<setw(8)<<hasym_acc1->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc1->GetRMS()<<endl;
     cout<<"elemental asymmetry mean:                     "<<setw(8)<<hasym_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;

     Double_t inte1 = hasym_origin->GetMaximum();
     Double_t inte2 = hasym_e_acc->GetMaximum();
     hasym_e_acc->Scale(inte1/inte2);

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

     TLegend *leg = new TLegend(0.25,0.65,0.4,0.8);
     leg->AddEntry(hasym_origin,"origin","L");
     leg->AddEntry(hasym_acc,"reproduce","L");
     leg->AddEntry(hasym_acc1,"reproduce (dp cut)","L");
     leg->AddEntry(hasym_e_acc,"elemental","L");
     leg->Draw();
    

}
