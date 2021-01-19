#include "LoadAsym.h"
#include "LoadACC.h"
#include "../../ryanscript/CollimatorL.C"
#include "../../ryanscript/UpPlane.C"
#include "../../ryanscript/DownPlane.C"


void kCheckAccL() {

  double sept = 0.0;
  int iver = 0;
  double pinch = 0.0023;
  double colpinch = 0.000;
  double E0 = 953.4;
  TString stub = Form("LHRS_s%.1f_v%i",sept,iver);
  TString accfun = "./accfunc_"+stub+".csv";
  int iprint =1;
  
     double ppb = 1e9;
     double amu_c2 = 0.931494028; // GeV

     TChain *T = new TChain("T");
     for(int i = 1; i < 6; i++){ T->Add(Form("../../sim953/Zero1_SandwichLHRS_PREX_%.1f_%i.root",sept,i)); }


     TString xcut = "x_fp_tr!=-333.";
     TString colcut = Form("CollimatorL(x_col_tr,y_col_tr+%f)",colpinch);
     TString upcut = "UpPlane(x_zup1,y_zup1,x_zup2,y_zup2,1)";
     TString dncut = Form(" DownPlane(x_zdown1,y_zdown1,x_zdown2,y_zdown2,x_zdown3-%f,y_zdown3,x_zdown4-%f,y_zdown4,x_zdown5,y_zdown5,x_zdown6,y_zdown6,x_zdown7,y_zdown7,x_zdown8,y_zdown8,x_zdown9,y_zdown9,1)",pinch,pinch);

     TCut XCUT = Form("%s",xcut.Data());//&&x_vdc_tr>0.0";
     TCut colCut = Form("%s",colcut.Data());
     TCut upCut = Form("%s",upcut.Data());
     TCut dnCut = Form("%s",dncut.Data());
          
     // lead cut
     TString ispb = "ev.nuclA==208";
     TCut isPb = Form("%s",ispb.Data());    

     double p_peak = E0-1.0;
     double simmomcut = 2.2;
    //   cut th_tg = 0 at threshold below the expected elastic peak
     TString pcut = Form("(p_ztarg_tr - 50*th_ztarg_tr) > (%5.1f - %5.1f)",p_peak,simmomcut);
     TString incpcut = Form("(1000*ev.beamp) > (%5.1f - 5)",p_peak);
     // cout << pcut.Data() << endl;
     // cout << incpcut.Data() << endl;
     // scmom = (thismom - 50*thisTh) > ((E0 - 0.001 - simMomCut)*1000);
    // scincp = (beamP - 50*thisTh) > ((E0 - 0.001 - simMomCut)*1000);

     TCut DP = Form("%s",pcut.Data());
     TCut ECUT = Form("%s",incpcut.Data());


     
     // TCanvas *tmp = new TCanvas();
     // TH1F *hpztarg = new TH1F("hpztarg","p_ztarg distribution with cut",200,945,954);
     // T->Draw("p_ztarg>>hpztarg",(colCut+isPb+XCUT)*"rate");
     // hpztarg->Draw("HIST");

     // Int_t pbin = hpztarg->GetMaximumBin();
     // Double_t p_peak = hpztarg->GetBinCenter(pbin);
     // cout<<"p peak:  "<<p_peak<<endl;

     // // dp cut for accepted events
     // TString dpcut = Form("(%f-p_ztarg)<2.2",p_peak);
     // TCut DP = Form("%s",dpcut.Data());

     // // dp cut for the incident events
     // double delta_p_percent = 2./100.;
     // double delta_p = p_peak*delta_p_percent;

     // final cuts
     TCut ACC = DP+isPb+colCut+upCut+dnCut+XCUT;

     // load acceptance table
     int nbins=100;
     double accp_angle[nbins], accp[nbins], accp_err[nbins];
     for (int i=0;i<nbins;i++) { accp_angle[i]=0; accp[i]=0; accp_err[i]=0;}
     TString acc_shape = "original";  // option: original, smear, box, gaus, shift

     int status = LoadACC(accfun.Data(),accp_angle, accp);
     if(status==0) exit(0);  // asymmetry table doesn't exist

     Double_t mean0=0, max0=0, rms0=0;
     GetACCshape(accfun.Data(),mean0,rms0,max0);

     // gaus shape
     TF1 *f_gaus = new TF1("f_gaus","gaus",3,8);
     f_gaus->SetParameters(max0, mean0, rms0);

     // box shape
     double th_1 = mean0 - sqrt(3.)*rms0;
     double th_2 = mean0 + sqrt(3.)*rms0;
     TF1 *f_box = new TF1("f_box",Form("%f*((x>=%f && x<=%f)? 1 : 0)",max0,th_1,th_2),3,8);

     // smearing: original + rms_percent
     double rms_percent = 3./100.;
     double_t accp_angle_p[nbins];
     for(int ii=0; ii<nbins; ii++){
	accp_angle_p[ii] = mean0+(accp_angle[ii]-mean0)*(1.+rms_percent); // +rms_percent
     }

     // shift
     double shift_th = 0.5; // deg
     double_t accp_angle_sf[nbins];
     for(int ii=0; ii<nbins; ii++){
	accp_angle_sf[ii] = accp_angle[ii]+shift_th; // +rms_percent
     }

     bool check_ACC = true; // check mean and rms of acceptance function
     TH1F *hACC_check = new TH1F("hACC_check","acceptance function distribution",nbins,3,8);
     if(check_ACC){
        int check_nbin_th = 100000;
        double dtheta = (8.0-3.0)/(check_nbin_th*1.0);  // delta theta in radius
        for(int ii=0; ii<check_nbin_th; ii++){
          double thisAngle = 3.0 + dtheta*ii;
          double thisACC = 0;
	  if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,nbins);  // original acceptance function
	  if(acc_shape == "box" ) thisACC = f_box->Eval(thisAngle);	
	  if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisAngle);	
	  if(acc_shape == "smear" ) thisACC = FindACC(thisAngle,accp_angle_p,accp,accp_err,nbins);  
	  if(acc_shape == "shift" ) thisACC = FindACC(thisAngle,accp_angle_sf,accp,accp_err,nbins);

          if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

          hACC_check->Fill(thisAngle,thisACC/1000.);
        }

        Double_t mean = hACC_check->GetMean();
        Double_t RMS = hACC_check->GetRMS();
        Double_t mean_bin = hACC_check->FindBin(mean);
        Double_t peak = hACC_check->GetBinContent(mean_bin);
        if(acc_shape=="smear") cout<<Form("smear %.3f",rms_percent)<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
        else {
	  if(acc_shape=="shift") cout<<Form("shift %.3f",shift_th)<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
          else cout<<acc_shape<<" "<<"ACC function peak, mean, rms"<<"  "<<peak<<"  "<<mean<<"  "<<RMS<<endl;
	}
     }

     //     return;    
     // load asymmetry table
     LoadTable("../../model/horpb.dat", 0);
     LoadTable("../../model/horpb1.dat", 1);

     int nabins=250;
     TH1F *hasym_origin = new TH1F("hasym_origin","vertex asymmetry in g4hrs with ACC cut",nabins,200,1000);
     TH1F *hasym_acc1 = new TH1F("hasym_acc1","vertex asymmetry calculated from acceptance function with dp cut",nabins,200,1000);
     TH1F *hasym_e_acc = new TH1F("hasym_e_acc","elemental asymmetry with acceptance function applied",nabins,200,1000);
     TH1F *hasyms_e_acc = new TH1F("hasyms_e_acc","elemental asymmetry with acceptance function applied",nabins,200,1000);

     TH1F *hth_origin = new TH1F("hth_origin","vertex theta in g4hrs with ACC cut",150,3,8);
     TH1F *hth_e_acc = new TH1F("hth_e_acc","elemental theta with acceptance function applied",150,3,8);

     TH1F *hq2_origin = new TH1F("hq2_origin","vertex Q2 in g4hrs with ACC cut",150,0,0.015);
     TH1F *hq2_e_acc = new TH1F("hq2_e_acc","elemental Q2 with acceptance function applied",150,0,0.015);

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
	double thisACC = 0;
	if(acc_shape == "original" ) thisACC = FindACC(thisTh,accp_angle,accp,accp_err,nbins);  // original acceptance function
	if(acc_shape == "box" ) thisACC = f_box->Eval(thisTh);	
	if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisTh);	
	if(acc_shape == "smear" ) thisACC = FindACC(thisTh,accp_angle_p,accp,accp_err,nbins); 
	if(acc_shape == "shift" ) thisACC = FindACC(thisTh,accp_angle_sf,accp,accp_err,nbins); 

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	//if( thisA==208 && (p_peak-p_ztarg)<delta_p) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
	if( thisA==208 && abs(thisBeamE*1000.-p_peak)<2.2) hasym_acc1->Fill(thisAsym*ppb,thisRate*thisACC);
     }

     double setBeamE = E0/1000.; 
     int nbin_th = 10000;
     double dtheta = (8.0-3.0)/(nbin_th*1.0);  // delta theta in radius

     for(int ii=0; ii<nbin_th; ii++){
	double thisAngle = 3.0 + dtheta*ii;
	double thisAsym = Interpolate(setBeamE,thisAngle,0,1);  // elemental asymmetry
	double thisXS = Interpolate(setBeamE,thisAngle,0,0);  // elemental cross section
	double thisAsymS = Interpolate(setBeamE,thisAngle,1,1);
	//	double thisAsymS = 0;
	
	double deg_to_rad = TMath::Pi()/180.;
	double thisEp = setBeamE/(1.+setBeamE/(208.*amu_c2)*(1.-TMath::Cos(thisAngle*deg_to_rad)));
	double thisQ2 = 4.*setBeamE*thisEp*pow(TMath::Sin(thisAngle*deg_to_rad/2.),2);

	double thisACC=0;
	if(acc_shape == "original" ) thisACC = FindACC(thisAngle,accp_angle,accp,accp_err,nbins);  // original acceptance function
	if(acc_shape == "box" ) thisACC = f_box->Eval(thisAngle);	
	if(acc_shape == "gaus" ) thisACC = f_gaus->Eval(thisAngle);	
	if(acc_shape == "smear" ) thisACC = FindACC(thisAngle,accp_angle_p,accp,accp_err,nbins);  
	if(acc_shape == "shift" ) thisACC = FindACC(thisAngle,accp_angle_sf,accp,accp_err,nbins);  

	if(thisACC<0) {printf("Something wrong here: ACC= %f\n",thisACC); exit(0);}

	double omega = TMath::Sin(thisAngle*deg_to_rad)*dtheta*deg_to_rad;
	hasym_e_acc->Fill(thisAsym*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
	hasyms_e_acc->Fill(thisAsymS*ppb,thisXS*omega*thisACC);   // elemetal asymmetry with the acceptance function applied
	hth_e_acc->Fill(thisAngle,thisXS*omega*thisACC);   // elemetal angle with the acceptance function applied
	hq2_e_acc->Fill(thisQ2,thisXS*omega*thisACC);   // elemetal Q2 with the acceptance function applied
     } 

     cout<<"================================="<<endl;
     cout<<"orignal vertex asymmetry mean:     "<<setw(8)<<hasym_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_origin->GetRMS()<<endl;
     cout<<"reproduce vertex asymmetry mean:   "<<setw(8)<<hasym_acc1->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_acc1->GetRMS()<<endl;
     cout<<"elemental asymmetry mean:          "<<setw(8)<<hasym_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hasym_e_acc->GetRMS()<<endl;
     cout<<"elemental stretched asymmetry mean:          "<<setw(8)<<hasyms_e_acc->GetMean()<<"   "<<"Err:  "<<setw(8)<<hasyms_e_acc->GetMeanError()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex angle mean:   "<<setw(8)<<hth_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_origin->GetRMS()<<endl;
     cout<<"elemental angle mean:        "<<setw(8)<<hth_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hth_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;
     cout<<"orignal vertex Q2 mean:   "<<setw(8)<<hq2_origin->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_origin->GetRMS()<<endl;
     cout<<"elemental Q2 mean:        "<<setw(8)<<hq2_e_acc->GetMean()<<"   "<<"RMS:  "<<setw(8)<<hq2_e_acc->GetRMS()<<endl;
     cout<<"================================="<<endl;

     Double_t inte1 = hasym_origin->Integral();
     Double_t inte2 = hasym_e_acc->Integral();
     Double_t inte3 = hasym_acc1->Integral();
     hasym_e_acc->Scale(inte1/inte2);
     //hasym_acc1->Scale(inte1/inte3);

     TCanvas *c1 = new TCanvas("c1","c1",800,600);
     hasym_origin->SetLineColor(1);
     hasym_acc1->SetLineColor(6);
     hasym_e_acc->SetLineColor(2);
     hasyms_e_acc->SetLineColor(3);

     hasym_origin->SetLineWidth(2);
     hasym_acc1->SetLineWidth(2);
     hasym_e_acc->SetLineWidth(2);
     hasyms_e_acc->SetLineWidth(3);

     if(inte1>inte3){
       hasym_origin->Draw("HIST");    
       hasym_acc1->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_origin->SetTitle("LHRS Asymmetry with acceptance; asym(ppb)");
     }
     else{
       hasym_acc1->Draw("HIST");    
       hasym_origin->Draw("HIST sames");    
       hasym_e_acc->Draw("HIST sames");
       hasym_acc1->SetTitle("LHRS Asymmetry with acceptance; asym(ppb)");
     }

     TLegend *leg = new TLegend(0.15,0.65,0.35,0.8);
     leg->AddEntry(hasym_origin,"origin","L");
     //leg->AddEntry(hasym_acc1,Form("reproduce (dp_inc/p=%.3f)",delta_p_percent),"L");
     leg->AddEntry(hasym_acc1,Form("reproduce (beamE cut)"),"L");
     leg->AddEntry(hasym_e_acc,"elemental","L");
     leg->Draw();

   gPad->Update();
   TPaveStats* statA1 = (TPaveStats*)hasym_origin->FindObject("stats");
   TPaveStats* statA2 = (TPaveStats*)hasym_acc1->FindObject("stats");
   TPaveStats* statA3 = (TPaveStats*)hasym_e_acc->FindObject("stats");
   statA1->SetY1NDC(0.90);
   statA1->SetY2NDC(0.75);
   statA2->SetY1NDC(0.75);
   statA2->SetY2NDC(0.60);
   statA3->SetY1NDC(0.60);
   statA3->SetY2NDC(0.45);
   statA1->SetTextColor(1);
   statA2->SetTextColor(6);
   statA3->SetTextColor(2);
   gPad->Modified();

     TLatex tex1;
     tex1.SetTextAlign(11);
     tex1.SetTextSize(0.02);
     if(acc_shape=="smear") tex1.DrawLatexNDC(0.65,0.3,Form("smear %.3f",rms_percent));
     else {
	if(acc_shape=="shift") tex1.DrawLatexNDC(0.65,0.3,Form("shift %.3f",shift_th));
        else tex1.DrawLatexNDC(0.65,0.3,Form("%s acceptance",acc_shape.Data()));
     }
    
     inte1 = hth_origin->Integral();
     inte2 = hth_e_acc->Integral();
     hth_e_acc->Scale(inte1/inte2);

     TCanvas *c2 = new TCanvas("c2","c2",800,600);
     c2->Divide(2,1);
     c2->cd(1);
     hth_origin->SetLineColor(1);
     hth_e_acc->SetLineColor(2);

     hth_origin->SetLineWidth(2);
     hth_e_acc->SetLineWidth(2);

     Double_t hpeak1 = hth_origin->GetMaximum();
     Double_t hpeak2 = hth_e_acc->GetMaximum();
  
     if(hpeak2<hpeak1){
       hth_origin->Draw("HIST");    
       hth_e_acc->Draw("HIST sames");
       hth_origin->SetTitle("LHRS angle with acceptance; angle");
     }
     else{
       hth_e_acc->Draw("HIST");
       hth_origin->Draw("HIST sames");    
       hth_e_acc->SetTitle("LHRS angle with acceptance; angle");
     }


     TLegend *leg1 = new TLegend(0.5,0.65,0.7,0.8);
     leg1->AddEntry(hth_origin,"origin","L");
     leg1->AddEntry(hth_e_acc,"elemental","L");
     leg1->Draw();

   gPad->Update();
   TPaveStats* statB1 = (TPaveStats*)hth_origin->FindObject("stats");
   TPaveStats* statB2 = (TPaveStats*)hth_e_acc->FindObject("stats");
   statB1->SetY1NDC(0.90);
   statB1->SetY2NDC(0.75);
   statB2->SetY1NDC(0.75);
   statB2->SetY2NDC(0.60);
   statB1->SetTextColor(1);
   statB2->SetTextColor(6);
   gPad->Modified();
    
     c2->cd(2);
     inte1 = hq2_origin->Integral();
     inte2 = hq2_e_acc->Integral();
     hq2_e_acc->Scale(inte1/inte2);

     hq2_origin->SetLineColor(1);
     hq2_e_acc->SetLineColor(2);

     hq2_origin->SetLineWidth(2);
     hq2_e_acc->SetLineWidth(2);

     hpeak1 = hq2_origin->GetMaximum();
     hpeak2 = hq2_e_acc->GetMaximum();

     if(hpeak1>hpeak2){
       hq2_origin->Draw("HIST");    
       hq2_e_acc->Draw("HIST sames");
       hq2_origin->SetTitle("LHRS Q2 with acceptance; Q2");
     }
     else{
       hq2_e_acc->Draw("HIST");
       hq2_origin->Draw("HIST sames");    
       hq2_e_acc->SetTitle("LHRS Q2 with acceptance; Q2");
     }


     TLegend *leg2 = new TLegend(0.15,0.65,0.35,0.8);
     leg2->AddEntry(hq2_origin,"origin","L");
     leg2->AddEntry(hq2_e_acc,"elemental","L");
     leg2->Draw();

   gPad->Update();
   TPaveStats* statC1 = (TPaveStats*)hq2_origin->FindObject("stats");
   TPaveStats* statC2 = (TPaveStats*)hq2_e_acc->FindObject("stats");
   statC1->SetY1NDC(0.90);
   statC1->SetY2NDC(0.75);
   statC2->SetY1NDC(0.75);
   statC2->SetY2NDC(0.60);
   statC1->SetTextColor(1);
   statC2->SetTextColor(6);
   gPad->Modified();
    
     TLatex tex2;
     tex2.SetTextAlign(11);
     tex2.SetTextSize(0.03);
     if(acc_shape=="smear") tex2.DrawLatexNDC(0.65,0.3,Form("smear %.3f",rms_percent));
     else {
	if(acc_shape=="shift") tex2.DrawLatexNDC(0.65,0.3,Form("shift %.3f",shift_th));
        else tex2.DrawLatexNDC(0.5,0.5,Form("%s acceptance",acc_shape.Data()));
     }
   
     TCanvas *c3;
     if(check_ACC){
       c3 = new TCanvas("c3","c3",800,600);
       hACC_check->Draw();
     }



     
     if (iprint ==1) {
       c1->Print("fig/Accfun_"+stub+"_overlay.pdf");
       c2->Print("fig/Accfun_"+stub+"_overlay2.pdf");
       c3->Print("fig/Accfun_"+stub+"_summary.pdf");
     }
}
