#include "ReadFile.h"

void ACC_avg_asym(){
     
     Double_t angle_asym[240]={0.0}, asym[240]={0.0}, xs[240]={0.0};
     Double_t angle_acc[100]={0.0}, accp[100]={0.0};

     TString file1 = "elastic.dat"; 
     ReadAsym(file1,angle_asym,asym,xs);

     //TString file2 = "/w/halla-scifs17exp/parity/disk1/hanjie/ACC/Kent_code/accfunction_norm_1000.csv"; 
     TString file2 = "/w/halla-scifs17exp/parity/disk1/hanjie/ACC/Kent_code/accfunc_central_average.csv"; 
     ReadACC_Kent(file2,angle_acc,accp);

     Double_t weight=0.0;
     Double_t asym_w = 0.0;

     for(int ii=0; ii<240; ii++){
	Double_t tmp_accp=1.0;
	int found=0;
	if(angle_asym[ii]>angle_acc[99])continue;

	for(int jj=0; jj<100; jj++){
	  if( angle_asym[ii]<=angle_acc[jj] ){
	      //if(jj==0)tmp_accp=accp[jj];
	      if(jj==0)tmp_accp=0;
	      else{
		tmp_accp = ( accp[jj]-accp[jj-1] )/( angle_acc[jj]-angle_acc[jj-1] )*( angle_asym[ii]-angle_acc[jj-1] )+accp[jj-1];
//cout<<angle_acc[jj]<<"  "<<angle_asym[ii]<<"  "<<angle_acc[jj-1]<<"  "<<tmp_accp<<endl;
	      }
	      found=1 ;
	      break;
	  }
	}	

	if(found==1){
	   weight = weight + tmp_accp*xs[ii]*sin( angle_asym[ii]*TMath::Pi()/180.0);
	   asym_w = asym_w + asym[ii]*tmp_accp*xs[ii]*sin( angle_asym[ii]*TMath::Pi()/180.0);
//cout<<tmp_accp<<"  "<<asym[ii]<<"  "<<weight<<"  "<<asym_w<<endl;
	}

     }

     Double_t asym_final = asym_w/weight;
     cout<<"acceptance averaged asymmetry:  "<<endl;
     cout<<asym_final<<endl;
  
}
