void PlotVZP(){
	TChain *T = new TChain("T");
	T->Add("g4hrs_rootfiles/prex_c12_wIon.root");

	TH2F *hEz = new TH2F("hEz","c12 beamp vs. vertex z",100,-0.0015,0.0015,100,948,950);
	TH2F *hptgz = new TH2F("hptgz","c12 p_tg vs. vertex z",100,-0.0015,0.0015,100,948,950);
	TH2F *hptargz = new TH2F("hptargz","c12 p_ztarg vs. vertex z",100,-0.0015,0.0015,100,948,950);

	TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
	c1->Divide(2,2);
	c1->cd(1);
	T->Draw("ev.beamp*1000:ev.vz>>hEz","","COLZ");	
	c1->cd(2);
	T->Draw("p_tg:ev.vz>>hptgz","","COLZ");	
	c1->cd(3);
	T->Draw("p_ztarg:ev.vz>>hptargz","","COLZ");	

}
