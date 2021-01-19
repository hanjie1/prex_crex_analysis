void test(){
	TH1F *h1 = new TH1F("h1","h1",10,0,10);
	TH1F *h2 = new TH1F("h2","h2",10,0,10);
	TH1F *h3 = new TH1F("h3","h3",10,0,10);
	Double_t mean1 = 0.0;
	Double_t meanerr1 = 0.0;
	Double_t i1 = 0.0;
	Double_t mean2 = 0.0;
	Double_t meanerr2 = 0.0;
	Double_t i2 = 0.0;
	Double_t mean3 = 0.0;
	Double_t meanerr3 = 0.0;
	Double_t i3 = 0.0;
 	for(int ii=0; ii<5;ii++) {
	   h1->Fill(ii);
	   mean1 = mean1 + ii;
	   i1 = i1+1;
	   mean3 = mean3 + ii;
	   i3 = i3+1;
	}
	cout<<"1 mean:  "<<h1->GetMean()<<endl;
 	for(int ii=5; ii<10; ii++) {
	   h1->Fill(ii,0.5);
	   mean2 = mean2 + ii*0.5;
	   i2 = i2+(1*0.5);
	   mean3 = mean3 + ii*0.5;
	   i3 = i3+(1*0.5);
	}
 	for(int ii=0; ii<5;ii++) {
	   meanerr1 = meanerr1 + (pow(ii-mean1,2)/pow(1.0*i1,2));
	   meanerr3 = meanerr3 + (pow(ii-mean3,2)/pow(1.0*i3,1));
	}
 	for(int ii=5; ii<10;ii++) {
	   meanerr2 = meanerr2 + (pow(ii-mean2,2)/pow(i2,2));
	   meanerr3 = meanerr3 + 0.5*(pow(ii-mean3,2)/pow(i3,1));
	}
	std::cout << mean1 << " summed, average = " << 1.0*mean1/(1.0*i1) << " +- " << sqrt(meanerr1) << std::endl;
	std::cout << mean2 << " summed, average = " << 1.0*mean2/(1.0*i2) << " +- " << sqrt(meanerr2) << std::endl;
	std::cout << mean3 << " summed, average = " << 1.0*mean3/(1.0*i3) << " +- " << sqrt(meanerr3 *10./9.) << std::endl;

	h1->Draw();
	cout<<"mean:  "<<h1->GetMean()<<endl;
	cout<<"mean err:  "<<h1->GetMeanError()<<endl;

}
