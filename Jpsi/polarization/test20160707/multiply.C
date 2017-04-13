void multiply(){
	TFile* file1 = new TFile("efficiency.root","read");
	TFile* file2 = new TFile("histograms.root","read");

	TCanvas* c = new TCanvas();
	c->Divide(2,2);

	TH2F* h1 = (TH2F*)(file1->Get("eff_ratio_0_0_2_0"));
	TH2F* h2 = (TH2F*)(file2->Get("theta_20_phi_20"));
	TH2F* h3 = (TH2F*)(file1->Get("eff_ratio_0_0_2_0"));	
	c->cd(1);
	h1->Draw("colz");
	c->cd(2);
	h2->Draw("colz");

	c->cd(3);
//	h3->Multiply(h2);
//	h2->Multiply(h1);
	h2->Draw("colz");







}
