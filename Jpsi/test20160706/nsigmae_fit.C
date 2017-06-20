void nsigmae_fit(){
	TFile* infile = new TFile("Nsigma.root","read");
	TH1F *nsigma_mean_0, *nsgima_mean_plus, *nsigma_mean_minus;
	TH1F *nsigma_sigma_0, *nsigma_sigma_plus, *nsigma_sigma_minus;

	nsigma_mean_0 = (TH1F*)infile->Get("mh1mean");
	nsigma_mean_0->SetName("nsigma_mean_0");
	nsigma_mean_plus = (TH1F*)nsigma_mean_0->Clone("nsigma_mean_plus");
	nsigma_mean_plus->SetName("nsigma_mean_plus");
	nsigma_mean_minus = (TH1F*)nsigma_mean_0->Clone("nsigma_mean_minus");
	nsigma_mean_minus->SetName("nsigma_mean_minus");

	nsigma_sigma_0 = (TH1F*)infile->Get("mh1sigma");
	nsigma_sigma_0->SetName("nsigma_sigma_0");
	nsigma_sigma_plus = (TH1F*)nsigma_sigma_0->Clone("nsigma_mean_plus");
	nsigma_sigma_plus->SetName("nsigma_sigma_plus");
	nsigma_sigma_minus = (TH1F*)nsigma_sigma_0->Clone("nsigma_sigma_minus");
	nsigma_sigma_minus->SetName("nsigma_sigma_minus");	

	for(int i=0;i<nsigma_mean_0->GetSize()-2;i++){
		nsigma_mean_plus->SetBinContent(i+1,nsigma_mean_0->GetBinContent(i+1)+nsigma_mean_0->GetBinError(i+1));
		nsigma_mean_plus->SetBinError(i+1,nsigma_mean_0->GetBinError(i+1));
		nsigma_mean_minus->SetBinContent(i+1,nsigma_mean_0->GetBinContent(i+1)-nsigma_mean_0->GetBinError(i+1));
		nsigma_mean_minus->SetBinError(i+1,nsigma_mean_0->GetBinError(i+1));
	
		nsigma_sigma_plus->SetBinContent(i+1,nsigma_sigma_0->GetBinContent(i+1)+nsigma_sigma_0->GetBinError(i+1));
		nsigma_sigma_plus->SetBinError(i+1,nsigma_sigma_0->GetBinError(i+1));
		nsigma_sigma_minus->SetBinContent(i+1,nsigma_sigma_0->GetBinContent(i+1)-nsigma_sigma_0->GetBinError(i+1));
		nsigma_sigma_minus->SetBinError(i+1,nsigma_sigma_0->GetBinError(i+1));
	}

	TF1 *mean1, *mean2, *mean3, *mean4;
	mean1 = new TF1("mean1","[0]");
	mean2 = new TF1("mean2","[0]");
	mean3 = new TF1("mean3","[0]");
	mean4 = new TF1("mean4","[0]+[1]*x");
	
	TF1 *sigma1, *sigma2, *sigma3, *sigma4;
	sigma1 = new TF1("sigma1","[0]");
	sigma2 = new TF1("sigma2","[0]");
	sigma3 = new TF1("sigma3","[0]");
	sigma4 = new TF1("sigma4","[0]+[1]*x");

	TCanvas *c = new TCanvas("c","c",2000,1000);
	c->Divide(4,2);
	c->cd(1);
	nsigma_mean_0->Fit("mean1");
	nsigma_mean_0->GetYaxis()->SetRangeUser(-2.,0.5);
	nsigma_mean_0->SetTitle(";p_{T} GeV/c; #mu");
	nsigma_mean_0->Draw();
	c->cd(2);
	nsigma_mean_plus->Fit("mean2");
	nsigma_mean_plus->GetYaxis()->SetRangeUser(-2.,0.5);
	nsigma_mean_plus->SetTitle(";p_{T} GeV/c; #mu");
	nsigma_mean_plus->Draw();
	c->cd(3);
	nsigma_mean_minus->Fit("mean3");
	nsigma_mean_minus->GetYaxis()->SetRangeUser(-2.,0.5);
	nsigma_mean_minus->SetTitle(";p_{T} GeV/c; #mu");
	nsigma_mean_minus->Draw();
	c->cd(4);
	TH1F *nsigma_mean_linear = (TH1F*)nsigma_mean_0->Clone("nsigma_mean_linear");
	nsigma_mean_linear->Fit("mean4");
	nsigma_mean_linear->GetYaxis()->SetRangeUser(-2.,0.5);
	nsigma_mean_linear->SetTitle(";p_{T} GeV/c; #mu");
	nsigma_mean_linear->Draw();
	c->cd(5);
	nsigma_sigma_0->Fit("sigma1");
	nsigma_sigma_0->GetYaxis()->SetRangeUser(0.1,4.);
	nsigma_sigma_0->SetTitle(";p_{T} GeV/c; #sigma");
	nsigma_sigma_0->Draw();
	c->cd(6);
	nsigma_sigma_plus->Fit("sigma2");
	nsigma_sigma_plus->GetYaxis()->SetRangeUser(0.1,4.);
	nsigma_sigma_plus->SetTitle(";p_{T} GeV/c; #sigma");
	nsigma_sigma_plus->Draw();
	c->cd(7);
	nsigma_sigma_minus->Fit("sigma3");
	nsigma_sigma_minus->GetYaxis()->SetRangeUser(0.1,4.);
	nsigma_sigma_minus->SetTitle(";p_{T} GeV/c; #sigma");
	nsigma_sigma_minus->Draw();
	c->cd(8);
	TH1F *nsigma_sigma_linear = (TH1F*)nsigma_sigma_0->Clone("nsigma_sigma_linear");
	nsigma_sigma_linear->Fit("sigma4");
	nsigma_sigma_linear->GetYaxis()->SetRangeUser(0.1,4.);
	nsigma_sigma_linear->SetTitle(";p_{T} GeV/c; #sigma");
	nsigma_sigma_linear->Draw();

	TFile* outfile = new TFile("nsigmae_fit.root","recreate");
	outfile->cd();
	nsigma_mean_0->Write();
	nsigma_mean_plus->Write();
	nsigma_mean_minus->Write();
	nsigma_sigma_0->Write();
	nsigma_sigma_plus->Write();
	nsigma_sigma_minus->Write();
	
	mean1->Wrtie();
	mean2->Write();
	mean3->Write();
	mean4->Write();

	sigma1->Write();
	sigma2->Write();
	sigma3->Write();
	sigma4->Write();
	
	c->Write();
	c->SaveAs("nsigmae_fit.pdf");
	outfile->Close();

//	nsigma_mean_plus = (TH1F*)infile->Get();
//	nsigma_mean_minus = (TH1F*)infile->Get();

//	nsigma_sigma_plus = (TH1F*)infile->Get();
//	nsigma_sigma_minus = (TH1F*)infile->Get();

//	nsigma_mean_0->Fit("mean1");
//	nsigma_mean_plus->Fit("mean2");
//	nsigma_mean_minus->Fit("mean3");

//	nsigma_sigma_0->Fit("mean2");

}
