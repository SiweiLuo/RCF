
void invMass(){
	TFile* inputfile = new TFile("rootfile/ht1_trg4_1TrkPid_0_0.ana.root","read"); 

	TH1F* unlike = (TH1F*)((TH3F*)inputfile->Get("hJpsiCosThetaInvMPt"))->ProjectionY("h1");
	TH1F* like = (TH1F*)((TH3F*)inputfile->Get("hJpsiCosThetaInvMPtBG"))->ProjectionY("h2");

	signal = (TH1F*)unlike->Clone();
	signal->Add(like,-1);

	signal->Draw();

	int min = signal->GetXaxis()->FindBin(2.9);
	int max = signal->GetXaxis()->FindBin(3.2)-1;

	double Jpsi = signal->Integral(min,max);

	cout<<"=============>"<<Jpsi<<endl;

}
