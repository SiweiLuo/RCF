TString trigSet[4] = {"MB","HT0","HT1","HT2"};
TString ptbin[6] = {"0 < p_{T} < 2 GeV/c","2 < p_{T} < 3 GeV/c","3 < p_{T} < 4 GeV/c","4 < p_{T} < 6 GeV/c","6 < p_{T} < 8 GeV/c","8 < p_{T} < 12 GeV/c"};

void PWG_compare_invariant(int eid=2, int inm=30){
	gStyle->SetOptStat(false);
	gStyle->SetOptFit(false);
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);

	TLegend *leg[6];

	TFile* infile1 = new TFile(Form("~/polarization/invmass_%deid_inv%d/rootfile20170912/data_invariant_20171030.root",eid,inm),"read");
	TFile* infile2 = new TFile(Form("rootfile20170921/emb_invariant_inv%d.root",inm),"read");
	TFile* outfile = new TFile("rootfile20170921/PWG_compare_invariant.root","recreate");

	outfile->cd();
	TH1F* inv[6][2];	

	TCanvas* c[4];

	TCanvas* invc;
	invc = new TCanvas("invc","invc",1200,1200);
	invc->Divide(2,2);

	TCanvas* PWGinvc;
	PWGinvc = new TCanvas("PWGinvc","PWGinvc",1200,360);
	PWGinvc->Divide(3,1);

	TCanvas* ptc;
	ptc = new TCanvas("ptc","ptc",1200,1200);
	ptc->Divide(2,2);

	TLegend* invleg[4];
	TLegend* ptleg[4];

	TH2F* data2d[4];
	TH2F* emb2d[4];	

	TH1F* datainvm[4];
	TH1F* embinvm[4];
	TH1F* datapt[4];
	TH1F* embpt[4];	

	TF1 *fCB2 = new TF1("fCB2", CrystalBall2_POL, 2.5, 3.5, 9);
	fCB2->SetLineWidth(1.);
	fCB2->SetLineColor(kRed);
	fCB2->SetParameters(1.e5, 3.005, 0.09, 1.5, 1.5, 1.5, 1.5 , 1.5 , 1.5 );
	fCB2->SetParLimits(1, 2., 4.);
	fCB2->SetParLimits(2, 0., 1.2);
	fCB2->SetParLimits(3, 0., 20.);
	fCB2->SetParLimits(4, 0., 20.);
	fCB2->SetParLimits(5, 0., 20.);
	fCB2->SetParLimits(6, 0., 20.);

	TLatex pttxt[6];

	TCanvas *canvas = new TCanvas("canvas","canvas",1600,400);
	canvas->Divide(6,1,0,0);

	Int_t xmin,xmax;

	for(int pt=0,trig=0;pt<6;pt++){
		canvas->cd(1+pt);
		gPad->SetTickx(2);
		if(pt==0) trig=0;
		else if(pt>=1 && pt<=2) trig=1;
		else if(pt>=3 && pt<=5) trig=3;
		data2d[trig] = (TH2F*)infile1->Get(Form("InvMpT%d",trig));	
		emb2d[trig] = (TH2F*)infile2->Get(Form("h%sJpsiPtInvM",trigSet[trig].Data()));

		datainvm[trig] = (TH1F*)data2d[trig]->ProjectionY(Form("datainvm%d",trig-1));
		if(trig==0) datainvm[trig]->SetTitle(Form("MB Invariant mass comparison;m_{ee} GeV/c^{2};"));	
		else datainvm[trig]->SetTitle(Form("HT%d Invariant mass comparison;m_{ee} GeV/c^{2};",trig-1));
		cout<<"integral from 3.0 to 3.15 ==="<<datainvm[trig]->Integral(datainvm[trig]->GetXaxis()->FindBin(3.0),datainvm[trig]->GetXaxis()->FindBin(3.15))<<endl;
		cout<<"datainvm integral ====="<<datainvm[trig]->Integral(21,24)<<endl;
		embinvm[trig] = (TH1F*)emb2d[trig]->ProjectionY(Form("embinvm%d",trig-1));
		embinvm[trig]->SetLineColor(kRed);

		xmin = datainvm[trig]->GetXaxis()->FindBin(3.0);
		xmax = datainvm[trig]->GetXaxis()->FindBin(3.15);

		inv[pt][0] = (TH1F*)infile1->Get(Form("pro_trig%d_pT%d",trig,pt));
		inv[pt][0]->SetName(Form("pro_trig%d_pT%d",trig,pt));
		inv[pt][1] = (TH1F*)infile2->Get(Form("trig%d_pt%d",trig,pt));
		inv[pt][1]->SetName(Form("trig%d_pt%d",trig,pt));
		inv[pt][1]->Scale(inv[pt][0]->Integral(xmin,xmax-1)/inv[pt][1]->Integral(xmin,xmax-1));
		inv[pt][0]->SetMaximum(inv[pt][1]->GetMaximum()*1.5);
		inv[pt][0]->SetTitle(Form("Invariant mass spectra;m_{ee} GeV/c^{2}"));

		//		fCB2->SetParLimits(7, -10., 10.);
		//		fCB2->SetParLimits(8, -10., 10.);
		inv[pt][0]->GetYaxis()->SetRangeUser(-10,260);
		inv[pt][0]->GetXaxis()->SetLimits(2.5,3.5);
		inv[pt][1]->GetYaxis()->SetRangeUser(-10,260);
		inv[pt][1]->GetXaxis()->SetLimits(2.5,3.5);

		inv[pt][0]->Draw(); // data
		//		inv[pt][1]->Draw("HIST same"); // efficiency
		cout<<" bins ==== "<<inv[pt][0]->GetNbinsX()<<"   "<<inv[pt][1]->GetNbinsX()<<endl;


		inv[pt][0]->Fit("fCB2");

		pttxt[pt].SetNDC(1);
		if(pt==0){
			pttxt[pt].DrawLatex(0.15,0.9,Form("#J/#psi signal S = %.1f",inv[pt][0]->Integral(datainvm[trig]->GetXaxis()->FindBin(3.0),datainvm[trig]->GetXaxis()->FindBin(3.15))));
			//		pttxt[pt].DrawLatex(0.15,0.8,Form("#J/#psi background B = %.1f",inv[pt][0]->Integral(datainvm[trig]->GetXaxis()->FindBin(3.0),datainvm[trig]->GetXaxis()->FindBin(3.15))));
			//		pttxt[pt].DrawLatex(0.15,0.7,Form("#J/#psi S/#sqrt{S+2B} = %.2f",inv[pt][0]->Integral(datainvm[trig]->GetXaxis()->FindBin(3.0),datainvm[trig]->GetXaxis()->FindBin(3.15))));
		}
		else {
			pttxt[pt].DrawLatex(0.07,0.9,Form("#J/#psi signal S = %.1f",inv[pt][0]->Integral(datainvm[trig]->GetXaxis()->FindBin(3.0),datainvm[trig]->GetXaxis()->FindBin(3.15))));
		}


		leg[pt] = new TLegend(0.56,0.7,0.95,0.95);
		leg[pt]->SetBorderSize(0);
		leg[pt]->SetFillStyle(0);
		leg[pt]->SetTextSize(0.04);
		leg[pt]->AddEntry(inv[pt][0],"data","lep");
		leg[pt]->AddEntry(fCB2,"crystalball+POL1","l");
		leg[pt]->Draw("same");

	}
	canvas->SaveAs("~/WWW/Proposal/figure/PWGcanvas.pdf");
}

Double_t CrystalBall2_POL(Double_t *x, Double_t *par)
{
	Double_t N = par[0];
	Double_t mu = par[1];
	Double_t s = par[2];
	Double_t n1 = par[3];
	Double_t alpha1 = par[4];
	Double_t n2 = par[5];
	Double_t alpha2 = par[6];
	Double_t linea = par[7];
	Double_t lineb = par[8];

	Double_t A = TMath::Power(n1/fabs(alpha1), n1) * TMath::Exp(-alpha1*alpha1/2.);
	Double_t B = n1/fabs(alpha1) - fabs(alpha1);

	Double_t C = TMath::Power(n2/fabs(alpha2), n2) * TMath::Exp(-alpha2*alpha2/2.);
	Double_t D = n2/fabs(alpha2) - fabs(alpha2);

	Double_t norm = (x[0]-mu)/s;

	Double_t result=0;
	if(norm < -alpha1) {
		result = N * A * TMath::Power(B-norm, -n1);
	} else if(norm < alpha2) {
		result = N * TMath::Exp(-0.5*norm*norm);
	} else {
		result = N * C * TMath::Power(D+norm, -n2);
	}
	result += par[7]*x[0]+par[8];
	return result;
}

