
TString ptbin[6] = {"0 < p_{T} < 2 GeV/c","2 < p_{T} < 3 GeV/c","3 < p_{T} < 4 GeV/c","4 < p_{T} < 6 GeV/c","6 < p_{T} < 8 GeV/c","8 < p_{T} < 12 GeV/c"};

void compare_invariant(int eid=2, int inm=30){
	//	int eid=1,inm =30;
	gStyle->SetOptStat(false);

	TFile* infile1 = new TFile(Form("~/polarization/invmass_%deid_inv%d/rootfile20170912/data_invariant_20170912.root",eid,inm),"read");
//	TFile* infile2 = new TFile(Form("rootfile20170920/emb_invariant_inv%d.root",inm),"read");
	TFile* infile2 = new TFile(Form("rootfile20170921/emb_invariant_inv%d.root",inm),"read");
//	TFile* outfile = new TFile("rootfile20170920/compare_invariant.root","recreate");
	TFile* outfile = new TFile("rootfile20170921/compare_invariant.root","recreate");

	outfile->cd();
	TH1F* inv[5][6][2];	

	TCanvas* c[4];
	TLegend* leg[4];

	TCanvas* invc;
	invc = new TCanvas("invc","invc",1200,1200);
	invc->Divide(2,2);
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

	TLatex pttxt;

	for(int trig=0;trig<3;trig++){

		data2d[trig] = (TH2F*)infile1->Get(Form("InvMpT%d",trig));	
		data2d[trig]->Write();
	
		if(eid==1){
			emb2d[trig] = (TH2F*)infile2->Get(Form("hHT%dJpsiPtInvM",trig));
			emb2d[trig]->Write();
		}
		else if(eid==2){
			emb2d[trig] = (TH2F*)infile2->Get(Form("hHT%dJpsiPtInvM2",trig));
			emb2d[trig]->Write();
		}
		invc->cd(trig+1);
		datainvm[trig] = (TH1F*)data2d[trig]->ProjectionY(Form("datainvm%d",trig));	
		datainvm[trig]->SetTitle(Form("HT%d Invariant mass comparison;m_{ee} GeV/c^2;",trig));
		cout<<"datainvm integral ====="<<datainvm[trig]->Integral(21,24)<<endl;
		embinvm[trig] = (TH1F*)emb2d[trig]->ProjectionY(Form("embinvm%d",trig));
		embinvm[trig]->SetLineColor(kRed);
		Int_t xmin,xmax;

		if(inm==29){
			xmin = datainvm[trig]->GetXaxis()->FindBin(2.9);
			xmax = datainvm[trig]->GetXaxis()->FindBin(3.2);
		}
		else if(inm==30){
			xmin = datainvm[trig]->GetXaxis()->FindBin(3.0);
			xmax = datainvm[trig]->GetXaxis()->FindBin(3.15);
		}
		if(eid==2) embinvm[trig]->Rebin(10);
		embinvm[trig]->Scale(datainvm[trig]->Integral(xmin,xmax)/embinvm[trig]->Integral(xmin,xmax));
		cout<<"   bin   "<<embinvm[trig]->GetNbinsX()<<"    "<<datainvm[trig]->GetNbinsX()<<endl;
		//	embinvm[trig]->Scale(datainvm[trig]->GetMaximum()/embinvm[trig]->GetMaximum());
		datainvm[trig]->SetMaximum((embinvm[trig]->GetMaximum()*1.2));
		datainvm[trig]->Draw();
		datainvm[trig]->Write();
		embinvm[trig]->Draw("same HIST");
		embinvm[trig]->Write();
		invleg[trig] = new TLegend(0.7,0.7,0.9,0.9);
		invleg[trig]->AddEntry(datainvm[trig],"data","lep");
		invleg[trig]->AddEntry(embinvm[trig],"emb","l");
		invleg[trig]->Draw("same");

		ptc->cd(trig+1);
		Int_t ymin,ymax;
		ymin = data2d[trig]->GetYaxis()->FindBin(3.0);
		ymax = data2d[trig]->GetYaxis()->FindBin(3.2)-1;

		Int_t ptmin,ptmax;
		ptmin = data2d[trig]->GetXaxis()->FindBin(0);
		ptmax = data2d[trig]->GetXaxis()->FindBin(15);

		data2d[trig]->GetXaxis()->SetRangeUser(0,15);

		datapt[trig] = (TH1F*)data2d[trig]->ProjectionX(Form("datapt%d",trig),ymin,ymax);
		datapt[trig]->SetTitle(Form("HT%d p_{T} shape comparison;p_{T} GeV/c;",trig));
		cout<<"integral datapt ====="<<datapt[trig]->Integral()<<endl;
		//datapt[trig]->GetXaxis()->SetLimits(0,10);	
		//		datapt[trig]->RebinX(2);
		datapt[trig]->Write();
		embpt[trig] = (TH1F*)emb2d[trig]->ProjectionX(Form("embpt%d",trig));
		//		embpt[trig]->GetXaxis()->SetLimits(0,10);
		embpt[trig]->SetLineColor(kRed);
		embpt[trig]->Write();

		datapt[trig]->Draw();

		embpt[trig]->Scale(datapt[trig]->Integral(ptmin,ptmax)/embpt[trig]->Integral(ptmin,ptmax));
		//		if(trig==1) embpt[trig]->Scale(datapt[trig]->GetMaximum()/embpt[trig]->GetMaximum());
		if(trig==1) embpt[trig]->Scale(datapt[trig]->Integral()/embpt[trig]->Integral());
		//		embpt[trig]->Scale(datapt[trig]->GetMaximum()/embpt[trig]->GetMaximum());
		embpt[trig]->Draw("same HIST");
		ptleg[trig] = new TLegend(0.7,0.7,0.9,0.9);
		ptleg[trig]->AddEntry(datapt[trig],"data","lep");
		ptleg[trig]->AddEntry(embpt[trig],"emb","l");
		ptleg[trig]->Draw("same");
		pttxt.SetNDC();
		pttxt.SetTextSize(0.04);
		pttxt.DrawLatex(0.5,0.6,"3.0 < m_{ee} < 3.2 GeV/c^2");
		pttxt.DrawLatex(0.5,0.5,Form("No.Jpsi = %.0f",datapt[trig]->Integral()));

		c[trig] = new TCanvas(Form("c%d",trig),Form("c%d",trig),1200,800);
		c[trig]->Divide(3,2);

		//		Int_t xmin,xmax;	
		for(int pt=0;pt<6;pt++){
			inv[trig][pt][0] = (TH1F*)infile1->Get(Form("pro_trig%d_pT%d",trig,pt));
			inv[trig][pt][0]->SetName(Form("pro_trig%d_pT%d",trig,pt));
			inv[trig][pt][1] = (TH1F*)infile2->Get(Form("trig%d_pt%d",trig+1,pt));
			inv[trig][pt][1]->SetName(Form("trig%d_pt%d",trig+1,pt));
			//			inv[trig][pt][1]->Scale(inv[trig][pt][0]->GetMaximum()/inv[trig][pt][1]->GetMaximum());
			inv[trig][pt][1]->Scale(inv[trig][pt][0]->Integral(xmin,xmax-1)/inv[trig][pt][1]->Integral(xmin,xmax-1));
			inv[trig][pt][0]->SetMaximum(inv[trig][pt][1]->GetMaximum()*1.5);
			inv[trig][pt][0]->Write();
			inv[trig][pt][1]->Write();

			c[trig]->cd(pt+1);	
			inv[trig][pt][0]->SetTitle(Form("%s;m_{ee}",ptbin[pt].Data()));
			//			inv[trig][pt][0]->SetMaximum(inv[trig][pt][1]->GetMaximum()*1.5);
			inv[trig][pt][0]->Draw();
			inv[trig][pt][1]->Draw("HIST same");
			leg[trig] = new TLegend(0.7,0.7,0.9,0.9);
			leg[trig]->AddEntry(inv[trig][pt][1],Form("emb %.f J/#psi",inv[trig][pt][1]->Integral(xmin,xmax-1)),"l");
			leg[trig]->AddEntry(inv[trig][pt][0],Form("dat %.f J/#psi",inv[trig][pt][0]->Integral(xmin,xmax-1)),"lep");
			leg[trig]->Draw("same");
		}
		c[trig]->SaveAs(Form("~/WWW/invMass/%deid_inv%d/pdf20170921/canvas_invariant_ht%d.pdf",eid,inm,trig));
		c[trig]->Write();
	}
	invc->Write();
	invc->SaveAs(Form("~/WWW/invMass/%deid_inv%d/pdf20170921/invc.pdf",eid,inm));
	ptc->Write();
	//	ptc->SaveAs("~/WWW/invMass/1eid_inv30/pdf20170912/ptc.pdf");
}
