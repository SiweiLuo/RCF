#define NFRAME 2
#define NPT 6

TString trigSet[3] = {"MB","HT0","HT2"};
TString histType[4] = {"unlikeSign","likeSign","signal","eff"};
//TString pTbin[3] = {"0 < p_{T} < 2 GeV/c","2 < p_{T} < 3 GeV/c","6 < p_{T} < 8 GeV/c"};
TString pTbin[6] = {"0 < p_{T} < 2GeV/c","2 < p_{T} < 3GeV/c","3 < p_{T} < 4GeV/c","4 < p_{T} < 6GeV/c","6 < p_{T} < 8GeV/c","8 < p_{T} < 14GeV/c"};
TString FRM[2] = {"HX","CS"};

void PWGproposal(){

	gStyle->SetOptStat(false);
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetPadLeftMargin(0.15);

	//	gStyle->SetTitleSize(1.2);
//	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFontSize(0.08);
//	gStyle->SetTitleFillColor(-1);

	gStyle->SetErrorX(0);

	TFile *infilept[NPT][NFRAME];
	/*
	   infilept[0] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_0_0_0_0_20171025.root","read");
	   infilept[1] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_1_1_0_0_20171025.root","read");
	   infilept[2] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_3_4_0_0_20171025.root","read");
	   */
	for(int frame =0;frame<2;frame++){
		for(int pt=0,trig=0;pt<6;pt++){
			if(pt==0) trig=0;
			else if(pt>=1 && pt<=2) trig=1;
			else if(pt>=3 && pt<=5) trig=3;
			infilept[pt][frame] = new TFile(Form("~/polresults/20160707/functional/splot_3D_20171116_2eid_inv30/splot_3D_%d_%d_%d_0_20171116.root",trig,pt,frame),"read");
		}
	}

	TFile* infile[3];
	/*   	infile[0] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_0_0_0_0_20171025.root","read");
			infile[1] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_1_1_0_0_20171025.root","read");
			infile[2] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_3_4_0_0_20171025.root","read");
			*/

	infile[0] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht-1_sys0.ana.root","read");
	infile[1] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht0_sys0.ana.root","read");
	infile[2] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht2_sys0.ana.root","read");

	TH2F *unlike[3],*like[3];
	//unlike[0] = infile[0]->Get("tothist");
	//unlike[0]->SetName("MB_unlike_sign");

	for(int i=0;i<3;i++){
		unlike[i] = (TH2F*)((TH3F*)infile[i]->Get("hJpsiCosThetaPhiPt"))->Project3D("xy");
		unlike[i]->SetName(Form("%s_unlike_sign",trigSet[i].Data()));
		unlike[i]->SetTitle(Form("%s unlike sign distribution",trigSet[i].Data()));
		unlike[i]->RebinX(4);
		unlike[i]->RebinY(4);
		unlike[i]->SetLabelSize(0.08);	
		like[i] = (TH2F*)((TH3F*)infile[i]->Get("hJpsiCosThetaPhiPtBG"))->Project3D("xy");
		like[i]->SetName(Form("%s_like_sign",trigSet[i].Data()));
		like[i]->SetTitle(Form("%s like sign distribution",trigSet[i].Data()));
		like[i]->RebinX(4);
		like[i]->RebinY(4);
		like[i]->SetLabelSize(.2);

		//		efficiency[i] = (TH2F*)infile[i]->Get("effhist");
		//		efficiency[i]->SetName(Form("%s_eff",trigSet[i].Data()));
	}

	TCanvas *dataeff[2];//unlike sign ; like sign  
	for(int i=0;i<2;i++){
		dataeff[i] = new TCanvas(Form("dataeff_%d",i),Form("dataeff_%d",i),1200,400);
		dataeff[i]->Divide(3,1,0,0);
		for(int j=0;j<3;j++){
			dataeff[i]->cd(j+1);	
			if(i==0){
				unlike[j]->Scale(1./unlike[j]->Integral()); 
				unlike[j]->Draw("box");
				//				unlike[j]->Draw("colz");
			}	
			if(i==1){
				like[j]->Scale(1./like[j]->Integral());
				like[j]->Draw("box");
			}
			//			if(i==2) efficiency[j]->Draw("colz");
		}
		dataeff[i]->SaveAs(Form("~/WWW/Proposal/figure/%s.pdf",histType[i].Data()));
	}	

	TH2F *unlikept[6][2],*likept[NPT][NFRAME],*efficiency[NPT][NFRAME],*signalpt[NPT][NFRAME];
	TH1F *signalpt_py[NPT][NFRAME],*signalpt_px[NPT][NFRAME];


	TCanvas *dataeffpt[4];
	/*
	   for(int i=0;i<3;i++){
	   dataeffpt[i] = new TCanvas(Form("dataeffpt_%d",i),Form("dataeffpt_%d",i),1200,400);
	   dataeffpt[i]->Divide(3,1);
	   for(int j=0;j<3;j++){
	   unlikept[j] = (TH2F*)infilept[j]->Get("tothist");
	   unlikept[j]->SetName("%s_unlike_sign_pt");
	   unlikept[j]->SetTitle(Form("%s %s unlike sign distribution;#phi;cos#theta",trigSet[j].Data(),pTbin[j].Data()));
	   likept[j] = (TH2F*)infilept[j]->Get("bkghist");
	   likept[j]->SetName("%s_like_sign_pt");
	   likept[j]->SetTitle(Form("%s %s like sign distribution;#phi;cos#theta",trigSet[j].Data(),pTbin[j].Data()));
	   efficiency[j] = (TH2F*)infilept[j]->Get("effhist");
	   efficiency[j]->SetName("%s_efficiency_pt");
	   efficiency[j]->SetTitle(Form("%s %s efficiency;#phi;cos#theta",trigSet[j].Data(),pTbin[j].Data()));
	   dataeffpt[i]->cd(j+1);
	   if(i==0)unlikept[j]->Draw("colz");
	   if(i==1)likept[j]->Draw("colz");
	   if(i==2)efficiency[j]->Draw("box");
	   }
	   dataeffpt[i]->SaveAs(Form("~/WWW/Proposal/figure/%s_pT.pdf",histType[i].Data()));	
	   }
	   */


	for(int frame=0;frame<2;frame++){
		for(int pt=0,trig=0;pt<6;pt++){
			if(pt==0) trig=0;
			else if(pt>=1 && pt<=2) trig=1;
			else if(pt>=3 && pt<=5) trig=2;
			unlikept[pt][frame] = (TH2F*)infilept[pt][frame]->Get("tothist");
			unlikept[pt][frame]->SetName(Form("%d_%d_unlike",pt,frame));
			unlikept[pt][frame]->SetTitle(Form("%s %s in %s;#phi;cos#theta",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			//			unlikept[pt][frame]->SetTitleSize(1.,"t");

			likept[pt][frame] = (TH2F*)infilept[pt][frame]->Get("bkghist");
			likept[pt][frame]->SetName(Form("%d_%d_like",pt,frame));
			likept[pt][frame]->SetTitle(Form("%s %s in %s;#phi;cos#theta",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			//			likept[pt][frame]->SetTitleSize(1.,"t");

			signalpt[pt][frame] = (TH2F*)infilept[pt][frame]->Get("datahist");
			signalpt[pt][frame]->SetName(Form("%d_%d_signal",pt,frame));
//			signalpt[pt][frame]->SetTitle(Form("%s %s in %s;#phi;cos#theta",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			signalpt[pt][frame]->SetTitle(Form(";#phi;cos#theta"));
			signalpt[pt][frame]->GetXaxis()->SetTitleOffset(0.4);
			signalpt[pt][frame]->GetXaxis()->SetTitleSize(0.1);
			signalpt[pt][frame]->GetYaxis()->SetTitleOffset(0.4);
			signalpt[pt][frame]->GetYaxis()->SetTitleSize(0.1);


			//signalpt[pt][frame]->SetLabelSize(0.06);
			//			signalpt[pt][frame]->SetTitleSize(1.,"t");

			signalpt_py[pt][frame] = (TH1F*)signalpt[pt][frame]->ProjectionY(Form("singalpt_%d_%d_py",pt,frame));
//			signalpt_py[pt][frame]->SetTitle(Form("%s %s in %s;cos#theta;",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			signalpt_py[pt][frame]->SetTitle(Form(";cos#theta;"));
			signalpt_py[pt][frame]->GetXaxis()->SetTitleOffset(0.4);
			signalpt_py[pt][frame]->GetXaxis()->SetTitleSize(0.1);
			//			signalpt_py[pt][frame]->SetTitleSize(1.,"t");
			signalpt_px[pt][frame] = (TH1F*)signalpt[pt][frame]->ProjectionX(Form("singalpt_%d_%d_px",pt,frame));
//			signalpt_px[pt][frame]->SetTitle(Form("%s %s in %s;#phi;",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			signalpt_px[pt][frame]->SetTitle(Form(";#phi;"));
			signalpt_px[pt][frame]->GetXaxis()->SetTitleOffset(0.4);
			signalpt_px[pt][frame]->GetXaxis()->SetTitleSize(0.1);
			//			signalpt_px[pt][frame]->SetTitleSize(1.,"t");

			efficiency[pt][frame] = (TH2F*)infilept[pt][frame]->Get("effhist");
			efficiency[pt][frame]->SetName(Form("%d_%d_eff",pt,frame));
//			efficiency[pt][frame]->SetTitle(Form("%s %s in %s;#phi;cos#theta",trigSet[trig].Data(),pTbin[pt].Data(),FRM[frame].Data()));
			efficiency[pt][frame]->SetTitle(Form(";#phi;cos#theta"));
			efficiency[pt][frame]->GetXaxis()->SetTitleOffset(0.4);
			efficiency[pt][frame]->GetXaxis()->SetTitleSize(0.1);
			efficiency[pt][frame]->GetYaxis()->SetTitleOffset(0.5);
			efficiency[pt][frame]->GetYaxis()->SetTitleSize(0.1);

			
			
			//			efficiency[pt][frame]->SetLabelSize(0.04);
			//			efficiency[pt][frame]->SetTitleSize(1.,"t");
		}
	}

	for(int type=0;type<4;type++){
		dataeffpt[type] = new TCanvas(Form("dataeffpt_%d",type),Form("dataeffpt_%d",type),1200,400);
		dataeffpt[type]->Divide(6,2,0,0);
		for(int frame=0;frame<2;frame++){
			for(int pt=0;pt<6;pt++){	
				dataeffpt[type]->cd(6*frame+pt+1);
				gPad->SetTickx(2);
				if(type==0)unlikept[pt][frame]->Draw("box");
				if(type==1)likept[pt][frame]->Draw("box");
				if(type==2){
					//					signalpt[pt][frame]->GetZaxis()->SetRangeUser(-5,30);
					//	signalpt[pt][frame]->Draw("colz");
					signalpt[pt][frame]->Scale(1./signalpt[pt][frame]->Integral());
					//					cout<<" minimum = "<<signalpt[pt][frame]->GetMinimum()<<endl;
					//					cout<<" maximum = "<<signalpt[pt][frame]->GetMaximum()<<endl;
					//					signalpt[pt][frame]->GetZaxis()->SetRangeUser(-0.012,0.08);
					signalpt[pt][frame]->Draw("box z");
					//					signalpt[pt][frame]->Draw("z");
					//					signalpt[pt][frame]->Draw("colz");
					//					cout<<" maximum ======== "<<signalpt[pt][frame]->GetMaximum()<<endl;
				}
				if(type==3)efficiency[pt][frame]->Draw("box");
			}
		}
		dataeffpt[type]->SaveAs(Form("~/WWW/Proposal/figure/%s.pdf",histType[type].Data()));
	}

	TFile *comparisonfile[6][2];
	TCanvas *allcomparison = new TCanvas("allcomparison","allcomparison",1200,400);
	//	allcomparison->Divide(6,2,0.5,0.5);
	allcomparison->Divide(6,2,0,0);

	TCanvas *allcomparison_costheta = new TCanvas("allcomparison_costheta","allcomparison_costheta",1200,400);
	allcomparison_costheta->Divide(6,2,0,0);


	TLegend *allcomparisonleg;
	TCanvas *subcanvas[6][2][2];
	//	TCanvas *subcanvas[6][2];
	TPad *subpad[6][2][2];

	TCanvas *comparison1d = new TCanvas("comparison1d","comparison1d",1200,800);
	comparison1d->Divide(6,4,0,0);
	TH1F *th1f[6][2][2][3];


	for(int frame=0;frame<2;frame++){
		for(int pt=0,trig=0;pt<6;pt++){
			if(pt==0) trig=0;
			else if(pt>=1 && pt<=2) trig=1;
			else if(pt>=3 && pt<=5) trig=3;
			//			comparisonfile[pt][frame] = new TFile(Form("~/polresults/20160707/functional/splot_3D_20171221_2eid_inv30_scan_migrad_hess/splot_3D_%d_%d_%d_0_20171221.root",trig,pt,frame),"read");
			comparisonfile[pt][frame] = new TFile(Form("~/polresults/20160707/functional/splot_3D_20180102_2eid_inv30_scan_migrad_hess/splot_3D_%d_%d_%d_0_20180102.root",trig,pt,frame),"read");

			subcanvas[pt][frame][0] = (TCanvas*)comparisonfile[pt][frame]->Get("comparison1");	
			subcanvas[pt][frame][0]->SetName(Form("subcanvas1_%d_%d",pt,frame));

			subcanvas[pt][frame][1] = (TCanvas*)comparisonfile[pt][frame]->Get("comparison2");	
			subcanvas[pt][frame][1]->SetName(Form("subcanvas2_%d_%d",pt,frame));

			th1f[pt][frame][0][0] = (TH1F*)comparisonfile[pt][frame]->Get("py");
			th1f[pt][frame][0][1] = (TH1F*)comparisonfile[pt][frame]->Get("hsigthetaplus_py");			
			th1f[pt][frame][0][2] = (TH1F*)comparisonfile[pt][frame]->Get("hsigthetaminus_py");			
			th1f[pt][frame][0][1]->SetLineStyle(3);
			th1f[pt][frame][0][2]->SetLineStyle(3);

			th1f[pt][frame][1][0] = (TH1F*)comparisonfile[pt][frame]->Get("px");
			th1f[pt][frame][1][1] = (TH1F*)comparisonfile[pt][frame]->Get("hsigphiplus_px");			
			th1f[pt][frame][1][2] = (TH1F*)comparisonfile[pt][frame]->Get("hsigphiminus_px");			
			th1f[pt][frame][1][1]->SetLineStyle(3);
			th1f[pt][frame][1][2]->SetLineStyle(3);

			//			subcanvas[pt][frame] = (TCanvas*)comparisonfile[pt][frame]->Get("comparison");
			//			cout<<" pointer = "<<subcanvas[pt][frame]<<endl;
			//			subcanvas[pt][frame]->SetName(Form("subcanvas_%d_%d",pt,frame));
			//			subpad[pt][frame][0] = subcanvas[pt][frame]->GetPadPointer(1)->Clone();
			//			subpad[pt][frame][0]->SetName(Form("pad1_%d_%d",pt,frame));

			allcomparison_costheta->cd(pt+1+6*frame);
			//			TPad *allcomparison_costheta_pad = new TPad(Form("pad1_%d_%d",pt,frame),"",0,0,100,100);
			//			allcomparison_costheta_pad->Range(-2,-4,10,10);
			//			allcomparison_costheta_pad->cd();
			//			gPad->SetTickx(2); 
			//			subcanvas[pt][frame]->Draw();
			//			subcanvas[pt][frame][0]->DrawClonePad();			
			//			subpad[pt][frame][0]->DrawClonePad();

			signalpt_py[pt][frame]->Scale(1./signalpt_py[pt][frame]->Integral());
			signalpt_py[pt][frame]->GetYaxis()->SetRangeUser(0,0.35);	
			signalpt_py[pt][frame]->SetMarkerStyle(29);
			signalpt_py[pt][frame]->Draw("");

			th1f[pt][frame][0][0]->Scale(1./th1f[pt][frame][0][0]->Integral());
			th1f[pt][frame][0][1]->Scale(1./th1f[pt][frame][0][1]->Integral());
			th1f[pt][frame][0][2]->Scale(1./th1f[pt][frame][0][2]->Integral());

			th1f[pt][frame][0][0]->Draw("hist same");
			th1f[pt][frame][0][1]->Draw("hist same");
			th1f[pt][frame][0][2]->Draw("hist same");
			if(pt==5 && frame==0){
				allcomparisonleg = new TLegend(0.1,0.55,0.6,0.94);
				allcomparisonleg->SetBorderSize(0);
				allcomparisonleg->AddEntry(signalpt_px[pt][frame],"raw data","lep");
				allcomparisonleg->AddEntry(th1f[pt][frame][0][0],"default fit","l");
				allcomparisonleg->AddEntry(th1f[pt][frame][0][1],"#lambda_{#theta} = +1","l");
				allcomparisonleg->AddEntry(th1f[pt][frame][0][2],"#lambda_{#theta} = -1","l");
				allcomparisonleg->Draw("same");
			}

			allcomparison->cd(pt+1+6*frame);

			//			gPad->SetTickx(2); 
			//			subcanvas[pt][frame][1]->DrawClonePad();
			//			subpad[pt][frame][0]->Draw();
			signalpt_px[pt][frame]->Scale(1./signalpt_px[pt][frame]->Integral());	
			signalpt_px[pt][frame]->GetYaxis()->SetRangeUser(0,0.35);
			signalpt_px[pt][frame]->SetMarkerStyle(29);
			//			signalpt_px[pt][frame]->Draw("ap e1");
			signalpt_px[pt][frame]->Draw("e1");

			th1f[pt][frame][1][0]->Scale(1./th1f[pt][frame][1][0]->Integral());
			th1f[pt][frame][1][1]->Scale(1./th1f[pt][frame][1][1]->Integral());
			th1f[pt][frame][1][2]->Scale(1./th1f[pt][frame][1][2]->Integral());

			th1f[pt][frame][1][0]->Draw("hist same");
			th1f[pt][frame][1][1]->Draw("hist same");
			th1f[pt][frame][1][2]->Draw("hist same");
			if(pt==5 && frame==0){
				allcomparisonleg = new TLegend(0.1,0.55,0.6,0.94);
				allcomparisonleg->SetBorderSize(0);
				allcomparisonleg->AddEntry(signalpt_px[pt][frame],"raw data","lep");
				allcomparisonleg->AddEntry(th1f[pt][frame][1][0],"default fit","l");
				allcomparisonleg->AddEntry(th1f[pt][frame][1][1],"#lambda_{#phi} = +1","l");
				allcomparisonleg->AddEntry(th1f[pt][frame][1][2],"#lambda_{#phi} = -1","l");
				allcomparisonleg->Draw("same");
			}

			comparison1d->cd(pt+7+12*frame);
			th1f[pt][frame][0][0]->Draw();
			th1f[pt][frame][0][1]->Draw("same hist");
			th1f[pt][frame][0][2]->Draw("same hist");

			comparison1d->cd(pt+1+12*frame);
			th1f[pt][frame][1][0]->Draw();
			th1f[pt][frame][1][1]->Draw("same hist");
			th1f[pt][frame][1][2]->Draw("same hist");
		}
	}

	allcomparison->SaveAs("~/WWW/Proposal/figure/allcomparison.pdf");
	allcomparison_costheta->SaveAs("~/WWW/Proposal/figure/allcomparison_costheta.pdf");
	comparison1d->SaveAs("~/WWW/Proposal/figure/comparison1d.pdf");

}


void read_eff(){








}
