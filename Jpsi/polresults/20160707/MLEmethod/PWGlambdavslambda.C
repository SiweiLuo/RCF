#include <algorithm>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFormula.h"
#include "TRandom.h"
#include "TDatime.h"

#define NEXP     1000
//#define NSIGNAL  500
//#define NBKG     200
//#define LBTHESIG -0.5
//#define LBPHISIG +0.4
#define LBTHEBKG +0.6
#define LBPHIBKG -0.3
#define XMIN    -3.1415927
#define XMAX    +3.1415927
//#define YMIN    -3.1415927
//#define YMAX    +3.1415927
#define YMIN -1
#define YMAX 1 
#define NBIN    10

Int_t NSIGNAL = 500;
Int_t NBKG = 200;

//TF2* fsig=new TF2("fsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);
//TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);

TH2F *hunl, *hlik, *hsig, *hbkg, *hsub;
TH2F *hlambda, *hlambda2, *hlambda_err, *hlambda2_err;
TGraph *gc, *gc2;

TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);

#define NFILE 100
#define NTRIG 4
#define NPT 6
#define NPHASE 2 
#define NFRAME 2
#define NREBIN 4
#define XNBIN 40
#define YNBIN 40
#define CHIXBIN 100
#define CHIYBIN 100
//#define rebin 0 //1

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigSet[NTRIG] = {"mb","ht0","ht1","ht2"};

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME];
TH2F *bkgdata2D[NFILE][NTRIG][NPT][NFRAME];
TH2F *totdata2D[NFILE][NTRIG][NPT][NFRAME];
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][2]; // pass , total ;
TH2F *chi2result[NFILE][NTRIG][NPT][NFRAME];

TFile *efficiencyfile;
TFile *rawdatafile[NFILE][NTRIG];
TFile *rawdatafile1[NFILE][NTRIG];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker

TH2F* mle[NTRIG][NPT][NFRAME];
TFile* histograms;
TFile* histograms_copy;
TFile* likelihoodhistograms[NTRIG][NPT][NFRAME];

//float lambda[2][2]; // theta, phi; initial, terminal;

TF2* sigma; // recover from observation 
TF2* sigmafcn;// generate mc data 
TH2F* effhist;
TH2F* datahist;
TH2F* bkghist;
TH2F* tothist;

TFile* outputfile;
TFile* lambdafile;

void PWGlambdavslambda(){

	TString trigName[NTRIG] = {"MB","HT0","HT1","HT2"};
	TString pTName[NPT] = {"0 < p_{T} < 2","2 < p_{T} < 3","3 < p_{T} < 4","4 < p_{T} < 6","6 < p_{T} < 8","8 < p_{T} < 12"};
	TString frameName[NFRAME] = {"HX frame","CS frame"};
	int lthetaphi=0;

	TDatime* time = new TDatime();

	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptDate(0);
	gStyle->SetOptStat();
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);

	TFile *infile, *infilep, *infilen;
	TH2F* h2;
	TH2F* h2p;
	TH2F* h2n;
	TH2F* h2err;
	TH2F* h2perr;
	TH2F* h2nerr;

	TF1* linear1 = new TF1("linear1","x",-1,1);
	linear1->SetLineColor(kRed);
	TF1* fit1 = new TF1("fit1","[0]*x+[1]",-1,1);

	TF1* linear2 = new TF1("linear2","x",-1.,1.);
	linear2->SetLineColor(kRed);
	TF1* fit2 = new TF1("fit2","[0]*x+[1]",-1,1);

	TF1* linear3 = new TF1("linear3","x",0,0.6);
	linear3->SetLineColor(kRed);
	TF1* fit3 = new TF1("fit3","[0]*x+[1]",-1,1);

	TF1* linear4 = new TF1("linear4","x",0,0.6);
	linear4->SetLineColor(kRed);
	TF1* fit4 = new TF1("fit4","[0]*x",-1,1);

	TF1* linear5 = new TF1("linear5","x",-1.,1.);
	linear5->SetLineColor(kRed);
	TF1* fit5 = new TF1("fit5","[0]*x",-1,1);

	TF1* linear6 = new TF1("linear6","x",0,0.6);
	linear6->SetLineColor(kRed);
	TF1* fit6 = new TF1("fit6","[0]*x",-1,1);

	Double_t ltheta[5],ltheta_err[5],mltheta[5],mltheta_err[5],sigltheta_err[5],sigltheta_err_err[5];
	Double_t lphi[5],lphi_err[5],mlphi[5],mlphi_err[5],siglphi_err[5],siglphi_err_err[5];

	Double_t mltheta1p[5],mltheta1p_err[5],sigltheta1p_err[5],sigltheta1p_err_err[5];
	Double_t lphi1p[5],lphi1p_err[5],mlphi1p[5],mlphi1p_err[5],siglphi1p_err[5],siglphi1p_err_err[5];

	Double_t mltheta1n[5],mltheta1n_err[5],sigltheta1n_err[5],sigltheta1n_err_err[5];
	Double_t lphi1n[5],lphi1n_err[5],mlphi1n[5],mlphi1n_err[5],siglphi1n_err[5],siglphi1n_err_err[5];

	Double_t lthetaphii[5],lthetaphii_err[5],mlthetaphi[5],mlthetaphi_err[5],siglthetaphi_err[5],siglthetaphi_err_err[5];

	TCanvas* meanc = new TCanvas("meanc","meanc",2400,2400);
	meanc->Divide(6,6,0,0);
	TCanvas* sigmac = new TCanvas("sigmac","sigmac",2400,2400);
	sigmac->Divide(6,6,0,0);

	Double_t LBTHESIG,LBPHISIG,LBTHEPHISIG;
	TH3F* lambda;
	TGraphErrors *thetavstheta[6][2], *phivsphi[6][2],*thetaphivsthetaphi[6][2],*sigmatheta[6][2],*sigmaphi[6][2],*sigmathetaphi[6][2];
	for(int frame=0;frame<2;frame++){
		for(int pt=0,trig=0;pt<6;pt++){
			if(pt==0) trig=0;
			else if(pt>=1 && pt<=2) trig=1;
			else if(pt>=3 && pt<=5) trig=3;
			
			infile = new TFile(Form("/star/u/siwei/polresults/20160707/functional/splot_3D_20171116_2eid_inv30/splot_3D_%d_%d_%d_0_20171116.root",trig,pt,frame),"read");
			lambda = (TH3F*)infile->Get("hlambda");
			lambda->SetName(Form());

			if(lambda==0x0) return;
			else {
				LBTHESIG = lambda->GetMean(1);
				LBPHISIG = lambda->GetMean(2);
				LBTHEPHISIG = lambda->GetMean(3);
			}

			delete infile;

			for(int i=-2;i<=2;i++){
				infile = new TFile(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171213_2eid_inv30_wofactor_scan_hess/splot_3D_%d_%d_%d_%d_0_0_20171213.root",trig,pt,frame,i),"read");
				h2 = (TH2F*)infile->Get("hlambda");
				if(h2==0){
					delete infile;
					continue;
				}
				mltheta[i+2] = h2->GetMean(1);
				mltheta_err[i+2] = h2->GetMeanError(1);
				ltheta[i+2] = LBTHESIG+0.2*i;
				ltheta_err[i+2] = h2->GetRMS(1);
				cout<<"ltheta_err = "<<ltheta_err[i+2]<<endl;
				h2err = (TH2F*)infile->Get("hlambda_err");
				sigltheta_err[i+2] = h2err->GetMean(1);
				sigltheta_err_err[i+2] = h2err->GetMeanError(1);
				delete h2;
				delete infile;
			}

			for(int i=-2;i<=2;i++){
				infile = new TFile(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171213_2eid_inv30_wofactor_scan_hess/splot_3D_%d_%d_%d_0_%d_0_20171213.root",trig,pt,frame,i),"read");
				h2 = (TH2F*)infile->Get("hlambda");
				if(h2==0) {
					delete infile;
					continue;
				}
				mlphi[i+2] = h2->GetMean(2);
				mlphi_err[i+2] = h2->GetMeanError(2);
				lphi[i+2] = LBPHISIG+0.2*i;
				lphi_err[i+2] = h2->GetRMS(2);

				h2err = (TH2F*)infile->Get("hlambda_err");
				siglphi_err[i+2] = h2err->GetMean(2);
				siglphi_err_err[i+2] = h2err->GetMeanError(2);
				delete h2;
				delete infile;
			}
			
			for(int i=-2;i<=2;i++){
				infile = new TFile(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171213_2eid_inv30_wofactor_scan_hess/splot_3D_%d_%d_%d_0_0_%d_20171213.root",trig,pt,frame,i),"read");
				cout<<"file ====="<<Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171213_2eid_inv30_wofactor_scan_hess/splot_3D_%d_%d_%d_0_0_%d_20171213.root",trig,pt,frame,i)<<endl;
				h2 = (TH2F*)infile->Get("hlambda");
				if(h2==0){
					cout<<"*** no histogram ***"<<endl;
					delete infile;
					continue;
				}
				mlthetaphi[i+2] = h2->GetMean(3);
				mlthetaphi_err[i+2] = h2->GetMeanError(3);

				cout<<" lambda theta phi ===="<<mlthetaphi[i+2]<<endl;
				lthetaphii[i+2] = LBTHEPHISIG+0.2*i;
				lthetaphii_err[i+2] = h2->GetRMS(3);

				h2err = (TH2F*)infile->Get("hlambda_err");
				siglthetaphi_err[i+2] = h2err->GetMean(3);
				siglthetaphi_err_err[i+2] = h2err->GetMeanError(3);
				cout<<" sigma lambda theta phi ="<<siglthetaphi_err[i+2]<<endl;
				cout<<" i=== "<<i<<endl;
				delete h2;
				delete infile;
			}
		cout<<"111111111111111111111111"<<endl;
			thetavstheta[pt][frame] = new TGraphErrors(5,ltheta,mltheta,0,mltheta_err);
			thetavstheta[pt][frame]->SetName(Form("theta_%d_%d",pt,frame));
			//	thetavstheta[pt][frame]->SetTitle("measurements of #lambda_{#theta};true #lambda_{#theta};measured #lambda_{#theta}");
			thetavstheta[pt][frame]->SetTitle(Form("%s @ %s GeV/c in %s",trigName[trig].Data(),pTName[pt].Data(),frameName[frame].Data()));
			thetavstheta[pt][frame]->GetXaxis()->SetLimits(-3,3);
			thetavstheta[pt][frame]->GetYaxis()->SetRangeUser(-3,3);
			thetavstheta[pt][frame]->SetMarkerStyle(29);
			thetavstheta[pt][frame]->SetMarkerSize(5);
			thetavstheta[pt][frame]->SetMarkerColor(kBlue);
			thetavstheta[pt][frame]->SetLineColor(kBlue);
			phivsphi[pt][frame] = new TGraphErrors(5,lphi,mlphi,0,mlphi_err);
			phivsphi[pt][frame]->SetName(Form("phi_%d_%d",pt,frame));
			phivsphi[pt][frame]->SetTitle("measurements of #lambda_{#phi};true #lambda_{#phi};measured #lambda_{#phi}");
			phivsphi[pt][frame]->GetXaxis()->SetLimits(-3,3);
			phivsphi[pt][frame]->GetYaxis()->SetRangeUser(-3,3);
			phivsphi[pt][frame]->SetMarkerStyle(29);
			phivsphi[pt][frame]->SetMarkerSize(5);
			phivsphi[pt][frame]->SetMarkerColor(kBlue);
			phivsphi[pt][frame]->SetLineColor(kBlue);
			
			thetaphivsthetaphi[pt][frame] = new TGraphErrors(5,lthetaphii,mlthetaphi,0,mlthetaphi_err);
			thetaphivsthetaphi[pt][frame]->GetXaxis()->SetLimits(-3,3);
			thetaphivsthetaphi[pt][frame]->GetYaxis()->SetRangeUser(-3,3);
			thetaphivsthetaphi[pt][frame]->SetMarkerStyle(29);
			thetaphivsthetaphi[pt][frame]->SetMarkerSize(5);
			thetaphivsthetaphi[pt][frame]->SetMarkerColor(kBlue);
			thetaphivsthetaphi[pt][frame]->SetLineColor(kBlue);
			thetaphivsthetaphi[pt][frame]->SetTitle("measurements of #lambda_{#theta#phi};true #lambda_{#theta#phi};measured #lambda_{#theta#phi}");

			sigmatheta[pt][frame] = new TGraphErrors(5,ltheta_err,sigltheta_err,0,sigltheta_err_err);
			sigmatheta[pt][frame]->SetName(Form("sigmatheta_%d_%d",pt,frame));
			sigmatheta[pt][frame]->GetXaxis()->SetLimits(0,1.5);
			sigmatheta[pt][frame]->GetYaxis()->SetRangeUser(0,1.5);
			sigmatheta[pt][frame]->SetMarkerStyle(29);
			sigmatheta[pt][frame]->SetMarkerSize(5);
			sigmatheta[pt][frame]->SetMarkerColor(kBlue);
			sigmatheta[pt][frame]->SetTitle("measurements of #sigma_{#lambda_{#theta}};true #sigma_{#lambda_{#theta}};measured #sigma_{#lambda_{#theta}}");


			sigmaphi[pt][frame] = new TGraphErrors(5,lphi_err,siglphi_err,0,siglphi_err_err);
			//	sigmaphi->GetXaxis()->SetLimits(0,1.2*siglphi_err[4]);
			sigmaphi[pt][frame]->GetXaxis()->SetLimits(0,1.5);
			//	sigmaphi->GetYaxis()->SetRangeUser(0,1.2*siglphi_err[4]);
			sigmaphi[pt][frame]->GetYaxis()->SetRangeUser(0,1.5);
			sigmaphi[pt][frame]->SetMarkerStyle(29);
			sigmaphi[pt][frame]->SetMarkerSize(5);
			sigmaphi[pt][frame]->SetMarkerColor(kBlue);
			sigmaphi[pt][frame]->Print("all");
			sigmaphi[pt][frame]->SetTitle("measurements of #sigma_{#lambda_{#phi}};true #sigma_{#lambda_{#phi}};measured #sigma_{#lambda_{#phi}}");

			sigmathetaphi[pt][frame] = new TGraphErrors(5,lthetaphii_err,siglthetaphi_err,0,siglthetaphi_err_err);
			sigmathetaphi[pt][frame]->GetXaxis()->SetLimits(0,1.5);
			sigmathetaphi[pt][frame]->GetYaxis()->SetRangeUser(0,1.5);
			sigmathetaphi[pt][frame]->SetMarkerStyle(29);
			sigmathetaphi[pt][frame]->SetMarkerSize(5);
			sigmathetaphi[pt][frame]->SetMarkerColor(kBlue);
			sigmathetaphi[pt][frame]->SetTitle("measurements of #sigma_{#lambda_{#theta#phi}};true #sigma_{#lambda_{#theta#phi}};measured #sigma_{#lambda_{#theta#phi}}");
			sigmathetaphi[pt][frame]->Print();

			lambdafile = new TFile(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171213_2eid_inv30_wofactor_scan_hess/PWGthetavstheta_%d_%d_%d_%d.root",trig,pt,frame,time->GetDate()),"recreate");
			//	TCanvas* thetacanvas = new TCanvas("thetacanvas","thetacanvas");
			//thetacanvas->Divide(3,2);
			//thetacanvas->cd(1);
			meanc->cd(18*frame+pt+1);
			gPad->SetTickx(2);
			thetavstheta[pt][frame]->Draw("ap");
			linear1->Draw("same");
			thetavstheta[pt][frame]->Draw("psame");
			thetavstheta[pt][frame]->Fit("fit1");
			meanc->cd(18*frame+pt+7);
			gPad->SetTickx(2);
			//thetacanvas->cd(2);
			phivsphi[pt][frame]->Draw("ap");
			phivsphi[pt][frame]->Fit("fit2");
			linear2->Draw("same");
			phivsphi[pt][frame]->Draw("psame");

			meanc->cd(18*frame+pt+13);
			gPad->SetTickx(2);
			//thetacanvas->cd(3);
			thetaphivsthetaphi[pt][frame]->Draw("ap");
			thetaphivsthetaphi[pt][frame]->Fit("fit3");
			linear5->Draw("same");
			thetaphivsthetaphi[pt][frame]->Draw("psame");
			sigmac->cd(18*frame+pt+1);
			gPad->SetTickx(2);
			//thetacanvas->cd(4);
			sigmatheta[pt][frame]->Draw("ap");
			sigmatheta[pt][frame]->Fit("fit4");
			linear3->Draw("same");
			sigmatheta[pt][frame]->Draw("psame");

			sigmac->cd(18*frame+pt+7);
			gPad->SetTickx(2);
			//thetacanvas->cd(5);
			sigmaphi[pt][frame]->Draw("ap");
			sigmaphi[pt][frame]->Fit("fit5");
			linear4->Draw("same");
			sigmaphi[pt][frame]->Draw("psame");

			sigmac->cd(18*frame+pt+13);
			gPad->SetTickx(2);
			//thetacanvas->cd(6);
			sigmathetaphi[pt][frame]->Draw("ap");
			sigmathetaphi[pt][frame]->Fit("fit6");
			linear6->Draw("same");
			sigmathetaphi[pt][frame]->Draw("psame");	

			//	//thetacanvas->SaveAs(Form("/star/u/siwei/polresults/20160707/MLEmethod/lambdavslambda/%d/thetacanvas_%d_%d_%d_%d.pdf",time->GetDate(),trig,pt,frame,time->GetDate()));
			//thetacanvas->SaveAs(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171116_2eid_inv30_wofactor/thetacanvas_%d_%d_%d_20171116.pdf",trig,pt,frame));
			//	thetacanvas->SaveAs(Form("~/WWW/MonteCarlo/MC/2eid_inv30_lambdasvslambdas/%d/thetacanvas_%d_%d_%d_20171020.pdf",time->GetDate(),trig,pt,frame));
			lambdafile->cd();
						thetavstheta[pt][frame]->Write("",TObject::kOverwrite);
			//thetacanvas[pt][frame]->Write("",TObject::kOverwrite);
			//			fit1[pt][frame]->Write("",TObject::kOverwrite);
			//			fit2[pt][frame]->Write("",TObject::kOverwrite);
			//			fit3[pt][frame]->Write("",TObject::kOverwrite);

			//			fit4[pt][frame]->Write("",TObject::kOverwrite);
			//			fit5[pt][frame]->Write("",TObject::kOverwrite);
			//			fit6[pt][frame]->Write("",TObject::kOverwrite);
			lambdafile->Close();

		}// pt loop
	}// frame loop
		meanc->SaveAs(Form("~/WWW/Proposal/figure/PWGmeanc_%d.pdf",time->GetDate()));
		sigmac->SaveAs(Form("~/WWW/Proposal/figure/PWGsigmac_%d.pdf",time->GetDate()));
}


