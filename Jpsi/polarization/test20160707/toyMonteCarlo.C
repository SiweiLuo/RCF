#include "TH2F.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TError.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <algorithm>

#define NFILE 40
#define NTRIG 3
#define NPT 6
#define NPHASE 2 
#define NFRAME 2
#define NPLOT 6
#define NREBIN 4
#define	nsample 2500
#define XNBIN 10
#define YNBIN 10
#define CHIXBIN 100
#define CHIYBIN 100
#define NEVENT 2500
#define CHI2DIFF 1.

Int_t gfile,gtrig,gpt,gframe;

Float_t z1[40/NREBIN],z2[40/NREBIN],errorz1[40/NREBIN],errorz2[40/NREBIN],x[40/NREBIN],y[40/NREBIN];
Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht0","ht1","ht2"};

TFile* jpsifile;
TFile* crosscheckfile;
TFile* lambdas;

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME][2];// raw data, corrected data  // marker remove the last dimension
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH1F *costheta_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to costheta pass , total 
TH1F *phi_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to phi pass , total 
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH2F *chi2result[NFILE][NTRIG][NPT][NFRAME];

TH1F *costheta[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TH1F *phi[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TGraphAsymmErrors* efficiency1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];
TH1F * eff1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];

TGraphAsymmErrors *efficiency;
TFile *efficiencyfile[NFILE];
TFile *rawdatafile[NFILE][NTRIG];

TF1 *fit1[NPT][NFRAME][2]; //npt frame default/revised
TF1 *fit2[NPT][NFRAME][2];

double lambda_theta[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_theta_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_invariant[NFILE][NTRIG][NPT][NFRAME][3];// central value , statistic error , systematic error

double fit2dtheta[NFILE][NTRIG][NPT][NFRAME];
double fit2dtheta_err[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi_err[NFILE][NTRIG][NPT][NFRAME];
double fitparameters[NFILE][NTRIG][NPT][NFRAME][4];//theta , theta_err , phi , phi_err
double ptsmearing[NTRIG][NPT][NFRAME][NPHASE],weight[NTRIG][NPT][NFRAME][NPHASE],nhitsfit[NTRIG][NPT][NFRAME][NPHASE],nsigma[NTRIG][NPT][NFRAME][NPHASE],dsmadc[NTRIG][NPT][NFRAME][NPHASE],beta[NTRIG][NPT][NFRAME][NPHASE],poe[NTRIG][NPT][NFRAME][NPHASE],pol[NTRIG][NPT][NFRAME][NPHASE];
double systematic_error[NTRIG][NPT][NFRAME][NPHASE];

Double_t matrix[4][4];
Double_t covariant[NTRIG][NPT][NFRAME];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker
TCanvas *systematic[NTRIG];
TCanvas *check;

double fParamVal[NFILE][NTRIG][NPT][NFRAME][4];
double fParamErr[NFILE][NTRIG][NPT][NFRAME][4];
double arglist[10];
Int_t ierflg=0;
Int_t fitflag[4][NPT][2][4];
TH1F* fitdata[NTRIG][NPT][NFRAME][NPHASE];// theta , phi

TGraphErrors* lambda_theta_hx;// marker write them as an array
TGraphErrors* lambda_phi_hx; 
TGraphErrors* lambda_theta_cs;
TGraphErrors* lambda_phi_cs; 
TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda_parameters_sys[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant

TLegend* legend;

TFile *outputfile;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

TH2F* sample1;
TH2F* sample2;
TF2* sigma1[CHIXBIN][CHIXBIN];
TH2F* random1[CHIXBIN][CHIXBIN];
TF2* sigma2;
TH2F* random2;
TH2F* chi2;
TH2F* result_mean;
TH2F* result_error;
TF1* fitx;
TF1* fity;

TH2F* chiiii;

void toyMonteCarlo(int ifile = 0,int itrig = 0, int method = 0, int xsect = 1,int ysect = 1){
	gStyle->SetPadRightMargin(0.2);
	outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	
	gStyle->SetOptStat(false);
	correcteddata(ifile);
	fitx=new TF1("fitx","pol3",-1,1);
	fity=new TF1("fity","pol3",-1,1);

	result_mean = new TH2F("result_mean","",CHIXBIN,-1,1,CHIYBIN,-1,1);
	result_error = new TH2F("result_error","",CHIXBIN,0,1,CHIYBIN,0,1);

	//	static Double_t vstart[2];
	//	static Double_t step[2];
	TCanvas* canchiiii = new TCanvas("canvas","canvas",1200,1200);
	canchiiii->Divide(2,2);
	chiiii = new TH2F("chiiii","chiiii",CHIXBIN,-1,1,CHIYBIN,-1,1);

//	for(;itrig<NTRIG;itrig++){
		for(int ipt=0;ipt<NPT;ipt++){
			for(int iframe=0;iframe<NFRAME;iframe++){
				if(method==0)calchi2(ifile,itrig,ipt,iframe,xsect,ysect);
				if(method==1)calculatechi2(ifile,itrig,ipt,iframe);
			}
		}
//	}
	canchiiii->cd(1);
	chiiii->Draw("colz");
	canchiiii->SaveAs("~/polresults/20160707/figures/canchiii.pdf");

//	lambdaparameters(ifile,itrig);
}

void calculatechi2(int ifile,int itrig,int ipt,int iframe){
	static Double_t vstart[2];
	static Double_t step[2];

	gfile=ifile,gtrig=itrig,gpt=ipt,gframe=iframe;

	vstart[0] = 0.1;
	vstart[1] = 0.1;

	step[0] = 0.5;
	step[1] = 0.5;

	arglist[0] = 1;
	TMinuit *gMinuit = new TMinuit(3);
	gMinuit->SetFCN(fcn);
	gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

	gMinuit->mnparm(0,"lambda1",vstart[0],step[0],0,0,ierflg);
	gMinuit->mnparm(1,"lambda2",vstart[1],step[1],0,0,ierflg);

	arglist[0] = 500.;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	gMinuit->GetParameter(0,lambda_theta[ifile][itrig][ipt][iframe],lambda_theta_err[ifile][itrig][ipt][iframe]);
	gMinuit->GetParameter(1,lambda_phi[ifile][itrig][ipt][iframe],lambda_phi_err[ifile][itrig][ipt][iframe]);

}

void calchi2(int ifile,int itrig, int ipt,int iframe,int xsect,int ysect){
	if(rawdata2D[ifile][itrig][ipt][iframe][0]==0x0 || eff2D[ifile][itrig][ipt][iframe][2]==0x0) return;
	
	std::vector<double> chi2vector;
	chi2 = new TH2F("chi2","#chi^{2};#lambda_{#theta};#lambda_{#phi}",CHIXBIN,-1,1,CHIYBIN,-1,1);
	chi2->Sumw2();


	TRandom3 rndeff;
	sample2 = (TH2F*)rawdata2D[ifile][itrig][ipt][iframe][0]->Clone();
	sample2->Sumw2();
	sample2->Scale(1./sample2->Integral());

	TH2F* effhist = (TH2F*)eff2D[ifile][itrig][ipt][iframe][2]->Clone();
	effhist->SetTitle("efficiency;#phi;cos#theta");
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",ifile,itrig,ipt,iframe));

	TFile* lambdafile = new TFile(Form("~/polresults/20160707/rootfiles/file%d_trig%d_pt%d_frame%d_x%d_y%d.root",ifile,itrig,ipt,iframe,xsect,ysect),"recreate");
	lambdafile->cd();

	for(int xnbin=xsect;xnbin<xsect+10;xnbin++){
		for(int ynbin=ysect;ynbin<ysect+10;ynbin++){
			double lambda1 = chi2->GetXaxis()->GetBinCenter(xnbin);
			double lambda2 = chi2->GetYaxis()->GetBinCenter(ynbin);
			double chisquare = 0.;

			sigma1[xnbin-1][ynbin-1] = new TF2(Form("sigma_%d_%d",xnbin,ynbin),"1+[0]*y*y+[1]*(1-y*y)*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);
			sigma1[xnbin-1][ynbin-1]->SetTitle(Form("1+%.2f*cos^{2}#theta+%.2f*sin^{2}#thetacos(2#phi);#phi;cos#theta",lambda1,lambda2));
			sigma1[xnbin-1][ynbin-1]->SetParameters(lambda1,lambda2);
			sigma1[xnbin-1][ynbin-1]->Write();
			random2 = new TH2F(Form("theta%dphi%d_%d_%d_%d_%d",xnbin,ynbin,ifile,itrig,ipt,iframe),"template;#phi;cos#theta",YNBIN,-TMath::Pi(),TMath::Pi(),XNBIN,-1,1);
			random2->Sumw2();
			for(int i=0;i<1000*nsample;i++){
				double costheta_toy,phi_toy;
				sigma1[xnbin-1][ynbin-1]->GetRandom2(phi_toy,costheta_toy);
				double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi_toy),effhist->GetYaxis()->FindBin(costheta_toy));// marker
				random2->Fill(phi_toy,costheta_toy,eff);
			}
			random2->Scale(1./random2->Integral());
			random2->Write();
			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					if(sample2->GetBinContent(xbin,ybin)>=1e-3 && sample2->GetBinError(xbin,ybin)>=1e-6)chisquare += TMath::Power((random2->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin))/(sample2->GetBinError(xbin,ybin)),2);
				}
			}
			chi2->Fill(lambda1,lambda2,chisquare);
		}
	}
	int xx=1,yy=1,zz;
	chi2->GetMinimumBin(xx,yy,zz);
	result_mean->Fill(chi2->GetXaxis()->GetBinCenter(xx),chi2->GetYaxis()->GetBinCenter(yy));
	int xx1=1,xx2=CHIXBIN,yy1=1,yy2=CHIYBIN;
	double chi2min,deltachi2;

/*
	chi2->ProjectionX("xx",yy,yy)->Fit(fitx);
	chi2min=fitx->Eval(chi2->GetXaxis()->GetBinCenter(xx));
	deltachi2=100;

	for(int ix=1;ix<=xx;ix++) {
		double chi2scan=fitx->Eval(-1+2./CHIXBIN*ix);
		if(fabs(chi2scan-chi2min-CHI2DIFF)<deltachi2) {
			xx1=ix;
			deltachi2=fabs(chi2scan-chi2min-CHI2DIFF);
		}
	}
	deltachi2=100;
	for(int ix=CHIXBIN;ix>=xx;ix--) {
		double chi2scan=fitx->Eval(-1+2./CHIXBIN*ix);
		if(fabs(chi2scan-chi2min-CHI2DIFF)<deltachi2) {
			xx2=ix;
			deltachi2=fabs(chi2scan-chi2min-CHI2DIFF);
		}
	}

	chi2->ProjectionY("yy",xx,xx)->Fit(fity);
	chi2min=fity->Eval(chi2->GetYaxis()->GetBinCenter(yy));
	deltachi2=100;
	for(int iy=1;iy<=yy;iy++) {
		double chi2scan=fity->Eval(-1+2./CHIYBIN*iy);
		if(fabs(chi2scan-chi2min-CHI2DIFF)<deltachi2) {
			yy1=iy;
			deltachi2=fabs(chi2scan-chi2min-CHI2DIFF);
		}
	}
	deltachi2=100;
	for(int iy=CHIYBIN;iy>=yy;iy--) {
		double chi2scan=fity->Eval(-1+2./CHIYBIN*iy);
		if(fabs(chi2scan-chi2min-CHI2DIFF)<deltachi2) {
			yy2=iy;
			deltachi2=fabs(chi2scan-chi2min-CHI2DIFF);
		}
	}

	cout<<chi2->GetXaxis()->GetBinCenter(xx)<<" "<<chi2->GetYaxis()->GetBinCenter(yy)<<" "<<chi2->GetBinContent(xx,yy)<<endl;
	cout<<chi2->GetXaxis()->GetBinCenter(xx1)<<" "<<chi2->GetYaxis()->GetBinCenter(yy)<<" "<<chi2->GetBinContent(xx1,yy)<<endl;
	cout<<chi2->GetXaxis()->GetBinCenter(xx2)<<" "<<chi2->GetYaxis()->GetBinCenter(yy)<<" "<<chi2->GetBinContent(xx2,yy)<<endl;
	cout<<chi2->GetXaxis()->GetBinCenter(xx)<<" "<<chi2->GetYaxis()->GetBinCenter(yy1)<<" "<<chi2->GetBinContent(xx,yy1)<<endl;
	cout<<chi2->GetXaxis()->GetBinCenter(xx)<<" "<<chi2->GetYaxis()->GetBinCenter(yy2)<<" "<<chi2->GetBinContent(xx,yy2)<<endl;
	result_error->Fill(0.5*(chi2->GetXaxis()->GetBinCenter(xx2)-chi2->GetXaxis()->GetBinCenter(xx1)),0.5*(chi2->GetYaxis()->GetBinCenter(yy2)-chi2->GetYaxis()->GetBinCenter(yy1)));
*/
//	lambda_theta[0][itrig][ipt][iframe] = chi2->GetXaxis()->GetBinCenter(xx);
//	lambda_theta_err[0][itrig][ipt][iframe] = chi2->GetXaxis()->GetBinCenter(xx2)-chi2->GetXaxis()->GetBinCenter(xx1);
//	lambda_phi[0][itrig][ipt][iframe] = chi2->GetYaxis()->GetBinCenter(yy);
//	lambda_phi_err[0][itrig][ipt][iframe] = chi2->GetYaxis()->GetBinCenter(yy2)-chi2->GetYaxis()->GetBinCenter(yy1);

	check = new TCanvas("check","check",1200,1200);
	check->Divide(3,3);

	check->cd(1);        
		cout<<"xx==============="<<xx<<"yy===================="<<yy<<endl;
	if(xx<=xsect-1) xx=xsect;
	if(yy<=ysect-1) yy=ysect;
	if(xx>=10+xsect) xx=9+xsect;
	if(yy>=10+ysect) yy=9+ysect;
//	   	&& yy>ysect-1 && xx<10+xsect && yy<10+ysect) 
		sigma1[xx-1][yy-1]->Draw("colz");

	check->cd(2);
	effhist->Draw("colz");

	check->cd(3);
	random2->Draw("colz");

	check->cd(4);
	random2->Draw("colz");

	check->cd(5);
	sample2->Draw("colz");

	check->cd(6);
	chi2->Draw("colz");

	check->cd(7);
//	chi2->ProjectionX("pro_x",yy,yy)->Fit(fitx);

	check->cd(8);	
//	chi2->ProjectionY("pro_y",xx,xx)->Fit(fity);

	check->SaveAs(Form("~/polresults/20160707/figures/check_%d_trig%d_pt%d_frame%d_x%d_y%d.pdf",ifile,itrig,ipt,iframe,xsect,ysect));
	//	check->Clear();
	sigma1[xx-1][yy-1]->Write();
	effhist->Write();
	sample2->Write();
	result_mean->Write();
	chi2->Write();
	lambdafile->Close();

	for(int i=xsect-1;i<xsect+9;i++){
		for(int j=ysect-1;j<ysect+9;j++)sigma1[i][j]->Clear();
	}

	sample2->Clear();
	chi2->Clear();
}

void correcteddata(int file = 0){
	int selectedfile,file;
	if(file>=0 && file<=NFILE) selectedfile=file+1;
	else if(file==100) file=0,selectedfile=NFILE;
	else return;

	for(;file<selectedfile;file++){
		efficiencyfile[file] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/OutFile_cent_0_9_%d.root",file),"read"); 
		for(int trig=0;trig<NTRIG;trig++){//marker
			if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,file),"read");//marker
			else rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,0),"read");//marker
			for(int pt=0;pt<NPT;pt++){
				for(int frame=0;frame<NFRAME;frame++){
					canvas[trig][pt][frame] = new TCanvas(Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),3200,1600);
					canvas[trig][pt][frame]->Divide(4,3);
					canvas[trig][pt][frame]->cd(1);
					if(frame==0){
						rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPt"));
						rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtBG"));
						rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_hx");
						rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_hx");
					}
					else{
						rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtCS"));
						rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtCSBG"));
						rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_cs");
						rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_cs");
					}
					rawdata3D[file][trig][pt][frame][0]->Sumw2();
					rawdata3D[file][trig][pt][frame][0]->Draw();

					canvas[trig][pt][frame]->cd(2);
					rawdata3D[file][trig][pt][frame][1]->Sumw2();
					rawdata3D[file][trig][pt][frame][1]->Draw();

					canvas[trig][pt][frame]->cd(3);
					rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike-like");
					if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_hx");
					else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_cs");
					rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);
					rawdata3D[file][trig][pt][frame][2]->Sumw2();
					rawdata3D[file][trig][pt][frame][2]->Draw();

					int max = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
					int min = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	

					if(frame==0){
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPt1");
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
					}
					else{	
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPtCS1");
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
					}

					rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);
					rawdata3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					rawdata3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));

					canvas[trig][pt][frame]->cd(4);
					rawdata2D[file][trig][pt][frame][0] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
					rawdata2D[file][trig][pt][frame][0]->SetName(Form("rawdata_%s_pT%d_%d_frame%d",trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame));
					rawdata2D[file][trig][pt][frame][0]->RebinX(NREBIN)->RebinY(NREBIN);	
					rawdata2D[file][trig][pt][frame][0]->Draw("colz");

					canvas[trig][pt][frame]->cd(5);
					eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
					eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);

					eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
					eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);

					eff2D[file][trig][pt][frame][2] = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
					eff2D[file][trig][pt][frame][2]->SetName(Form("eff_ratio_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][2]->Divide(eff2D[file][trig][pt][frame][1]);
					eff2D[file][trig][pt][frame][2]->Draw("colz");

//					if(pt==0) canvas[trig][pt][frame]->Print(Form("~/polresults/20160707/figures/%s_frame%d.pdf(",trigName[trig].Data(),frame));			
//					else if(pt==NPT-1) canvas[trig][pt][frame]->Print(Form("~/polresults/20160707/figures/%s_frame%d.pdf)",trigName[trig].Data(),frame));			
//					else canvas[trig][pt][frame]->Print(Form("~/polresults/20160707/figures/%s_frame%d.pdf",trigName[trig].Data(),frame));			
				}
			}
		}
	}
}

void compare_lambda(int trig,int frame, int phase){
	if(phase==0)	comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
	if(phase==1) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));
	if(phase==2) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_inv_%d_1eID",trigName[trig].Data(),frame));

	if(trig==0){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kRed);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kRed);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(24);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==1){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlue);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlue);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(26);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==2){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlack);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlack);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(27);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	//	if(phase!=2) 
	comparisonlambda[trig][phase+3*frame]->Draw("samep");
}

void crosscheck(){
	TFile *oldfile = new TFile("~/jpsi/test20160210_Barbara/efficiency1eID.root","read");
	TH1F *oldhist[NTRIG][NPT];
	TH1F *h1[6],*h2[NPT][3];

	TFile *defaultfile[3];// embedding ht0_1eID ht2_1eID
	TH1F *defaulthist[2],*defaultrawdata[3];

	TCanvas *canvas2 = new TCanvas("canvas2","canvas2",1200,800);
	canvas2->Divide(6,2);

	gStyle->SetOptStat(0);
	defaultfile[0] = new TFile("~/jpsi/test20160210_0/rootfile/OutFile_cent_0_9.root","read");
	defaultfile[1] = new TFile("~/jpsi/test20160210_0/rootfile/ht0_trg3_1TrkPid_0.ana.root","read");
	defaultfile[2] = new TFile("~/jpsi/test20160210_0/rootfile/ht2_trg5_1TrkPid_0.ana.root","read");

	int xmin,xmax,ymin,ymax,zmin,zmax;

	for(int trig=0;trig<1;trig++){
		for(int pt=0;pt<NPT;pt++){
			if((TH1F*)oldfile->Get(Form("hEff%sHXThetaPt%d",trigName[trig].Data(),pt))!=0x0) oldhist[trig][pt] = (TH1F*)oldfile->Get(Form("hEff%sHXThetaPt%d",trigName[trig].Data(),pt));	
			ymin = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(-TMath::Pi());
			ymax = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(TMath::Pi());
			zmin = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt]);//marker
			zmax = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt+1])-1;

			canvas2->cd(pt+1);
			defaulthist[0] = (TH1F*)((TH2F*)defaultfile[0]->Get(Form("hHt%dJpsiCosThetaPt1",trig)))->ProjectionX(Form("defaulthist_%d_%d",trig,pt),zmin,zmax);
			defaulthist[0]->SetTitle("comparison of pass between 3D projection and 2D projection");
			defaulthist[0]->GetXaxis()->SetTitle("cos#theta");

			defaulthist[1] = (TH1F*)((TH2F*)defaultfile[0]->Get("hMcJpsiCosThetaPt"))->ProjectionX(Form("Mcdefaulthist_%d_%d",trig,pt),zmin,zmax);
			defaulthist[1]->SetTitle("comparison of total between 3D projection and 2D projection");
			defaulthist[1]->GetXaxis()->SetTitle("cos#theta");
			h1[0] = (TH1F*)eff3D[0][0][1][0][0]->ProjectionX("h1",ymin,ymax,zmin,zmax);
			h1[0]->RebinX(4);
			h1[0]->SetMarkerStyle(20);

			h1[1] = (TH1F*)eff3D[0][1][1][0][1]->ProjectionX("h2",ymin,ymax,zmin,zmax);
			h1[1]->RebinX(4);
			h1[1]->SetMarkerStyle(20);

			TGraphAsymmErrors* test = new TGraphAsymmErrors(h1[0],h1[1],"N");
			oldhist[trig][pt]->SetTitle("efficiency comparison");
			oldhist[trig][pt]->GetXaxis()->SetTitle("cos#theta");
			oldhist[trig][pt]->Draw("p");
			test->Draw("samep");

			canvas2->cd(7+pt);
			defaultrawdata[0] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPt"))->ProjectionX("rawdata",zmin,zmax);
			defaultrawdata[0]->SetTitle("Unlike Sign Comparison");
			h2[pt][0] = (TH1F*)rawdata3D[0][0][2][0][0]->ProjectionX("h4",ymin,ymax,zmin,zmax);
			h2[pt][0]->RebinX(4);
			h2[pt][0]->SetMarkerStyle(24);

			defaultrawdata[1] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPtBG"))->ProjectionX("rawdatabkg",1,8);
			defaultrawdata[1]->SetTitle("Like Sign Comparison");
			h2[pt][1] = (TH1F*)rawdata3D[0][0][2][0][1]->ProjectionX("h5",1,41,1,8);
			h2[pt][1]->RebinX(4);
			h2[pt][1]->SetMarkerStyle(24);
			h2[pt][2] = (TH1F*)h2[pt][0]->Clone();
			h2[pt][2]->Add(h2[pt][1],-1);
			h2[pt][2]->SetTitle("raw data");
			h2[pt][2]->Draw();

			defaultrawdata[2] = (TH1F*)defaultrawdata[0]->Clone();
			defaultrawdata[2]->Add(defaultrawdata[1],-1);	
			defaultrawdata[2]->Draw("same");
		}
	}
}

double getaverage(std::vector<double> inputvector){
	double sum = 0.;
	for(int i=0;i<inputvector.size();i++) sum += inputvector[i];
	sum = sum/inputvector.size();
	return sum;
}

double getmaximum(std::vector<double> inputvector){
	double max;
	max = *std::max_element(&inputvector[0],&inputvector[0]+(inputvector.size()));// marker
	return max;
}

double getminimum(std::vector<double> inputvector){
	double min;
	min = *std::min_element(&inputvector[0],&inputvector[0]+(inputvector.size()));
	return min;
}

void lambdaparameters(int sys3,int trig){// plot and write lambda parameters 

	int selectedfile,file;
	if(sys3>=0 && sys3<=NFILE) file=sys3,selectedfile=sys3+1;
	else if(sys3==100) {
		file=0,selectedfile=NFILE;
	}
	else return;

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT],y_sys[NPT],theta[NPT],thetaerr[NPT],phi[NPT],phierr[NPT];

	for(;file<selectedfile;file++){
//		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){ 
				for(int phase=0;phase<NPHASE+1;phase++){
					for(int pt=0;pt<NPT;pt++) {
						x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.+0.2*trig;
						x_err[pt] = (PtEdge[pt+1]-PtEdge[pt])/2.;
						theta[pt] = lambda_theta[0][trig][pt][frame];
						thetaerr[pt] = lambda_theta_err[0][trig][pt][frame];
						phi[pt] = lambda_phi[0][trig][pt][frame];
						phierr[pt] = lambda_phi_err[0][trig][pt][frame];

						if(phase==0){
							y[pt] = lambda_theta[0][trig][pt][frame];
							y_err[pt] = lambda_theta_err[0][trig][pt][frame];
						}
						else if (phase==1){
							y[pt] = lambda_phi[0][trig][pt][frame];
							y_err[pt] = lambda_phi_err[0][trig][pt][frame];
						}
						else if (phase==2){//calculate lambda_invariant and its statistic error
							y[pt] = (theta[pt]+3*phi[pt])/(1-phi[pt]);
							y_err[pt] = TMath::Sqrt(theta[pt]*theta[pt]/((1-phi[pt])*(1-phi[pt]))+phi[pt]*phi[pt]*TMath::Power((3+theta[pt]/((1-phi[pt])*(1-phi[pt]))),2)+2*(3+theta[pt])/TMath::Power((1-phi[pt]),3)*covariant[trig][pt][frame]);
							//							y_sys[pt] = y[pt]-(fitparameters[0][trig][pt][frame][0]+3*fitparameters[0][trig][pt][frame][2])/(1-fitparameters[0][trig][pt][frame][2]);
							cout<<"invariant ==================="<<y[pt]<<"trig"<<trig<<"   "<<"frame"<<frame<<endl;
						}
					}
					lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
					lambda_parameters[trig][frame][phase]->SetName(Form("lambdas_%d_%d_%d",trig,frame,phase));
					lambda_parameters_sys[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_sys);
					lambda_parameters_sys[trig][frame][phase]->SetName(Form("lambdas_sys_%d_%d_%d",trig,frame,phase));
					if(sys3==100) {
						//						lambda_parameters[trig][frame][phase]->Write();		
						//			lambda_parameters_sys[trig][frame][phase]->Write();		
					}
				}
			}
//		}
	}
	drawlambdas(0,trig);
}

void drawlambdas(int drawoptions = 0,int trig){
	TCanvas* lambdacanvas; 
	lambdas = new TFile("~/polresults/20160707/rootfiles/lambdas.root","recreate");
	lambdas->cd();
	if(drawoptions==1){
		lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
		lambdacanvas->Divide(3,2);
//		for(int trig=0;trig<NTRIG;trig++){
			for(int phase=0;phase<NPHASE+1;phase++){
				for(int frame=0;frame<NFRAME;frame++){ 
					lambdacanvas->cd(frame*3+phase+1);
					legend = new TLegend(0.7,0.7,0.89,0.89);
					lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
					if(trig==0){
						lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
						lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
						lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
						lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
						lambda_parameters[trig][frame][phase]->SetMarkerStyle(20);
					}
					else if(trig==1){
						lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
						lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
						lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
						lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
						lambda_parameters[trig][frame][phase]->SetMarkerStyle(22);
					}
					else lambda_parameters[trig][frame][phase]->SetMarkerStyle(23);
					plotaxissetting();
					lambda_parameters[trig][frame][phase]->Draw("ap");
					lambda_parameters[trig][frame][phase]->Write();
					compare_lambda(trig,frame,phase);
					legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
					legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
					legend->Draw("same");
				}
			}
//		}
		lambdacanvas->SaveAs("~/polresults/20160707/figures/lambdas_comparison_all.pdf");
	}

//	for(int trig=0;trig<NTRIG;trig++){
		lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
		lambdacanvas->Divide(3,2);
		for(int phase=0;phase<NPHASE+1;phase++){
			for(int frame=0;frame<NFRAME;frame++){ 
				lambdacanvas->cd(frame*3+phase+1);
				legend = new TLegend(0.7,0.7,0.89,0.89);
				lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
				if(trig==0){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(20);
				}
				else if(trig==1){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(22);
				}
				else lambda_parameters[trig][frame][phase]->SetMarkerStyle(23);
				plotaxissetting();
				lambda_parameters[trig][frame][phase]->Draw("ap");
				lambda_parameters[trig][frame][phase]->Write();
				//				if(phase!=2)
				compare_lambda(trig,frame,phase);
				legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
				legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
				legend->Draw("same");
				//				lambda_parameters_sys[trig][frame][phase]->Draw("samep[]");
			}
		}
		lambdacanvas->SaveAs(Form("~/polresults/20160707/figures/lambdas_comparison_%d.pdf",trig));
//	}
}

void plotaxissetting(){
	lambda_parameters[0][0][0]->SetTitle("#lambda_{#theta} in HX frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][1][0]->SetTitle("#lambda_{#theta} in CS frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][0][1]->SetTitle("#lambda_{#phi} in HX frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][1][1]->SetTitle("#lambda_{#phi} in CS frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][0][2]->SetTitle("#lambda_{inv} in HX frame;J/#psi p_{T};#lambda_{inv}");
	lambda_parameters[0][1][2]->SetTitle("#lambda_{inv} in CS frame;J/#psi p_{T};#lambda_{inv}");
}

void simultaneousfit(int sys1 = 0){// marker add frame 
	TCanvas* chi2resultplot[NFILE][NTRIG][NPT][NFRAME];
	int file,selectedfile;
	if(sys1>=0 && sys1<=NFILE) file=sys1,selectedfile=sys1+1;
	else if(sys1==100) file=0,selectedfile=NFILE;
	for(;file<selectedfile;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){
				for(int pt=0;pt<NPT;pt++){
					chi2resultplot[file][trig][pt][frame] = new TCanvas(Form("chi2result_%d_%d_%d_%d",file,trig,pt,frame));
					chi2result[file][trig][pt][frame] = new TH2F(Form("chi2result_%d_%d_%d_%d",file,trig,pt,frame),"#chi^{2}; #lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1);
					pull(file,trig,pt,frame);
					calchi2(file,trig,pt,frame);
					chi2result[file][trig][pt][frame]->Fill(lambda_theta[file][trig][pt][frame],lambda_phi[file][trig][pt][frame]);
					chi2resultplot[file][trig][pt][frame]->cd();
					chi2result[file][trig][pt][frame]->Draw("colz");
					chi2resultplot[file][trig][pt][frame]->SaveAs(Form("~/polresults/20160707/figures/chi2resultplot_%d_%d_%d_%d.pdf",file,trig,pt,frame));
				}
			}
		}
	}
}

void systematic(int sys = 0){

	correcteddata(0);
	simultaneousfit(0);

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT];

	systematic[0] = new TCanvas("systematic","systematic",1200,800);
	systematic[0]->Divide(4,2);
	for(int i=0;i<NPT;i++){
		x[i] = (PtEdge[i]+PtEdge[i+1])/2.;
		x_err[i] = (PtEdge[i+1]-PtEdge[i])/2.;
		y[i] = fit2dtheta[0][0][i][0];
		y_err[i] = fit2dtheta_err[0][0][i][0];
		cout<<"lambda_theta======"<<y[i]<<endl;	
	}

	lambda_theta_hx =  new TGraphErrors(NPT,x,y,x_err,y_err);
	systematic[0]->cd(1);
	lambda_theta_hx->GetYaxis()->SetRangeUser(-2,2);
	lambda_theta_hx->Draw("ape");
	//	systematic[0]->SaveAs("systematic1.pdf");	

	for(int i=0;i<NPT;i++){
		x[i] = (PtEdge[i]+PtEdge[i+1])/2.;
		x_err[i] = (PtEdge[i+1]-PtEdge[i])/2.;
		y[i] = fit2dphi[0][0][i][0];
		y_err[i] = fit2dphi_err[0][0][i][0];

		//		y[i] = lambda_phi[0][0][i][0];
		//		y_err[i] = lambda_phi_err[0][0][i][0];	
		cout<<"lambda_phi======"<<y[i]<<endl;	
	}

	lambda_phi_hx =  new TGraphErrors(NPT,x,y,x_err,y_err);
	systematic[0]->cd(2);
	lambda_phi_hx->GetYaxis()->SetRangeUser(-2,2);
	lambda_phi_hx->Draw("ape");
	systematic[0]->SaveAs("~/polresults/20160707/figures/systematic1.pdf");	

	for(int file=0;file<NFILE;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){
				for(int phase=0;phase<NPHASE;phase++){
					for(int pt=0;pt<NPT;pt++) {
						y[pt] = fitparameters[file][trig][pt][frame][phase*2];
						y_err[pt] = fitparameters[file][trig][pt][frame][phase*2+1];
					}
					lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);
				}
			}
		}
	}

	if(sys>=1 && sys<=NFILE){
		correcteddata(sys);
		simultaneousfit(sys);
		int trig=0;
		for(int pt=0;pt<5;pt++) cout<<"lambda_theta in "<<sys<<"   "<<trig<<"  pt "<<pt<<"    "<<lambda_theta[sys][trig][pt][0]<<"&&&&&&&&&"<<lambda_theta[0][trig][pt][0]<<endl;
		//		cout<<"systematic uncertainty from "<<sys<<endl;//marker
	}
	else if(sys==100){
		//		for(sys=1;sys<NFILE;sys++){
		correcteddata(sys);
		simultaneousfit(sys);
		//		}
		std::vector<double> sysvector_theta_hx,sysvector_theta_cs,sysvector_phi_hx,sysvector_phi_cs;

		int file,trig,pt,frame;

		for(trig=0;trig<NTRIG;trig++){//marker add the definition of systematic 
			for(pt=0;pt<NPT;pt++){
				for(file=1;file<12;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);// revise the frame 
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);// revise the frame 
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				ptsmearing[trig][pt][0][0] = getaverage(sysvector_theta_hx);// marker getmaximum()
				ptsmearing[trig][pt][1][0] = getaverage(sysvector_theta_cs);// marker getmaximum()
				ptsmearing[trig][pt][0][1] = getaverage(sysvector_phi_hx);// marker getmaximum()
				ptsmearing[trig][pt][1][1] = getaverage(sysvector_phi_cs);// marker getmaximum()
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"ptsmearing = "<<ptsmearing[trig][pt][0][0]<<endl;

				for(file=12;file<20;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				weight[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				weight[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				weight[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				weight[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"weight = "<<weight[trig][pt][0][0]<<endl;

				for(file=20;file<24;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				nhitsfit[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
				nhitsfit[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
				nhitsfit[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
				nhitsfit[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"nhitsfit = "<<nhitsfit[trig][pt][0][0]<<endl;

				for(file=24;file<29;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				nsigma[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				nsigma[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				nsigma[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				nsigma[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"nsigma = "<<nsigma[trig][pt][0][0]<<endl;

				for(file=29;file<31;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				dsmadc[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				dsmadc[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				dsmadc[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				dsmadc[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"dsmadc = "<<dsmadc[trig][pt][0][0]<<endl;

				for(file=31;file<35;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				beta[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
				beta[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
				beta[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
				beta[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"beta = "<<beta[trig][pt][0][0]<<endl;

				for(file=35;file<36;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				poe[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				poe[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				poe[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"poe = "<<poe[trig][pt][0][0]<<endl;

				for(file=36;file<40;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				pol[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				pol[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				pol[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				pol[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"pol = "<<pol[trig][pt][0][0]<<endl;

				cout<<"combining systematic uncertainty"<<endl;//marker

				for(int frame=0;frame<NFRAME;frame++){		
					for(int phase=0;phase<NPHASE;phase++){
						systematic_error[trig][pt][frame][phase] = TMath::Sqrt(fabs(ptsmearing[trig][pt][frame][phase])*fabs(ptsmearing[trig][pt][frame][phase])+fabs(weight[trig][pt][frame][phase])*fabs(weight[trig][pt][frame][phase])+fabs(nhitsfit[trig][pt][frame][phase])*fabs(nhitsfit[trig][pt][frame][phase])+fabs(nsigma[trig][pt][frame][phase])*fabs(nsigma[trig][pt][frame][phase])+fabs(dsmadc[trig][pt][frame][phase])*fabs(dsmadc[trig][pt][frame][phase])+fabs(beta[trig][pt][frame][phase])*fabs(beta[trig][pt][frame][phase])+fabs(poe[trig][pt][frame][phase])*fabs(poe[trig][pt][frame][phase])+fabs(pol[trig][pt][frame][phase])*fabs(pol[trig][pt][frame][phase]));
					}
				}
				cout<<"combining systematic uncertainty is done"<<endl;
			}
		}
	}
}

double RosenBrock(const double *xx )
{
	if(!chi2) {
		TFile* outputfile = new TFile("outputfile.root","read");
		chi2=(TH2F*)outputfile->Get("chi2");
	}
	int xbin = chi2->GetXaxis()->FindBin(xx[0]-1e-5);
	int ybin = chi2->GetYaxis()->FindBin(xx[1]-1e-5);

	if(fabs(xx[0])<1&&fabs(xx[1])<1) return chi2->GetBinContent(xbin,ybin);
	return 1e8;
}

double RosenBrock2(const double *xx)
{
	if(fabs(xx[0])>1||fabs(xx[1])>1) return 1e8;

	TRandom3 rndeff;
	double lambda_theta = xx[0];
	double lambda_phi = xx[1];
	sigma2 = new TF2("sigma2","1+[0]*x*x+[1]*(1-x*x)*TMath::Cos(2*y)",-1,1,-TMath::Pi(),TMath::Pi());
	sigma2->SetParameters(lambda_theta,lambda_phi);

	random2=new TH2F("random2","random2",XNBIN,-1,1,YNBIN,-TMath::Pi(),TMath::Pi());
	random2->Sumw2();

	if(eff2D[gfile][gtrig][gpt][gframe][2]==0x0 || rawdata2D[gfile][gtrig][gpt][gframe][0]==0x0) return 1e8;
	TH2F* effhist = (TH2F*)eff2D[gfile][gtrig][gpt][gframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",gfile,gtrig,gpt,gframe));

	for(int i=0;i<100*NEVENT;i++){
		double costheta,phi;
		sigma2->GetRandom2(costheta,phi);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi),effhist->GetYaxis()->FindBin(costheta));
		if(rndeff.Uniform()<eff) random2->Fill(costheta,phi);
	}
	random2->Scale(1./random2->Integral());

	double chisquare = 0.;
	for(int xbin=1;xbin<XNBIN+1;xbin++){
		for(int ybin=1;ybin<YNBIN+1;ybin++){
			chisquare += TMath::Power((random2->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin))/(sample2->GetBinError(xbin,ybin)),2);
		}
	}

	sigma2->Clear();
	random2->Clear();
	return chisquare;
}

double fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	//	TCanvas* canvas1 = new TCanvas("canvas1","canvas1",1200,1200);
	//	canvas1->Divide(2,2);

	double chisquare = 0.;
	double lambda1 = par[0];
	double lambda2 = par[1];

	if(rawdata2D[gfile][gtrig][gpt][gframe][0]==0x0 || eff2D[gfile][gtrig][gpt][gframe][2]==0x0 ) {
		cout<<"something is wrong"<<endl;
		f = chisquare;	
		return 0.;
	}

	TRandom3 rndeff;
	sample2 = (TH2F*)rawdata2D[gfile][gtrig][gpt][gframe][0]->Clone();
	sample2->Scale(1./sample2->Integral());

	TH2F* effhist = (TH2F*)eff2D[gfile][gtrig][gpt][gframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",gfile,gtrig,gpt,gframe));

	sigma2 = new TF2(Form("sigma2"),"1+[0]*y*y+[1]*(1-y*y)*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);
	sigma2->SetParameters(lambda1,lambda2);
	random2 = new TH2F(Form("thetaphi_%d_%d_%d_%d",gfile,gtrig,gpt,gframe),"random2",YNBIN,-TMath::Pi(),TMath::Pi(),XNBIN,-1,1);
	random2->Sumw2();
	for(int i=0;i<1000*nsample;i++){
		double costheta_toy,phi_toy;
		sigma2->GetRandom2(phi_toy,costheta_toy);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi_toy),effhist->GetYaxis()->FindBin(costheta_toy));// marker
		if(rndeff.Uniform()<eff) random2->Fill(phi_toy,costheta_toy);
	}
	random2->Scale(1./random2->Integral());

	for(int xbin=1;xbin<XNBIN+1;xbin++){
		for(int ybin=1;ybin<YNBIN+1;ybin++){
			if(sample2->GetBinContent(xbin,ybin)>=1e-6 && sample2->GetBinError(xbin,ybin)>=1e-6)chisquare += TMath::Power((random2->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin))/(sample2->GetBinError(xbin,ybin)),2);

			cout<<"sample2============="<<(sample2->GetBinContent(xbin,ybin))<<endl;
			cout<<"sample2err============="<<(sample2->GetBinError(xbin,ybin))<<endl;
			cout<<"random22============="<<(random2->GetBinContent(xbin,ybin))<<endl;
		}
	}

	//	canvas1->cd(1);
	//	sample2->Draw("colz");
	//	canvas1->cd(2);
	//	effhist->Draw("colz");
	//	canvas1->cd(3);
	//	random2->Draw("colz");

	//	canvas1->SaveAs("figures/canvas1.pdf");
	//	canvas1->Clear();	

	//	chiiii->Fill(lambda1,lambda2,chisquare);

	sigma2->Clear();
	random2->Clear();
	f = chisquare;

	cout<<"chisquare===================="<<chisquare<<endl;
	cout<<"ff===========================>"<<f<<endl;
	return chisquare;
}
