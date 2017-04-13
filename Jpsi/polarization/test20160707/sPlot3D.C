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

#define NEXP     1
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
#define ZMIN -1
#define ZMAX 1
#define NBIN    10

Int_t NSIGNAL = 500;
Int_t NBKG = 200;

//TF2* fsig=new TF2("fsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
//TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*cos(x)*y*(1-y**2)**0.5)",XMIN,XMAX,YMIN,YMAX);
TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
//TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
//TF3* lsig=new TF3("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);

TH2F *hunl, *hlik, *hsig, *hbkg, *hsub;
TH3F *hlambda, *hlambda2, *hlambda_err, *hlambda2_err;
TGraph *gc, *gc2;

TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);

#define NFILE 100
#define NTRIG 3
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

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht0","ht1","ht2"};

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

TH1F* st;
TH1F* method;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	f = likelihoodfcn(par);
}

float likelihoodfcn(Double_t *x){
	lsig->SetParameter(0,1);
	lsig->SetParameter(1,x[0]);
	lsig->SetParameter(2,x[1]);
	lsig->SetParameter(3,x[2]);

	double phi,costheta;
	double scale=0.;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			scale += effhist->GetBinContent(i,j)*lsig->Eval(phi,costheta);
		}
	}
	lsig->SetParameter(0,1./scale); //renormalization	

	double eff;
	float result=0;
	double sum=0;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			if(effhist->GetBinContent(i,j)>1e-6) eff = effhist->GetBinContent(i,j);
		//	result += -1*hsub->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)));
			result += -1*datahist->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)));
			sum += eff*lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)); // check the normalization
		}
	}
	result*=datahist->Integral()/(tothist->Integral()+bkghist->Integral());
	return result;
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

Int_t Minimize(int n=0,int trig,int pt,int frame){
//	hunl=new TH2F("hunl","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
//	hlik=new TH2F("hlik","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
//	hsig=new TH2F("hsig","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
//	hbkg=new TH2F("hbkg","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
//	hsub=new TH2F("hsub","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);

	gRandom = new TRandom3();
	gRandom->SetSeed();

	NSIGNAL = datahist->Integral();	
	NBKG = bkghist->Integral();
/*
	for(int i=0;i<gRandom->Poisson(NSIGNAL);) {
		Double_t xsig,ysig;
		fsig->GetRandom2(xsig,ysig);
		if(gRandom->Uniform()<effhist->GetBinContent(effhist->FindBin(xsig,ysig))){
			hsig->Fill(xsig,ysig);
			hunl->Fill(xsig,ysig);
			hsub->Fill(xsig,ysig);
			i++;
		}
	}
	for(int i=0;i<gRandom->Poisson(NBKG);) {
		Double_t xbkg,ybkg;
		fbkg->GetRandom2(xbkg,ybkg);
		if(gRandom->Uniform()<effhist->GetBinContent(effhist->FindBin(xbkg,ybkg))){
			hbkg->Fill(xbkg,ybkg);
			hunl->Fill(xbkg,ybkg);
			hsub->Fill(xbkg,ybkg);
			i++;
		}
	}
	for(int i=0;i<gRandom->Poisson(NBKG);) {
		Double_t xbkg,ybkg;
		fbkg->GetRandom2(xbkg,ybkg);
		if(gRandom->Uniform()<effhist->GetBinContent(effhist->FindBin(xbkg,ybkg))){
			hlik->Fill(xbkg,ybkg);
			i++;
		}
	}
	hsub->Add(hlik,-1);
*/

	TMinuit *gMinuit = new TMinuit(3);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,10,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[4] = {0.01, 0.01 , 0.01 , 0.01};
	static Double_t step[4] = {0.01 , 0.01 , 0.01 , 0.001};
	gMinuit->mnparm(0, "p1", vstart[0], step[0], -1,1,ierflg);
	gMinuit->mnparm(1, "p2", vstart[1], step[1], -1,1,ierflg);
	gMinuit->mnparm(2, "p3", vstart[2], step[2], -1,1,ierflg);

	Double_t scan1list[4]={0,1000,-1,1};
	Double_t scan2list[4]={1,1000,-1,1};
	Double_t scan3list[4]={2,1000,-1,1};
	gMinuit->mnexcm("SCAN",scan1list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan2list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan3list,4,ierflg);	

	// Now ready for minimization step
	arglist[0] = 5000;
	arglist[1] = 10.;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	gMinuit->mnexcm("HESSE",arglist,2,ierflg);
	Double_t minoslist[4] = {5000,0,1,2};

	gMinuit->mnexcm("MINOS",minoslist,4,ierflg);

	Double_t getlambdatheta,getlambdathetaerr,getlambdaphi,getlambdaphierr,getlambdathetaphi,getlambdathetaphierr;

	gMinuit->GetParameter(0,getlambdatheta,getlambdathetaerr);
	gMinuit->GetParameter(1,getlambdaphi,getlambdaphierr);
	gMinuit->GetParameter(2,getlambdathetaphi,getlambdathetaphierr);

	st->Fill(ierflg);

	cout<<"status ====="<<st<<"=============="<<ierflg<<endl;

	cout<<"lambdatheta    = "<<getlambdatheta<<"  +/-   "<<getlambdathetaerr<<endl;
	cout<<"lambdaphi      = "<<getlambdaphi<<"  +/-  "<<getlambdaphierr<<endl;
	cout<<"lambdathetaphi = "<<getlambdathetaphi<<"  +/-  "<<getlambdathetaphierr<<endl;

	double lbthesig,lbphisig,lbthetaphi,lbthesig_err,lbphisig_err,lbthetaphi_err;
	double deltacontour = 0.5;
	if(ierflg==0){
		hlambda->Fill(getlambdatheta,getlambdaphi,getlambdathetaphi);
		hlambda_err->Fill(getlambdathetaerr,getlambdaphierr,getlambdathetaphierr);
	}
	return ierflg;
}

void sPlot3D(int trig=2,int pt=3,int frame=0,int sys=0) {

	gStyle->SetOptStat("eMR");
	gStyle->SetStatFontSize(0.08);

	TDatime* time = new TDatime();

//	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/splot/sys_0/functional_0_%d_%d_%d.root",trig,pt,frame));
//	TGraphErrors* lambda = (TGraphErrors*)infile->Get("lambda");
//	Double_t LBTHESIG,LBPHISIG;
//	if(lambda==0x0) return;
//	else lambda->GetPoint(0,LBTHESIG,LBPHISIG); 

//	fsig->SetParameters(1,LBTHESIG+deltatheta*0.2,LBPHISIG+deltaphi*0.2,-0.5);
//	fbkg->SetParameters(1,LBTHEBKG,LBPHIBKG);
	hlambda=new TH3F("hlambda","lambdas;#lambda_{#theta};#lambda_{#phi};#lambda_{#theta#phi}",100,-1,1,100,-1,1,100,-1,1);
	hlambda2=new TH3F("hlambda2","parameters;#lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1,100,-1,1);
	hlambda_err=new TH3F("hlambda_err","lambdas uncertainty;#sigma#lambda_{#theta};#sigma#lambda_{#phi};#sigma#lambda_{#theta#phi}",100,0,2,100,0,2,100,0,2);
	hlambda2_err=new TH3F("hlambda2_err","",100,0,2,100,0,2,100,0,2);

	TCanvas *c1=new TCanvas;
	c1->Divide(4,2);
	TCanvas *c2=new TCanvas;
	c2->Divide(4,2);

	outputfile = new TFile(Form("~/polresults/20160707/functional/splot_3D_%d/splot_3D_%d_%d_%d_%d_%d.root",time->GetDate(),trig,pt,frame,sys,time->GetDate()),"recreate");
//	int trig=2,pt=3,frame=0;
	correcteddata(sys,trig,pt,frame,1);

	st = new TH1F("st","st",10,-0.5,9.5);
	method = new TH1F("method","method",10,-0.5,9.5);

	if(Minimize(sys,trig,pt,frame)!=0){
		cout<<"*** Ill-posed problem ***"<<endl;
		return ;
	}
	c1->cd(1);
	hlambda->Draw();
	c1->cd(2);
	hlambda->ProjectionX("lambda_px")->Draw();
	c1->cd(3);
	hlambda->ProjectionY("lambda_py")->Draw();
	c1->cd(4);
	hlambda->ProjectionZ("lambda_pz")->Draw();
	c1->cd(5);
	hlambda_err->Draw();
	c1->cd(6);
	hlambda_err->ProjectionX("lambda_err_px")->Draw();
	c1->cd(7);
	hlambda_err->ProjectionY("lambda_err_py")->Draw();
	c1->cd(8);
	hlambda_err->ProjectionZ("lambda_err_pz")->Draw();
	//	c1->SaveAs(Form("~/polresults/20160707/MLEmethod/splot/c1_%d.pdf",time->GetDate()));

	c2->cd(1);
	hlambda2->Draw();
	c2->cd(2);
	hlambda2->ProjectionX()->Draw();
	c2->cd(3);
	hlambda2->ProjectionY()->Draw();
	c2->cd(4);
	hlambda2_err->Draw();
	c2->cd(5);
	hlambda2_err->ProjectionX()->Draw();
	c2->cd(6);
	hlambda2_err->ProjectionY()->Draw();
	c2->SaveAs(Form("~/polresults/20160707/MLEmethod/splot/c2_%d.pdf",time->GetDate()));

	outputfile->cd();
	hlambda->Write("",TObject::kOverwrite);
	hlambda_err->Write("",TObject::kOverwrite);
	hlambda2->Write("",TObject::kOverwrite);
	hlambda2_err->Write("",TObject::kOverwrite);
	//hsub->Write("",TObject::kOverwrite);
	c1->Write("",TObject::kOverwrite);
	st->Write("",TObject::kOverwrite);
	method->Write("",TObject::kOverwrite);
	cout<<outputfile->GetName()<<endl;
}

void correcteddata(int file = 0,int trig,int pt,int frame,int rebin = 0){
	int selectedfile,file;
	if(file>=0 && file<=NFILE) selectedfile=file+1;
	else if(file==100) file=0,selectedfile=NFILE;
	else return;

	for(;file<selectedfile;file++){
		efficiencyfile = new TFile(Form("rootfile/OutFile_sys%d.root",file),"read"); 
		if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),file),"read");//marker
		else rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),0),"read");//marker
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
		rawdata3D[file][trig][pt][frame][1]->Sumw2();

		rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike_like");
		if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_hx");
		else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_cs");
		rawdata3D[file][trig][pt][frame][2]->Sumw2();
		rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);
		int max = (((TH3F*)efficiencyfile->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
		int min = (((TH3F*)efficiencyfile->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	
		if(frame==0){
			eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
			eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
			eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPt1");
			eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
		}
		else{	
			eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
			eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
			eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPtCS1");
			eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
		}
		eff3D[file][trig][pt][frame][0]->Sumw2();
		eff3D[file][trig][pt][frame][1]->Sumw2();

		rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);
		rawdata3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
		rawdata3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);

		eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
		eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
		eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
		eff3D[file][trig][pt][frame][0]->Sumw2();
		eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));
		eff3D[file][trig][pt][frame][1]->Sumw2();

		rawdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
		//					rawdata2D[file][trig][pt][frame]->SetName(Form("rawdata_%d_%s_pT%d_frame%d",file,trigName[trig].Data(),pt,frame));
		rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
		rawdata2D[file][trig][pt][frame]->Sumw2();

		bkgdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][1]->Project3D("xy");
		bkgdata2D[file][trig][pt][frame]->SetName(Form("bkg_%d_%d_%d_%d",file,trig,pt,frame));
		bkgdata2D[file][trig][pt][frame]->Sumw2();

		totdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][0]->Project3D("xy");
		totdata2D[file][trig][pt][frame]->SetName(Form("tot_%d_%d_%d_%d",file,trig,pt,frame));
		totdata2D[file][trig][pt][frame]->Sumw2();

		if(rebin==1){
			rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
			rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
			bkgdata2D[file][trig][pt][frame]->RebinX(NREBIN);
			bkgdata2D[file][trig][pt][frame]->RebinY(NREBIN);
			totdata2D[file][trig][pt][frame]->RebinX(NREBIN);
			totdata2D[file][trig][pt][frame]->RebinY(NREBIN);
		}

		eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
		eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
		eff2D[file][trig][pt][frame][0]->Sumw2();
		if(rebin==1){
			eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);//10*10
			eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);//10*10
		}
		eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
		eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
		eff2D[file][trig][pt][frame][1]->Sumw2();
		if(rebin==1){					
			eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);//10*10
			eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);//10*10
		}
		eff2D[file][trig][pt][frame][2] = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
		eff2D[file][trig][pt][frame][2]->SetName(Form("eff_ratio_%d_%d_%d_%d",file,trig,pt,frame));
		eff2D[file][trig][pt][frame][2]->Sumw2();

		eff2D[file][trig][pt][frame][2]->Divide(eff2D[file][trig][pt][frame][1]);

		datahist = (TH2F*)rawdata2D[file][trig][pt][frame]->Clone("datahist");
		bkghist = (TH2F*)bkgdata2D[file][trig][pt][frame]->Clone("bkghist");
		tothist = (TH2F*)totdata2D[file][trig][pt][frame]->Clone("tothist");
		effhist = (TH2F*)eff2D[file][trig][pt][frame][2]->Clone("effhist");

		outputfile->cd();
		datahist->Write();
		bkghist->Write();
		tothist->Write();
		effhist->Write();
	}
}


