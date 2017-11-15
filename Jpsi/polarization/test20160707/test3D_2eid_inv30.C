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

#define NEXP     100
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
TF2* fsig=new TF2("fsig","max(0,[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x)))",XMIN,XMAX,YMIN,YMAX);
//TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* fbkg=new TF2("fbkg","max(0,[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* lsig=new TF2("lsig","max([0]*0.0001,[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x)))",XMIN,XMAX,YMIN,YMAX);
//TF3* lsig=new TF3("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);

TH2F *hunl, *hlik, *hsig, *hsigbefore, *hbkg, *hsub;
TH3F *hlambda, *hlambda2, *hlambda_err, *hlambda2_err;
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

TString trigName[NTRIG] = {"MB","HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht-1","ht0","ht1","ht2"};

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

//float likelihoodfcn(Double_t *x, Double_t *par){
float likelihoodfcn(Double_t *x){
	lsig->SetParameter(0,1);
	lsig->SetParameter(1,x[0]);
	lsig->SetParameter(2,x[1]);
	lsig->SetParameter(3,x[2]);

	double phi,costheta;
	double scale=0.,datatot=0.;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			scale += effhist->GetBinContent(i,j)*lsig->Eval(phi,costheta);
			datatot += hsub->GetBinContent(i,j);
		}
	}
	lsig->SetParameter(0,1./scale); //renormalization	
//	lsig->SetParameter(0,datatot/scale); //renormalization	

	double eff;
	float result=0;
	double sum=0;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			//			result += -1*hsub->GetBinContent(i,j)*TMath::Log(1./scale*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)));
			if(effhist->GetBinContent(i,j)>1e-8) eff = effhist->GetBinContent(i,j);
			else continue;
			//			if(hsub->GetBinContent(i,j)<0) continue;
//			if(lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j))<0) continue;
//			result += (hsub->GetBinContent(i,j)-eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)))**2;
//			result += -1*hsub->GetBinContent(i,j)*eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j));
			//eff=effhist->GetBinContent(i,j);
			result += -1*hsub->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)));
			sum += eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j));
		}
	}
//	result*=hsub->Integral()/(hunl->Integral()+hlik->Integral()); // "sPlot technique" 
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

Int_t pseudo(int n=0,int trig,int pt,int frame){
//	hunl->Clear();
//	hlik->Clear();
//	hsig->Clear();
//	hbkg->Clear();
//	hsub->Clear();

	hunl=new TH2F("hunl","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	hlik=new TH2F("hlik","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	hsig=new TH2F("hsig","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	hsigbefore=new TH2F("hsigbefore","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	hbkg=new TH2F("hbkg","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	hsub=new TH2F("hsub","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);

	gRandom = new TRandom3();
	gRandom->SetSeed();

	NSIGNAL = datahist->Integral();	
	NBKG = bkghist->Integral();

	for(int i=0;i<gRandom->Poisson(NSIGNAL);) {
		Double_t xsig,ysig;
		fsig->GetRandom2(xsig,ysig);
		hsigbefore->Fill(xsig,ysig);
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

	outputfile->cd();
	hsub->Write("",TObject::kOverwrite);
	hsigbefore->Write("",TObject::kOverwrite);


	TMinuit *gMinuit = new TMinuit(3);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,10,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[4] = {0.1, 0.1, 0.1, 0.01};
	static Double_t step[4] = {0.2, 0.2, 0.2, 0.001};

	gMinuit->mnparm(0, "p1", vstart[0], step[0], -5.,5.,ierflg);
	gMinuit->mnparm(1, "p2", vstart[1], step[1], -5.,5.,ierflg);
	gMinuit->mnparm(2, "p3", vstart[2], step[2], -5.,5.,ierflg);
	//   gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);

	Double_t scan1list[4]={1,1000,-5,5};
	Double_t scan2list[4]={2,1000,-5,5};
	Double_t scan3list[4]={3,1000,-5,5};
	gMinuit->mnexcm("SCAN",scan1list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan2list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan3list,4,ierflg);	

	// Now ready for minimization step
	arglist[0] = 5000;
	arglist[1] = 10.;
	gMinuit->mnexcm("MIGRAD", arglist ,4,ierflg);
	//	gMinuit->mnexcm("Minimize", arglist ,2,ierflg);
	if(ierflg==0) method->Fill(1);

	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//gMinuit->mnprin(3,amin);
	//		Int_t st = gMinuit->Command("MIGRAD");
	//		Int_t st = gMinuit->Command("MInimze");
	//		Int_t st = gMinuit->Command("Scan");

	//	if(ierflg!=0) gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
	/*	if(ierflg!=0) {
		gMinuit->mnexcm("Minimize",arglist,2,ierflg);
		if(ierflg==0) method->Fill(2);
		}
		if(ierflg!=0) {
		gMinuit->mnexcm("Simplex",arglist,2,ierflg);
		if(ierflg==0) method->Fill(3);
		}*/
//	ierflg=4;
//	if(ierflg!=0){
		gMinuit->mnexcm("HESSE",arglist,4,ierflg);
		if(ierflg==0) method->Fill(4);
//	}
//	ierflg=4;
	Double_t minoslist[4] = {5000,0,1,2};

//	if(ierflg!=0){
		gMinuit->mnexcm("MINOS",minoslist,4,ierflg);
		if(ierflg==0) method->Fill(5);
//	}

		Double_t getlambdatheta,getlambdathetaerr,getlambdaphi,getlambdaphierr,getlambdathetaphi,getlambdathetaphierr;

	gMinuit->GetParameter(0,getlambdatheta,getlambdathetaerr);
	gMinuit->GetParameter(1,getlambdaphi,getlambdaphierr);
	gMinuit->GetParameter(2,getlambdathetaphi,getlambdathetaphierr);

	st->Fill(ierflg);

	cout<<"status =================="<<ierflg<<endl;

	cout<<"lambdatheta = "<<getlambdatheta<<"    "<<getlambdathetaerr<<endl;
	cout<<"lambdatheta = "<<getlambdaphi<<"    "<<getlambdaphierr<<endl;
	cout<<"lambdatheta = "<<getlambdathetaphi<<"    "<<getlambdathetaphierr<<endl;

	double lbthesig,lbphisig,lbthetaphi,lbthesig_err,lbphisig_err,lbthetaphi_err;
	double deltacontour = 0.5;
	if(ierflg==0){
		hlambda->Fill(getlambdatheta,getlambdaphi,getlambdathetaphi);
		hlambda_err->Fill(getlambdathetaerr,getlambdaphierr,getlambdathetaphierr);
	}
	return ierflg;
}

void test3D_2eid_inv30(int trig=3,int pt=4,int frame=0,int deltatheta=0,int deltaphi=0,int deltathetaphi=0) {

	gStyle->SetOptStat("eMR");
	gStyle->SetStatFontSize(0.08);

	TDatime* time = new TDatime();

//	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/splot/sys_0/functional_0_%d_%d_%d.root",trig,pt,frame));
//	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/functional/splot_3D_20170430/splot_3D_%d_%d_%d_0_20170430.root",trig,pt,frame),"read");
	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/functional/splot_3D_20171114_2eid_inv30/splot_3D_%d_%d_%d_0_20171114.root",trig,pt,frame),"read");
//	TGraphErrors* lambda = (TGraphErrors*)infile->Get("lambda");
	TH3F* lambda = (TH3F*)infile->Get("hlambda");
	Double_t LBTHESIG,LBPHISIG,LBTHEPHISIG;
	if(lambda==0x0) return;
//	else lambda->GetPoint(0,LBTHESIG,LBPHISIG); 
	else {
		LBTHESIG = lambda->GetMean(1);
		LBPHISIG = lambda->GetMean(2);
		LBTHEPHISIG = lambda->GetMean(3);
	}

	fsig->SetParameters(1,LBTHESIG+deltatheta*0.2,LBPHISIG+deltaphi*0.2,LBTHEPHISIG+deltathetaphi*0.2);
	fbkg->SetParameters(1,LBTHEBKG,LBPHIBKG);
	hlambda=new TH3F("hlambda","lambdas central value;#lambda_{#theta};#lambda_{#phi};#lambda_{#theta#phi}",400,-5.,5.,400,-5.,5.,400,-5.,5.);
	hlambda2=new TH3F("hlambda2","parameters;#lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1,100,-1,1);
	hlambda_err=new TH3F("hlambda_err","lambdas uncertainty;#sigma#lambda_{#theta};#sigma#lambda_{#phi};#sigma#lambda_{#theta#phi}",100,0,2,100,0,2,100,0,2);
	hlambda2_err=new TH3F("hlambda2_err","",100,0,2,100,0,2,100,0,2);

	TCanvas *c1=new TCanvas;
	c1->Divide(3,2);
	TCanvas *c2=new TCanvas;
	c2->Divide(4,2);

	outputfile = new TFile(Form("~/polresults/20160707/MLEmethod/splot_3D_%d_2eid_inv30_wofactor/splot_3D_%d_%d_%d_%d_%d_%d_%d.root",time->GetDate(),trig,pt,frame,deltatheta,deltaphi,deltathetaphi,time->GetDate()),"recreate");
//	int trig=2,pt=3,frame=0;
	correcteddata(0,trig,pt,frame,1);

	st = new TH1F("st","st",10,-0.5,9.5);
	method = new TH1F("method","method",10,-0.5,9.5);

	for(int n=0;n<NEXP;) {
		if(pseudo(n,trig,pt,frame)==0) n++;
	}
//	c1->cd(1);
//	hlambda->Draw();
	c1->cd(1);
	hlambda->ProjectionX("lambda_px")->Draw();
	c1->cd(2);
	hlambda->ProjectionY("lambda_py")->Draw();
	c1->cd(3);
	hlambda->ProjectionZ("lambda_pz")->Draw();
//	c1->cd(5);
//	hlambda_err->Draw();
	c1->cd(4);
	hlambda_err->ProjectionX("lambda_err_px")->Draw();
	c1->cd(5);
	hlambda_err->ProjectionY("lambda_err_py")->Draw();
	c1->cd(6);
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
//	c2->SaveAs(Form("~/polresults/20160707/MLEmethod/splot/c2_%d.pdf",time->GetDate()));

	outputfile->cd();
	hlambda->Write("",TObject::kOverwrite);
	hlambda_err->Write("",TObject::kOverwrite);
	hlambda2->Write("",TObject::kOverwrite);
	hlambda2_err->Write("",TObject::kOverwrite);
//	hsub->Write("",TObject::kOverwrite);
	c1->SaveAs(Form("~/WWW/MonteCarlo/MC/2eid_inv30_lambdasvslambdas/%d/projection_%d_%d_%d_%d_%d_%d.pdf",time->GetDate(),trig,pt,frame,deltatheta,deltaphi,deltathetaphi));
	c1->SaveAs(Form("~/polresults/20160707/MLEmethod/splot_3D_%d_2eid_inv30_wofactor/projection_%d_%d_%d_%d_%d_%d.pdf",time->GetDate(),trig,pt,frame,deltatheta,deltaphi,deltathetaphi));
	c1->Write("",TObject::kOverwrite);
	st->Write("",TObject::kOverwrite);
	method->Write("",TObject::kOverwrite);
	cout<<outputfile->GetName()<<endl;
}

void correcteddata(int file = 0,int trig,int pt,int frame,int rebin = 0){
	int eid=2;
	int selectedfile,file;
	if(file>=0 && file<=NFILE) selectedfile=file+1;
	else if(file==100) file=0,selectedfile=NFILE;
	else return;

	for(;file<selectedfile;file++){
//		efficiencyfile = new TFile(Form("rootfile20170714/OutFile_sys%d.root",file),"read"); 
//		efficiencyfile = new TFile(Form("rootfile20170818/OutFile_sys%d.root",file),"read"); 
//		efficiencyfile = new TFile("/star/u/siwei/efficiency/test20160706/rootfile20170825/OutFile_sys0.root","read");// 2.9 < mee < 3.2 GeV/c^2 2eid;	
	//	efficiencyfile = new TFile("/star/u/siwei/efficiency/test20160706/rootfile20170911/OutFile_sys0_inv30.root","read");// 2.9 < mee < 3.2 GeV/c^2 2eid;	
	//	efficiencyfile = new TFile("/star/u/siwei/efficiency/test20160706/rootfile20170921/OutFile_sys0_inv30.root","read");// 2.9 < mee < 3.2 GeV/c^2 2eid;	
		efficiencyfile = new TFile("/star/u/siwei/efficiency/test20160706/rootfile20171017/OutFile_sys0_inv30.root","read");// 2.9 < mee < 3.2 GeV/c^2 2eid;	
//		efficiencyfile = new TFile("/star/u/siwei/efficiency/test20160706/rootfile20170828/OutFile_sys0_inv30.root","read");// 3.0 < mee < 3.15 GeV/c^2 2eid;	
	
	//	rawdatafile[file][trig] = new TFile(Form("../invmass_%deid_inv30/rootfile20170912/%s_sys%d.ana.root",eid,trigSet[trig].Data(),0),"read");// 2.9 < mee < 3.2 GeV/c^2 2eid;
		rawdatafile[file][trig] = new TFile(Form("../invmass_%deid_inv30/rootfile20170912/%s_sys%d.ana.root",eid,trigSet[trig].Data(),0),"read");// 2.9 < mee < 3.2 GeV/c^2 2eid;
//		rawdatafile[file][trig] = new TFile(Form("../invmass_%deid_inv30/rootfile20170912/%s_sys%d.ana.root",eid,trigSet[trig].Data(),0),"read");// 3. < mee < 3.15 GeV/c^2 2eid;
//		rawdatafile[file][trig] = new TFile(Form("../invmass_1eid_inv30/rootfile20170912/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");// 2.9 < mee < 3.2 GeV/c^2 1eid;
//		rawdatafile[file][trig] = new TFile(Form("../invmass3_0_1eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");// 3.0 < mee < 3.15 GeV/c^2 1eid;
//		rawdatafile[file][trig] = new TFile(Form("../invmass_2eid_inv30/rootfile20170912/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");// 3.0 < mee < 3.15 GeV/c^2 1eid;

/*	
		if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile20170818/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		else rawdatafile[file][trig] = new TFile(Form("rootfile20170818/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
*/		
		
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
	// 1eID 
if(eid==1){	
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
}
else if(eid==2){
	// 1eID 
//2eID
	
		if(frame==0){
			eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPt2",trigName[trig].Data()));
			eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
			eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPt2");
			eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt2_%s",trigName[trig].Data()));
		}
		else{	
			eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPtCS2",trigName[trig].Data()));
			eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
			eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPtCS2");
			eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS2_%s",trigName[trig].Data()));
		}

		//2eID
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


