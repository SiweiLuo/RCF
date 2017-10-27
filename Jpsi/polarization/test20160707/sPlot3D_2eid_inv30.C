#include <algorithm>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TGraph.h"
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

TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*TMath::Cos(2*x)+[3]*TMath::Sin(2*TMath::ACos(y))*TMath::Cos(x))",XMIN,XMAX,YMIN,YMAX);
TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*TMath::Cos(2*x)+2*[3]*TMath::Cos(x)*y*(1-y**2)**0.5)",XMIN,XMAX,YMIN,YMAX);

TH2F *hunl, *hlik, *hsig, *hbkg, *hsub;
TH3F *hlambda,  *hlambda_err,  *hlambda_bf, *hlambda_err_bf;
TH1F *hlambda_inv, *hlambda_inv_err, *hlambda_inv_bf, *hlambda_inv_err_bf;
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
TString pTName[] = {"0 < p_{T} < 2","2 < p_{T} < 3","3 < p_{T} < 4","4 < p_{T} < 6","6 < p_{T} < 8","8 < p_{T} < 14"};
TString frameName[] = {"HX frame","CS frame"};
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
TH2F* sighist;

TFile* outputfile;

TH1F* st;
TH1F* method;

TDatime* time;

Double_t eID = 1;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	f = likelihoodfcn(par);
}

float likelihoodfcn(Double_t *x){
	lsig->SetParameter(0,1);
	lsig->SetParameter(1,x[0]);
	lsig->SetParameter(2,x[1]);
	lsig->SetParameter(3,x[2]);

	cout<<" lsig parameters := "<<lsig->GetParameter(0)<<"   "<<lsig->GetParameter(1)<<"   "<<lsig->GetParameter(2)<<"   "<<lsig->GetParameter(3)<<endl;

	double phi,costheta;
	double scale=0.;
	double pdfv;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			if(effhist->GetBinContent(i,j)!=effhist->GetBinContent(i,j))continue;
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
//			cout<<" evaluation ===== "<<lsig->Eval(phi,costheta)<<endl;
			if(lsig->Eval(phi,costheta)!=lsig->Eval(phi,costheta)) {
				cout<<" phi ======== "<<phi<<"  costheta ======== "<<costheta<<endl;
				cout<<"  ACos() == "<<TMath::ACos(costheta)<<endl;
				cout<<" Sin ===== "<<TMath::Sin(TMath::ACos(costheta))<<endl;
				cout<<" !!!!!!!!!!!!!!!!!!!!!!! "<<endl;
			}
				//			pdfv = 0.5*(lsig->Eval(phi,costheta)+fabs(lsig->Eval(phi,costheta)));
//			if(lsig->Eval(phi,costheta)<0) pdfv=0.;
//			else pdfv = lsig->Eval(phi,costheta);
			scale += effhist->GetBinContent(i,j)*lsig->Eval(phi,costheta);
	//		scale += effhist->GetBinContent(i,j)*pdfv;
	//		if(lsig->Eval(phi,costheta)!=lsig->Eval(phi,costheta))cout<<"phi ==="<<phi<<"   costheta=== "<<costheta<<"  lsig->Eval() === "<<lsig->Eval(phi,costheta)<<endl;
		}
	}
	lsig->SetParameter(0,1./scale); //renormalization	

	double eff;
	float result=0;
	double sum=0;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			sum += eff*lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)); // check the normalization
			if(effhist->GetBinContent(i,j)>1e-6) eff = effhist->GetBinContent(i,j);
			else continue;
			if(datahist->GetBinContent(i,j)<1e-6) continue;
//			if(lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j))<0) cout<<" log negative == "<<lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j))<<endl;
	/*
			if(lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)<0)) pdfv=0;
			else pdfv = lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j));
	*/
			//			if(sighist->GetBinContent(i,j)<0){
			//				cout<<"sighist negative bins*********"<<endl;
			//			}
//							result += -1*sighist->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(sighist->GetXaxis()->GetBinCenter(i),sighist->GetYaxis()->GetBinCenter(j)));
			result += -1*datahist->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)));
//			result += -1*datahist->GetBinContent(i,j)*TMath::Log(eff*pdfv);
			//			sum += eff*lsig->Eval(sighist->GetXaxis()->GetBinCenter(i),sighist->GetYaxis()->GetBinCenter(j)); // check the normalization
//			if(lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j))!=lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j))) cout<<" !!!!!!!!!!!!!!!!!!!!! "<<endl;
		}
	}
	//	result*=datahist->Integral()/(tothist->Integral()+bkghist->Integral());
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

void compare(Double_t *par, Int_t *ptbin){
	gStyle->Reset();
	gStyle->SetOptStat(false);

	TCanvas* comparison = new TCanvas("comparison","comparison",1200,600);
	comparison->Divide(2,1);

	gRandom = new TRandom3();
	gRandom->SetSeed();

	hsig=new TH2F("hsig","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
	fsig->SetParameters(1,par[0],par[1],par[2]);

	Double_t eff;
	Double_t check1;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			eff = effhist->GetBinContent(i,j);
			hsig->Fill(effhist->GetXaxis()->GetBinCenter(i),effhist->GetYaxis()->GetBinCenter(j),eff*fsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)));
			//			sum += eff*lsig->Eval(datahist->GetXaxis()->GetBinCenter(i),datahist->GetYaxis()->GetBinCenter(j)); // check the normalization
		}
	}
	hsig->Scale(datahist->Integral()/hsig->Integral());


	/*
	   for(int i=0;i<datahist->Integral();) {
	   Double_t xsig,ysig;
	   fsig->GetRandom2(xsig,ysig);
	   if(gRandom->Uniform()<effhist->GetBinContent(effhist->FindBin(xsig,ysig))){
	   hsig->Fill(xsig,ysig);
	   i++;
	   }
	   }*/
	comparison->cd(1);
	datahist->ProjectionX("px");
	px->SetTitle(Form("%s @ %s GeV/c in %s",trigName[ptbin[0]].Data(),pTName[ptbin[1]].Data(),frameName[ptbin[2]].Data()));
	px->SetMinimum(0);	
	px->Draw();
	hsig->ProjectionX("hx")->Draw("same");
	comparison->cd(2);
	datahist->ProjectionY("py");
	py->SetTitle("");
	py->SetMinimum(0);
	py->Draw();
	hsig->ProjectionY("hy")->Draw("same");
	outputfile->cd();
	hsig->Write();
	fsig->Write();
	comparison->Write();
	comparison->SaveAs(Form("~/polresults/20160707/functional/splot_3D_%d_2eid_inv30/pdf/comparison_%d_%d_%d.pdf",time->GetDate(),ptbin[0],ptbin[1],ptbin[2]));
}


Int_t Minimize(int sys=0,int trig,int pt,int frame){
	TMinuit *gMinuit = new TMinuit(3);
	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,10,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[4] = {0.1, 0.1, 0.1, 0.01};
	static Double_t step[4] = {0.2, 0.2, 0.2, 0.001};
	//	static Double_t step[4] = {1, 1, 1, 0.001};
	gMinuit->mnparm(0, "p1", vstart[0], step[0], -5.,5.,ierflg);
	gMinuit->mnparm(1, "p2", vstart[1], step[1], -5.,5.,ierflg);
	gMinuit->mnparm(2, "p3", vstart[2], step[2], -5.,5.,ierflg);

	Double_t scan1list[4]={0,1000,-5.,5.};
	Double_t scan2list[4]={1,1000,-5.,5.};
	Double_t scan3list[4]={2,1000,-5.,5.};
	//	if(!(trig==1&&pt==1&&frame==1)){
	gMinuit->mnexcm("SCAN",scan1list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan2list,4,ierflg);	
	gMinuit->mnexcm("SCAN",scan3list,4,ierflg);	
	//	}
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
//	gMinuit->mnexcm("MINOS",minoslist,4,ierflg);

	Double_t getlambdatheta,getlambdathetaerr,getlambdaphi,getlambdaphierr,getlambdathetaphi,getlambdathetaphierr,getlambdainv,getlambdainverr;

	gMinuit->GetParameter(0,getlambdatheta,getlambdathetaerr);
	gMinuit->GetParameter(1,getlambdaphi,getlambdaphierr);
	gMinuit->GetParameter(2,getlambdathetaphi,getlambdathetaphierr);

	Double_t cov[3][3];
	gMinuit->mnemat(&cov[0][0],3);

	cout<<" before correction : "<<getlambdathetaerr<<"   "<<getlambdaphierr<<"   "<<getlambdathetaphierr<<endl;

	hlambda_bf = new TH3F("hlambda_bf","lambdas;#lambda_{#theta};#lambda_{#phi};#lambda_{#theta#phi}",400,-5,5,400,-5,5,400,-5,5);
	hlambda_err_bf = new TH3F("hlambda_err_bf","lambdas uncertainty;#sigma#lambda_{#theta};#sigma#lambda_{#phi};#sigma#lambda_{#theta#phi}",100,0,2,100,0,2,100,0,2);
	hlambda_inv_bf = new TH1F("hlambda_inv_bf","hlambda_inv_bf",400,-2,2);
	hlambda_inv_err_bf = new TH1F("hlambda_inv_err_bf","hlambda_inv_err_bf",400,0,5);

	hlambda_bf->Fill(getlambdatheta,getlambdaphi,getlambdathetaphi);
	hlambda_err_bf->Fill(getlambdathetaerr,getlambdaphierr,getlambdathetaphierr);

	hlambda_inv_bf->Fill((getlambdatheta+3*getlambdaphi)/(1-getlambdaphi));
	hlambda_inv_err_bf->Fill((cov[0][0]*(1/(1-getlambdaphi))**2+cov[1][1]*((3+getlambdatheta)/(1-getlambdaphi)**2)**2+2*(1/(1-getlambdaphi))*((3+getlambdatheta)/((1-getlambdaphi)**2))*cov[0][1])**0.5);

	if(trig!=1 && pt!=1){ //correction of bias
		TFile* mcvsdata = new TFile(Form("/star/u/siwei/polresults/20160707/MLEmethod/splot_3D_20171020_2eid_inv30_wofactor/thetavstheta_%d_%d_%d_20171023.root",trig,pt,frame),"read");	
		TF1 *thetaline = (TF1*)mcvsdata->Get("fit1");
		TF1 *philine = (TF1*)mcvsdata->Get("fit2");
		TF1 *thetaphiline = (TF1*)mcvsdata->Get("fit3");
		cout<<" get parameter initial : "<<getlambdatheta<<"      "<<getlambdaphi<<"    "<<getlambdathetaphi<<endl;
		TF1 *thetaget = new TF1("thetaget","[0]*x+[1]",-3,3);
		thetaget->SetParameters(1/thetaline->GetParameter(0),-1*thetaline->GetParameter(1)/thetaline->GetParameter(0));
		getlambdatheta = thetaget->Eval(getlambdatheta);
		TF1 *phiget = new TF1("phiget","[0]*x+[1]",-3,3);
		phiget->SetParameters(1/philine->GetParameter(0),-1*philine->GetParameter(1)/philine->GetParameter(0));
		getlambdaphi = thetaget->Eval(getlambdaphi);
		TF1 *thetaphiget = new TF1("thetaphiget","[0]*x+[1]",-3,3);
		thetaphiget->SetParameters(1/thetaphiline->GetParameter(0),-1*thetaphiline->GetParameter(1)/thetaphiline->GetParameter(0));
		getlambdathetaphi = thetaget->Eval(getlambdathetaphi);

		TF1 *thetaerr = (TF1*)mcvsdata->Get("fit4");
		TF1 *phierr = (TF1*)mcvsdata->Get("fit5");
		TF1 *thetaphierr = (TF1*)mcvsdata->Get("fit6");

		cout<<" get parameter 0 : "<<thetaerr->GetParameter(0)<<"    "<<phierr->GetParameter(0)<<"    "<<thetaphierr->GetParameter(0)<<endl;
		cout<<" get parameter 1 : "<<thetaerr->GetParameter(1)<<"    "<<phierr->GetParameter(1)<<"    "<<thetaphierr->GetParameter(1)<<endl;

		thetaerr->SetParameter(0,1/thetaerr->GetParameter(0));
		phierr->SetParameter(0,1/phierr->GetParameter(0));
		thetaphierr->SetParameter(0,1/thetaphierr->GetParameter(0));

		getlambdathetaerr = thetaerr->Eval(getlambdathetaerr);	
		getlambdaphierr = phierr->Eval(getlambdaphierr);
		getlambdathetaphierr = thetaphierr->Eval(getlambdathetaphierr);
	}

	cout<<" after correction : "<<getlambdathetaerr<<"   "<<getlambdaphierr<<"   "<<getlambdathetaphierr<<endl;
	st->Fill(ierflg);

	getlambdainv = (getlambdatheta+3*getlambdaphi)/(1-getlambdaphi);

	cout<<" covariant matrix entry ========"<<cov[0][1]<<endl;

	//	getlambdainverr = ((getlambdathetaerr/(1-getlambdaphi))**2+((3+getlambdatheta)*getlambdaphierr/(1-getlambdaphi)**2)**2+2*(1/(1-getlambdaphi))*((3+getlambdatheta)/(1-getlambdaphi)**2)*cov[0][1])**0.5; 
	getlambdainverr = (cov[0][0]*(1/(1-getlambdaphi))**2+cov[1][1]*((3+getlambdatheta)/(1-getlambdaphi)**2)**2+2*(1/(1-getlambdaphi))*((3+getlambdatheta)/((1-getlambdaphi)**2))*cov[0][1])**0.5; 

	cout<<"lambdatheta    = "<<getlambdatheta<<"  +/-   "<<getlambdathetaerr<<endl;
	cout<<"lambdaphi      = "<<getlambdaphi<<"  +/-  "<<getlambdaphierr<<endl;
	cout<<"lambdathetaphi = "<<getlambdathetaphi<<"  +/-  "<<getlambdathetaphierr<<endl;
	cout<<"lambdainv = "<<getlambdainv<<" +/- "<<getlambdainverr<<endl;

	if(getlambdatheta!=getlambdatheta) return ierflg;
	double lbthesig,lbphisig,lbthetaphi,lbthesig_err,lbphisig_err,lbthetaphi_err;
	double deltacontour = 0.5;

	if(ierflg==0 || getlambdatheta==getlambdatheta){
		hlambda=new TH3F("hlambda","lambdas;#lambda_{#theta};#lambda_{#phi};#lambda_{#theta#phi}",400,-5,5,400,-5,5,400,-5,5);
		hlambda_err=new TH3F("hlambda_err","lambdas uncertainty;#sigma#lambda_{#theta};#sigma#lambda_{#phi};#sigma#lambda_{#theta#phi}",100,0,2,100,0,2,100,0,2);
		hlambda_inv= new TH1F("hlambda_inv","hlambda_inv",400,-2,2);
		hlambda_inv_err = new TH1F("hlambda_inv_err","hlambda_inv_err",400,0,5);
		hlambda->Fill(getlambdatheta,getlambdaphi,getlambdathetaphi);
		hlambda_err->Fill(getlambdathetaerr,getlambdaphierr,getlambdathetaphierr);
		hlambda_inv->Fill(getlambdainv);
		hlambda_inv_err->Fill(getlambdainverr);
		cout<<"  histgroam ===>"<<hlambda->GetMean(1)<<"    "<<hlambda->GetMean(2)<<"     "<<hlambda->GetMean(3)<<endl;
	}
	else return 1;
	Double_t parameter[3] = {getlambdatheta,getlambdaphi,getlambdathetaphi}; 
	Int_t ptbin[4] = {trig,pt,frame,sys};
	compare(parameter,ptbin);
	return ierflg;
}

void sPlot3D_2eid_inv30(int trig=3,int pt=4,int frame=0,int sys=0) {

	gStyle->SetOptStat("eMR");
	gStyle->SetStatFontSize(0.08);

	time = new TDatime();

	//	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/splot/sys_0/functional_0_%d_%d_%d.root",trig,pt,frame));
	//	TGraphErrors* lambda = (TGraphErrors*)infile->Get("lambda");
	//	Double_t LBTHESIG,LBPHISIG;
	//	if(lambda==0x0) return;
	//	else lambda->GetPoint(0,LBTHESIG,LBPHISIG); 

	//	fsig->SetParameters(1,LBTHESIG+deltatheta*0.2,LBPHISIG+deltaphi*0.2,-0.5);
	//	fbkg->SetParameters(1,LBTHEBKG,LBPHIBKG);
	TCanvas *c1=new TCanvas;
	c1->Divide(4,2);
	TCanvas *c2=new TCanvas;
	c2->Divide(4,2);

	outputfile = new TFile(Form("~/polresults/20160707/functional/splot_3D_%d_2eid_inv30/splot_3D_%d_%d_%d_%d_%d.root",time->GetDate(),trig,pt,frame,sys,time->GetDate()),"recreate");
	correcteddata(sys,trig,pt,frame,1);

	st = new TH1F("st","st",10,-0.5,9.5);
	method = new TH1F("method","method",10,-0.5,9.5);

	if(Minimize(sys,trig,pt,frame)!=0){
		cout<<"*** Ill-posed problem ***"<<endl;
		return ;
	}

	if(hlambda==0x0) return;

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

	outputfile->cd();
	hlambda->Write("",TObject::kOverwrite);
	hlambda_err->Write("",TObject::kOverwrite);
	hlambda_inv->Write("",TObject::kOverwrite);
	hlambda_inv_err->Write("",TObject::kOverwrite);

	hlambda_bf->Write("",TObject::kOverwrite);
	hlambda_err_bf->Write("",TObject::kOverwrite);
	hlambda_inv_bf->Write("",TObject::kOverwrite);
	hlambda_inv_err_bf->Write("",TObject::kOverwrite);

	c1->Write("",TObject::kOverwrite);
	st->Write("",TObject::kOverwrite);
	method->Write("",TObject::kOverwrite);
}

void correcteddata(int file = 0,int trig,int pt,int frame,int rebin = 0){
	int eid=2;
	int selectedfile,file;
	//	trig += 1 ;
	if(file>=0 && file<=NFILE) selectedfile=file+1;
	else if(file==100) file=0,selectedfile=NFILE;
	else return;

	for(;file<selectedfile;file++){
		cout<<"  efficiency file ="<<Form("/star/u/siwei/efficiency/test20160706/rootfile20170817/OutFile_sys%d.root",file)<<endl;
		//	efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170510/OutFile_sys%d.root",file),"read"); 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170531/OutFile_sys%d.root",file),"read"); 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170629/OutFile_sys%d.root",file),"read"); 
		//	efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170817/OutFile_sys%d.root",file),"read"); 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170818/OutFile_sys%d.root",file),"read");// 1eID 2.9 < mee < 3.2 GeV/c^2  
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170825/OutFile_sys%d.root",file),"read");// 3.0 < mee < 3.15 GeV/c^2
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170825/OutFile_sys%d.root",file),"read");// 2.9 < mee < 3.20 GeV/c^2 2eid; 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170828/OutFile_sys%d_inv29.root",file),"read");// 2.9 < mee < 3.20 GeV/c^2 2eid; 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170913/OutFile_sys%d_inv30.root",file),"read");// 3.0 < mee < 3.15 GeV/c^2 2eid; 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170918/OutFile_sys%d_inv30.root",file),"read");// 3.0 < mee < 3.15 GeV/c^2 2eid; 
		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20171017/OutFile_sys%d_inv30.root",file),"read");// 3.0 < mee < 3.15 GeV/c^2 2eid; 
		//		efficiencyfile = new TFile(Form("/star/u/siwei/efficiency/test20160706/rootfile20170828/OutFile_sys%d_inv30.root",file),"read");// 3.0 < mee < 3.15 GeV/c^2 2eid; 

		//		if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),file),"read");//marker
		//		if((file>=20 && file<=23) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile20170521/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker

		//		if((file>=10 && file<=17) || file==22) rawdatafile[file][trig] = new TFile(Form("rootfile20170817/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		//		else rawdatafile[file][trig] = new TFile(Form("rootfile20170817/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker

		// 1eID data file
		/* 
		   if((file>=10 && file<=17) || file==22) rawdatafile[file][trig] = new TFile(Form("../invmass3_0_1eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		   else rawdatafile[file][trig] = new TFile(Form("../invmass3_0_1eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
		   */
		// 1eID data file 
		// 2eID data file 3.0 < mee < 3.15 GeV/c^2 

		if((file>=3 && file<=10) || file==15) rawdatafile[file][trig] = new TFile(Form("../invmass_2eid_inv30/rootfile20170922/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		//		else rawdatafile[file][trig] = new TFile(Form("../invmass3_0_2eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
		else rawdatafile[file][trig] = new TFile(Form("../invmass_2eid_inv30/rootfile20170922/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker

		// 2eID data file 


		// 2eID data file 
		/*
		   if((file>=10 && file<=17) || file==22) rawdatafile[file][trig] = new TFile(Form("rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		   else rawdatafile[file][trig] = new TFile(Form("rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
		   */
		// 2eID data file 

		// 2eID data file  2.9 < mee < 3.2 GeV/c^2
		/*
		   if((file>=10 && file<=17) || file==22) rawdatafile[file][trig] = new TFile(Form("../invmass2_9_2eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		   else rawdatafile[file][trig] = new TFile(Form("../invmass2_9_2eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
		   */
		// 2eID data file  2.9 < mee < 3.2 GeV/c^2

		// 1eID data file 2.9 < mee < 3.2 GeV/c^2
		/*
		   if((file>=10 && file<=17) || file==22) rawdatafile[file][trig] = new TFile(Form("../invmass2_9_1eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),file),"read");//marker
		   else rawdatafile[file][trig] = new TFile(Form("../invmass2_9_1eid/rootfile20170824/%s_sys%d.ana.root",trigSet[trig].Data(),0),"read");//marker
		   */
		// 1eID data file 2.9 < mee < 3.2 GeV/c^2 

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
		int max = (((TH3F*)efficiencyfile->Get("hHT0JpsiCosThetaPhiPt2"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
		int min = (((TH3F*)efficiencyfile->Get("hHT0JpsiCosThetaPhiPt2"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	

		// 1eID embedding 
		if(eid==1){
			if(frame==0){
				eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
				eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][0]->Print();
				eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPt1"); // the same as hJpsiCosThetaPhiPt2
				eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
			}
			else{	
				eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
				eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPtCS1");
				eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
			}
		}
		// 1eID embedding

		// 2eID embedding
		else if(eid==2){
			if(frame==0){
				eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPt2",trigName[trig].Data()));
				eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][0]->Print();
				eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPt2"); // the same as hJpsiCosThetaPhiPt2
				eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt2_%s",trigName[trig].Data()));
			}
			else{	
				eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPtCS2",trigName[trig].Data()));
				eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPtCS2");
				eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS2_%s",trigName[trig].Data()));
			}
		}
		// 2eID embedding


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
		/*
		   TMinuit *gMinuit = new TMinuit(101);
		   gMinuit->SetFCN(signalfcn);

		   Double_t arglist[10];
		   Int_t ierflg = 0;
		   arglist[0] = 1;
		   gMinuit->mnexcm("SET ERR", arglist,10,ierflg);

		   Double_t vstart[101];
		   Double_t step[101];

		   for(int i=0;i<=99;i++){
		   vstart[i] = 0.5;
		   step[i]= 0.01;
		   gMinuit->mnparm(i,Form("lambda%d",i),vstart[i],step[i],0,1,ierflg);
		   }
		   gMinuit->mnparm(100,"lambda100",-5.,0.01,-10.,0.,ierflg);

		   arglist[0]=500000;
		   arglist[1]=10.;
		   gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

		   Double_t amin,edm,errdef;
		   Int_t nvpar,nparx,icstat;
		   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);


		   Double_t lambda[101];
		   Double_t lambda_err[101];
		   Double_t par,par_err;	
		   for(int i=0;i<=99;i++) {
		   gMinuit->GetParameter(i,par,par_err);
		   lambda[i] = 1 - par;
		   cout<<"lambda ======="<<lambda[i]<<endl;	

		   }
		   sighist=new TH2F("signal","signal",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		   sighist->Sumw2();
		   for(int i=0;i<=9;i++){
		   for(int j=0;j<=9;j++){
		   sighist->SetBinContent(i+1,j+1,lambda[10*i+j]*tothist->GetBinContent(i+1,j+1));
		   }
		   }
		//		sighist->Scale(datahist->Integral()/sighist->Integral());
		*/
		outputfile->cd();
		datahist->Write();
		bkghist->Write();
		tothist->Write();
		effhist->Write();
		//		sighist->Write();

	}
}


void signalfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	f = likelisignal(par);
}

float likelisignal(Double_t *x)
{
	Double_t x2=0;
	Double_t sum=0;

	Int_t k=0;
	for(int i=0;i<=9;i++){
		for(int j=0;j<=9;j++){
			//			x2+=(tothist->GetBinContent(i,j)*x[i*10+j]-bkghist->GetBinContent(i,j))**2/bkghist->GetBinError(i,j)**2;
			//			if(i==9&&j==9) continue;
			if(tothist->GetBinContent(i+1,j+1)>=1e-6) ++k;
			x2+=(tothist->GetBinContent(i+1,j+1)*x[k]-bkghist->GetBinContent(i+1,j+1))**2;
			sum+=tothist->GetBinContent(i+1,j+1)*x[k];	
			//			x2+=(tothist->GetBinContent(i+1,j+1)*x[i*10+j]-bkghist->GetBinContent(i+1,j+1))**2;
			//			sum+=tothist->GetBinContent(i+1,j+1)*x[i*10+j];	
		}
	}
	k+=1;
	//	x2+=x[k]*(sum-bkghist->Integral());
	//	x2+=x[100]*(sum-bkghist->Integral());
	//	x2+=x[100]*(sum-bkghist->Integral());

	return x2;
}



