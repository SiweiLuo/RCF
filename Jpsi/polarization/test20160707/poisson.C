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

TF2* sigmab; // background pdf
TH2F* bkghist;
TH2F* tothist;

TF2* sfcn; //signal distribution
TF2* bfcn; //background distribution

//double ltheta[2],lphi[2];
double fs[10][10];
double fb[10][10];

TMatrix V(2,2);
TMatrix invV(2,2);

double Ns,Nb;

TRandom *sample;
TRandom *totsample;
TRandom *bkgsample;
//double ltheta[2],lphi[2];

Double_t likelihoodfcn(Double_t *x, Double_t *par)
{
	Double_t result=0.;

	int nbinxeff = effhist->GetNbinsX();
	int nbinyeff = effhist->GetNbinsY();

	double costheta=0.,phi=0.,eff=0.;
	double si = 0.;
	int Npositron = 0;

	sigma = new TF2("sigma","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",-TMath::Pi(),TMath::Pi(),-1,1);

	double scale = 0.;

	sigma->SetParameter(0,1);
	sigma->SetParameter(1,x[0]);
	sigma->SetParameter(2,x[1]);

	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){

			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);

			si = sigma->Eval(phi,costheta);

			eff = effhist->GetBinContent(i,j);
			if(si*eff<1e-16) continue;
			scale += si*eff;
		}
	}
	//			sigma->SetParameter(0,1);
	sigma->SetParameter(1,x[0]);
	sigma->SetParameter(2,x[1]);
	//			scale = sigma->Integral(-TMath::Pi(),TMath::Pi(),-1,1);
	sigma->SetParameter(0,1./scale);
	scale = sigma->Integral(-TMath::Pi(),TMath::Pi(),-1,1);

	double sum = 0.;
	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){

			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			si = sigma->Eval(phi,costheta);

			eff = effhist->GetBinContent(i,j);
			Npositron = datahist->GetBinContent(i,j);
			if(si*eff<1e-16) continue;
			sum += si*eff;
			result += -1*Npositron*TMath::Log(si*eff);
		}
	}
	return result;
}

void poisson(int file,int sys,int trig,int pt,int frame){

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);

	correcteddata(sys,trig,pt,frame,1);

	TFile* savefile = new TFile(Form("~/polresults/20160707/poisson/sys_%d/functional_%d_%d_%d_%d.root",sys,file,trig,pt,frame),"recreate");
	savefile->cd();
	datahist->Write();	
	TF2* likelihood;

	double ltheta,lphi,lthetaerr,lphierr;
	likelihood = new TF2("likelihood",likelihoodfcn,-1,1,-1,1,0);
	likelihood->SetTitle("L(#lambda_{#theta},#lambda_{#phi}) = -#sum^{1}_{cos#theta = -1}#sum^{#pi}_{#phi=-#pi}N*ln(#frac{#partial^{2}#sigma}{#partialcos#theta#partial#phi}*#varepsilon);#lambda_{#theta};#lambda_{#phi}");
	likelihood->Write("",TObject::kOverwrite);
	TCanvas* c = new TCanvas();
	c->cd();

	likelihood->Draw("colz");

	likelihood->GetMinimumXY(ltheta,lphi);
	double par0ini = sigma->GetParameter(0);
	sigma->SetParameter(1,ltheta);
	sigma->SetParameter(2,lphi);	
	double normalization = sigma->Integral(-TMath::Pi(),TMath::Pi(),-1,1);
	sigma->SetParameter(0,par0ini/normalization);

	sigma->Write("",TObject::kOverwrite);

	TGraphErrors* lambda = new TGraphErrors(1,&ltheta,&lphi,0,0);
	lambda->SetName("lambda");
	lambda->SetTitle("lambda parameter;#lambda_{#theta};#lambda_{#phi}");
	lambda->GetYaxis()->SetRangeUser(-1,1);
	lambda->GetXaxis()->SetLimits(-1,1);
	lambda->Write();

	/*
	double deltacontour = 0.5;
	double contour = likelihood->Eval(ltheta,lphi) + deltacontour;	

	likelihood->SetContour(1,&contour);
	likelihood->Draw("CONT Z LIST");
	c->Update();

	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

	TList* contlevel = NULL;
	TGraph* gc2 = NULL;
	contlevel = (TList*)conts->At(0);	
	if(contlevel==NULL) continue;
	contlevel->Write("",TObject::kOverwrite);

	gc2 = (TGraph*)contlevel->First();
	if(gc2==NULL) {
		cout<<"gc2 = NULL "<<endl;	
		return;
	}

	double xcont,ycont,xmin,xmax,ymin,ymax;
	std::vector<double> xvector,yvector;	
	while(1){
		for(int i=0;i<gc2->GetN();i++){
			gc2->GetPoint(i,xcont,ycont);
			xvector.push_back(xcont);
			yvector.push_back(ycont);
		}
		xmin = getminimum(xvector);
		xmax = getmaximum(xvector);
		ymin = getminimum(yvector);
		ymax = getmaximum(yvector);	
		if(ltheta<xmin || ltheta>xmax || lphi<ymin || lphi>ymax){//minimum point doesn't locate in the contour
			gc2 = (TGraph*)contlevel->After(gc2);
			xvector.clear();
			yvector.clear();
		}
		else break;
	}
	gc2->SetName("cont_min");
	gc2->SetTitle(Form("contour of minL+%.1f",deltacontour));
	gc2->SetLineWidth(1);
	gc2->Write("",TObject::kOverwrite);

	lthetaerr = 0.5*(xmax-xmin);
	lphierr = 0.5*(ymax-ymin);
*/
	//	}

	//	TGraphErrors* lambda = new TGraphErrors(1,&ltheta,&lphi,&lthetaerr,&lphierr);
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

void correcteddata(int file = 0, int trig,int pt,int frame,int rebin = 0){
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

		tothist = (TH2F*)totdata2D[file][trig][pt][frame]->Clone("tothist");
		datahist = (TH2F*)rawdata2D[file][trig][pt][frame]->Clone("datahist");
		bkghist = (TH2F*)bkgdata2D[file][trig][pt][frame]->Clone("bkghist");
		effhist = (TH2F*)eff2D[file][trig][pt][frame][2]->Clone("effhist");

		resample();

	}
}

void resample(){

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);

	sample = new TRandom3(0);
	totsample = new TRandom3(0);
	bkgsample = new TRandom3(0);

	Int_t totdata,bkgdata;
	for(int i=0;i<10;i++){
		for(int j=0;j<10;j++){
			totdata = sample->Poisson(tothist->GetBinContent(i+1,j+1));
			bkgdata = sample->Poisson(bkghist->GetBinContent(i+1,j+1));
			cout<<"tot = "<<totdata<<"  -  bkg" <<bkgdata<<"  ="<<totdata-bkgdata<<endl;
			datahist->SetBinContent(i+1,j+1,totdata-bkgdata);

		}
	}
}
/*
   Double_t P(Double_t *par)//poisson distribution 
   {





   }
   */


/*
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
//					rawdata3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
//					rawdata3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
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
if(rebin==1){
rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
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
effhist = (TH2F*)eff2D[file][trig][pt][frame][2]->Clone("effhist");
}
}
*/


