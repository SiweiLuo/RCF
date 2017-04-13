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
//TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*cos(x)*y*(1-y**2)**0.5)",XMIN,XMAX,YMIN,YMAX);
TF2* fsig=new TF2("fsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
TF2* fsigcompare=new TF2("fsigcompare","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);
//TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* fbkg=new TF2("fbkg","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",XMIN,XMAX,YMIN,YMAX);
//TF2* lsig=new TF2("lsig","[0]*(1+[1]*cos(x)**2+[2]*(1-cos(x)**2)*cos(2*y))",XMIN,XMAX,YMIN,YMAX);
TF2* lsig=new TF2("lsig","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x)+[3]*sin(2*acos(y))*cos(x))",XMIN,XMAX,YMIN,YMAX);

TH2F *hunl, *hlik, *hsig, *hbkg, *hsub;
TH2F *hlambda, *hlambda2, *hlambda_err, *hlambda2_err;
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

float likelihoodfcn(Double_t *x, Double_t *par){
	lsig->SetParameter(0,1);
	lsig->SetParameter(1,x[0]);
	lsig->SetParameter(2,x[1]);

	double phi,costheta;
	double scale=0.;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			scale += effhist->GetBinContent(i,j)*lsig->Eval(phi,costheta);
		}
	}

	//	lsig->SetParameter(1,x[0]);
	//	lsig->SetParameter(2,x[1]);
	lsig->SetParameter(0,1./scale);	

	double eff;
	float result=0;
	double sum=0;
	for(int i=1;i<=NBIN;i++){
		for(int j=1;j<=NBIN;j++){
			//			result += -1*hsub->GetBinContent(i,j)*TMath::Log(1./scale*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)));
			if(effhist->GetBinContent(i,j)>1e-6) eff = effhist->GetBinContent(i,j);
			result += -1*hsub->GetBinContent(i,j)*TMath::Log(eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j)));
			sum += eff*lsig->Eval(hsub->GetXaxis()->GetBinCenter(i),hsub->GetYaxis()->GetBinCenter(j));
		}
	}
	//	result*=hsub->GetEntries()/(hunl->GetEntries()+hlik->GetEntries());
	result*=hsub->Integral()/(hunl->Integral()+hlik->Integral());
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

void pseudo(int n=0,int trig,int pt,int frame){
	
hunl=new TH2F("hunl","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		hlik=new TH2F("hlik","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		hsig=new TH2F("hsig","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		hbkg=new TH2F("hbkg","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		hsub=new TH2F("hsub","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);
		hsubcompare=new TH2F("hsubcompare","",NBIN,XMIN,XMAX,NBIN,YMIN,YMAX);

//	hunl->Clear();
//		hlik->Clear();
//		hsig->Clear();
//		hbkg->Clear();
//		hsub->Clear();
//		hsubcompare->Clear();

	NSIGNAL = datahist->Integral();	
	NBKG = bkghist->Integral();

	for(int i=0;i<gRandom->Poisson(NSIGNAL);) {
		Double_t xsig,ysig,xsig1,ysig1;
		fsig->GetRandom2(xsig,ysig);
		fsigcompare->GetRandom2(xsig1,ysig1);
		//		effhist->FindBin(xsig,ysig);				
		if(gRandom->Uniform()<effhist->GetBinContent(effhist->FindBin(xsig,ysig))){
			hsig->Fill(xsig,ysig);
			hunl->Fill(xsig,ysig);
			hsub->Fill(xsig,ysig);
			hsubcompare->Fill(xsig1,ysig1);
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
			hsubcompare->Fill(xbkg,ybkg);
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
	hsubcompare->Add(hlik,-1);

	//	outputfile->cd();
	//	hsub->Write();
	//	hlik->Write();

	double lbthesig,lbphisig,lbthesig_err,lbphisig_err;
	double deltacontour = 0.5;
	TF2* flikelihood=new TF2("flikelihood",likelihoodfcn,-1,1,-1,1,2);
//	flikelihood->Write();
	flikelihood->GetMinimumXY(lbthesig,lbphisig);
	hlambda->Fill(lbthesig,lbphisig);
	const double contour[1] = {flikelihood->Eval(lbthesig,lbphisig) + deltacontour};
	flikelihood->SetContour(1,contour);
	c->cd();
	flikelihood->Draw("CONT Z LIST");
	c->Update();

	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	Int_t TotalConts = 0;
	if (conts == NULL){
		printf("*** No Contours Were Extracted!\n");
		TotalConts = 0;
		return 0;
	} else {
		TotalConts = conts->GetSize();
	}
	printf("TotalConts = %d\n", TotalConts);
	TList* contlevel = (TList*)conts->At(0);
	if(contlevel->First()==NULL) {
		cout<<"gc = NULL "<<endl;
	} else {
		gc = (TGraph*)contlevel->First()->Clone();
		double xcont,ycont,xmin,xmax,ymin,ymax;
		std::vector<double> xvector,yvector;
		for(int i=0;i<gc->GetN();i++){
			gc->GetPoint(i,xcont,ycont);
			xvector.push_back(xcont);
			yvector.push_back(ycont);
		}
		xmin = getminimum(xvector);
		xmax = getmaximum(xvector);
		ymin = getminimum(yvector);
		ymax = getmaximum(yvector);
		lbthesig_err = 0.5*(xmax-xmin);
		lbphisig_err = 0.5*(ymax-ymin);
		hlambda_err->Fill(lbthesig_err,lbphisig_err);
		cout<<n<<" "<<lbthesig_err<<" "<<lbphisig_err<<endl;

		double lbthesig2,lbphisig2,lbthesig2_err,lbphisig2_err;
		double deltacontour2 = 0.5;
		TF2* flikelihood2=new TF2("flikelihood2",likelihoodfcn,TMath::Max(-1.,lbthesig-3*lbphisig_err),TMath::Min(1.,lbthesig+3*lbphisig_err),TMath::Max(-1.,lbphisig-3*lbphisig_err),TMath::Min(1.,lbphisig+3*lbphisig_err),2);
		flikelihood2->GetMinimumXY(lbthesig2,lbphisig2);
		hlambda2->Fill(lbthesig2,lbphisig2);

		const double contour2[1] = {flikelihood2->Eval(lbthesig2,lbphisig2) + deltacontour2};
		flikelihood2->SetContour(1,contour2);
		c->cd();
		flikelihood2->Draw("CONT Z LIST");
		c->Update();

		TObjArray *conts2 = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
		Int_t TotalConts2 = 0;
		if (conts2 == NULL){
			printf("*** No Contours Were Extracted!\n");
			TotalConts2 = 0;
			return 0;
		} else {
			TotalConts2 = conts2->GetSize();
		}
		printf("TotalConts2 = %d\n", TotalConts2);
		TList* contlevel2 = (TList*)conts2->At(0);
		if(contlevel2->First()==NULL) {
			cout<<"gc2 = NULL "<<endl;
		} else {
			gc2 = (TGraph*)contlevel2->First()->Clone();
			double xcont2,ycont2,xmin2,xmax2,ymin2,ymax2;
			std::vector<double> xvector2,yvector2;
			for(int i=0;i<gc2->GetN();i++){
				gc2->GetPoint(i,xcont2,ycont2);
				xvector2.push_back(xcont2);
				yvector2.push_back(ycont2);
			}
			xmin2 = getminimum(xvector2);
			xmax2 = getmaximum(xvector2);
			ymin2 = getminimum(yvector2);
			ymax2 = getmaximum(yvector2);
			lbthesig2_err = 0.5*(xmax2-xmin2);
			lbphisig2_err = 0.5*(ymax2-ymin2);
			hlambda2_err->Fill(lbthesig2_err,lbphisig2_err);
		}
//		if(n%100==0) {
			TCanvas *c0=new TCanvas;
			c0->Divide(2,2);
			c0->cd(1);
			hsub->Draw("COLZ");
			c0->cd(2);
			hlik->Draw("COLZ");
			c0->cd(3);
			flikelihood->Draw("COLZ");

			gc->SetName("cont_min");
			gc->SetTitle(Form("contour of minL+%.1f",deltacontour));
			gc->SetLineWidth(1);
			gc->Draw("SAME");
			c0->cd(4);
			flikelihood2->Draw("COLZ");
			gc2->SetTitle(Form("contour of minL+%.1f",deltacontour2));
			gc2->SetLineWidth(1);
			gc2->Draw("SAME");
			c0->SaveAs("~/polresults/20160707/MLEmethod/splot/c0.pdf");

			outputfile->cd();
			hsub->Write("",TObject::kOverwrite);
			hsubcompare->Write("",TObject::kOverwrite);
//			hlik->Write();
			flikelihood->Write("",TObject::kOverwrite);
			flikelihood2->Write("",TObject::kOverwrite);
//		}
	}
}

//void test(Double_t LBTHESIG=-0.5,Double_t LBPHISIG=0.4) {
void test(int trig=2,int pt=3,int frame=0,int deltatheta=0,int deltaphi=0) {
	
	TDatime* time = new TDatime();

	TFile* infile = new TFile(Form("/star/u/siwei/polresults/20160707/splot/sys_0/functional_0_%d_%d_%d.root",trig,pt,frame));
	TGraphErrors* lambda = (TGraphErrors*)infile->Get("lambda");
	Double_t LBTHESIG,LBPHISIG;
	if(lambda==0x0) return;
	else lambda->GetPoint(0,LBTHESIG,LBPHISIG); 

	fsig->SetParameters(1,LBTHESIG+deltatheta*0.2,LBPHISIG+deltaphi*0.2,0);
	fsigcompare->SetParameters(1,LBTHESIG+deltatheta*0.2,LBPHISIG+deltaphi*0.2,-0.5);
	fbkg->SetParameters(1,LBTHEBKG,LBPHIBKG);
	hlambda=new TH2F("hlambda","parameters;#lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1);
	hlambda2=new TH2F("hlambda2","parameters;#lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1);
	hlambda_err=new TH2F("hlambda_err","",100,0,0.5,100,0,0.5);
	hlambda2_err=new TH2F("hlambda2_err","",100,0,0.5,100,0,0.5);

	TCanvas *c1=new TCanvas;
	c1->Divide(3,2);
	TCanvas *c2=new TCanvas;
	c2->Divide(3,2);

	outputfile = new TFile(Form("~/polresults/20160707/MLEmethod/splot_%d/splot_%d_%d_%d_%d_%d_%d.root",time->GetDate(),trig,pt,frame,deltatheta,deltaphi,time->GetDate()),"recreate");
	int trig=2,pt=3,frame=0;
	correcteddata(0,trig,pt,frame,1);

	for(int n=0;n<NEXP;n++) {
				pseudo(n,trig,pt,frame);
		if(n%100==0) {
			c1->cd(1);
			hlambda->Draw();
			c1->cd(2);
			hlambda->ProjectionX()->Draw();
			c1->cd(3);
			hlambda->ProjectionY()->Draw();
			c1->cd(4);
			hlambda_err->Draw();
			c1->cd(5);
			hlambda_err->ProjectionX()->Draw();
			c1->cd(6);
			hlambda_err->ProjectionY()->Draw();
			c1->SaveAs(Form("~/polresults/20160707/MLEmethod/splot/c1_%d.pdf",time->GetDate()));

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
		}
	}
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


