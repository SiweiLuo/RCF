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
#define XNBIN 10
#define YNBIN 10
#define CHIXBIN 100
#define CHIYBIN 100
#define pi 3.1415926

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigSet[NTRIG] = {"ht0","ht1","ht2"};
TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString frameName[NFRAME] = {"HX","CS"};
TString phaseName[NPHASE+1] = {"#lambda_{#theta}","#lambda_{#phi}","#lambda_{inv}"};

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
double systematic_error[NTRIG][NPT][NFRAME][NPHASE+1]; // theta, phi, inv.

Double_t matrix[4][4];
Double_t covariant[NTRIG][NPT][NFRAME];

TCanvas *canvas[NFILE][NTRIG][NPT][NFRAME];//marker
TCanvas *systematic[NTRIG];
TCanvas *check;

double fParamVal[NFILE][NTRIG][NPT][NFRAME][4];
double fParamErr[NFILE][NTRIG][NPT][NFRAME][4];
double arglist[10];
Int_t ierflg=0;
Int_t fitflag[4][NPT][2][4];
TH1F* fitdata[NTRIG][NPT][NFRAME][NPHASE];// theta , phi

TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda_parameters_sys[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant

TLegend* legend;

TFile *outputfile;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

TH2F* chi2;
TH2F* likelihood_contour;

double lambda[2][2];
TF2* mlefcn;
TF2* contfcn;

//double contourlevel[3][6]={
//	{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}
//};

double contourlevel[2][3][6]={
	{{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}},
	{{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}}
};

void fitlambda(int file = 0, int trig = 0, int pt = 0, int frame = 0){
	gStyle->SetPadRightMargin(0.2);	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	TFile* inputfile = new TFile(Form("~/polresults/20160707/rootcombined/lambda_file%d_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"update");

	chi2 = (TH2F*)inputfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));

	if(chi2==0x0 || chi2->GetMean()!=chi2->GetMean() || chi2->GetMean()==0 ){
		lambda_theta[file][trig][pt][frame] = -100;
		lambda_theta_err[file][trig][pt][frame] = 0.;
		lambda_phi[file][trig][pt][frame] = -100;
		lambda_phi_err[file][trig][pt][frame] = 0.;	
		return;
	}

	std::vector<double> lambda_theta_error,lambda_phi_error;

	double xx,yy,zz;
	mlefcn = new TF2("mlefcn","[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y+[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y+[21]*x*x*x*x*x*x+[22]*x*x*x*x*x*y+[23]*x*x*x*x*y*y+[24]*x*x*x*y*y*y+[25]*x*x*y*y*y*y+[26]*x*y*y*y*y*y+[27]*y*y*y*y*y*y+[28]*x*x*x*x*x*x*x+[29]*y*x*x*x*x*x*x+[30]*y*y*x*x*x*x*x+[31]*y*y*y*x*x*x*x+[32]*x*x*x*y*y*y*y+[33]*x*x*y*y*y*y*y+[34]*x*y*y*y*y*y*y+[35]*y*y*y*y*y*y*y+[36]*x*x*x*x*x*x*x*x+[37]*y*x*x*x*x*x*x*x+[38]*y*y*x*x*x*x*x*x+[39]*y*y*y*x*x*x*x*x+[40]*y*y*y*y*x*x*x*x+[41]*x*x*x*y*y*y*y*y+[42]*x*x*y*y*y*y*y*y+[43]*x*y*y*y*y*y*y*y+[44]*y*y*y*y*y*y*y*y",-1,1,-1,1);
//	mlefcn = new TF2("mlefcn","[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y",-1,1,-1,1);//[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y+[21]*x*x*x*x*x*x+[22]*x*x*x*x*x*y+[23]*x*x*x*x*y*y+[24]*x*x*x*y*y*y+[25]*x*x*y*y*y*y+[26]*x*y*y*y*y*y+[27]*y*y*y*y*y*y+[28]*x*x*x*x*x*x*x+[29]*y*x*x*x*x*x*x+[30]*y*y*x*x*x*x*x+[31]*y*y*y*x*x*x*x+[32]*x*x*x*y*y*y*y+[33]*x*x*y*y*y*y*y+[34]*x*y*y*y*y*y*y+[35]*y*y*y*y*y*y*y+[36]*x*x*x*x*x*x*x*x+[37]*y*x*x*x*x*x*x*x+[38]*y*y*x*x*x*x*x*x+[39]*y*y*y*x*x*x*x*x+[40]*y*y*y*y*x*x*x*x+[41]*x*x*x*y*y*y*y*y+[42]*x*x*y*y*y*y*y*y+[43]*x*y*y*y*y*y*y*y+[44]*y*y*y*y*y*y*y*y",-1,1,-1,1);
	mlefcn->SetTitle("mle fit function;#lambda_{#theta};#lambda_{#phi}");
	cout<<"Fit the whole area"<<endl;
	chi2->Fit("mlefcn","0");
	mlefcn->Write("",TObject::kOverwrite);
	mlefcn->GetMinimumXY(xx,yy);
	double minimum = mlefcn->GetMinimum();
	cout<<"minimum="<<minimum<<endl;
	lambda_theta[file][trig][pt][frame] = xx;
	lambda_phi[file][trig][pt][frame] = yy;

	canvas[file][trig][pt][frame] = new TCanvas(Form("mle_%d_%d_%d_%d",file,trig,pt,frame),Form("mle_%d_%d_%d_%d",file,trig,pt,frame),3200,600);
	canvas[file][trig][pt][frame]->Divide(5,1);
	canvas[file][trig][pt][frame]->cd(1);
	chi2->Draw("surf2");

	canvas[file][trig][pt][frame]->cd(2);
	double contour = minimum + contourlevel[frame][trig][pt];	
	cout<<"contour="<<contour<<endl;
	mlefcn->SetContour(1,&contour);
	mlefcn->Draw("CONT Z LIST");
	canvas[file][trig][pt][frame]->Update();

	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

	TList* contlevel = NULL;
//	TGraph* gc = NULL;
	TGraph* gc2 = NULL;
	contlevel = (TList*)conts->At(0);	
	contlevel->Write("",TObject::kOverwrite);
//	gc = (TGraph*)contlevel->First();
//	if(gc==NULL) return;
//	gc->SetName("cont_min_half");
//	gc->SetTitle("contour of minL+0.5");
//	gc->SetLineWidth(1);
//	gc->Write("",TObject::kOverwrite);	

	gc2 = (TGraph*)contlevel->First();
//	cout<<"contlevel->First() = "<<(TGraph*)contlevel->First()<<endl;
//	cout<<"gc = "<<gc<<endl;
//	gc2 = (TGraph*)contlevel->After((TGraph*)contlevel->First());
	if(gc2==NULL) return;
	gc2->SetName("cont_min_second");
	gc2->SetTitle("contour of minL+2");
	gc2->SetLineWidth(1);
	gc2->Write("",TObject::kOverwrite);

	canvas[file][trig][pt][frame]->Update();

	double xcont,ycont,xmin,xmax,ymin,ymax;
	std::vector<double> xvector,yvector;	
	for(int i=0;i<gc2->GetN();i++){
		gc2->GetPoint(i,xcont,ycont);
		xvector.push_back(xcont);
		yvector.push_back(ycont);
	}
	xmin = getminimum(xvector);
	xmax = getmaximum(xvector);
	ymin = getminimum(yvector);
	ymax = getmaximum(yvector);	
	cout<<"xmin="<<xmin<<"xmax="<<xmax<<"ymin="<<ymin<<"ymax="<<ymax<<endl;

	contfcn = new TF2("contfcn","[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y+[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y+[21]*x*x*x*x*x*x+[22]*x*x*x*x*x*y+[23]*x*x*x*x*y*y+[24]*x*x*x*y*y*y+[25]*x*x*y*y*y*y+[26]*x*y*y*y*y*y+[27]*y*y*y*y*y*y+[28]*x*x*x*x*x*x*x+[29]*y*x*x*x*x*x*x+[30]*y*y*x*x*x*x*x+[31]*y*y*y*x*x*x*x+[32]*x*x*x*y*y*y*y+[33]*x*x*y*y*y*y*y+[34]*x*y*y*y*y*y*y+[35]*y*y*y*y*y*y*y+[36]*x*x*x*x*x*x*x*x+[37]*y*x*x*x*x*x*x*x+[38]*y*y*x*x*x*x*x*x+[39]*y*y*y*x*x*x*x*x+[40]*y*y*y*y*x*x*x*x+[41]*x*x*x*y*y*y*y*y+[42]*x*x*y*y*y*y*y*y+[43]*x*y*y*y*y*y*y*y+[44]*y*y*y*y*y*y*y*y",xmin,xmax,ymin,ymax);
//	contfcn = new TF2("contfcn","[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y",xmin,xmax,ymin,ymax);//[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y+[21]*x*x*x*x*x*x+[22]*x*x*x*x*x*y+[23]*x*x*x*x*y*y+[24]*x*x*x*y*y*y+[25]*x*x*y*y*y*y+[26]*x*y*y*y*y*y+[27]*y*y*y*y*y*y+[28]*x*x*x*x*x*x*x+[29]*y*x*x*x*x*x*x+[30]*y*y*x*x*x*x*x+[31]*y*y*y*x*x*x*x+[32]*x*x*x*y*y*y*y+[33]*x*x*y*y*y*y*y+[34]*x*y*y*y*y*y*y+[35]*y*y*y*y*y*y*y+[36]*x*x*x*x*x*x*x*x+[37]*y*x*x*x*x*x*x*x+[38]*y*y*x*x*x*x*x*x+[39]*y*y*y*x*x*x*x*x+[40]*y*y*y*y*x*x*x*x+[41]*x*x*x*y*y*y*y*y+[42]*x*x*y*y*y*y*y*y+[43]*x*y*y*y*y*y*y*y+[44]*y*y*y*y*y*y*y*y",xmin,xmax,ymin,ymax);
	cout<<"Fit the contour area"<<endl;
//	chi2->GetXaxis()->SetRange(xmin,xmax);
//	chi2->GetYaxis()->SetRange(ymin,ymax);
//	chi2->Fit("contfcn","","",xmin,xmax,ymin,ymax);
	chi2->Fit("contfcn","R0");
//	chi2->Fit("contfcn");


	int xxbin = (xmax-xmin)/0.02+1;
	int yybin = (ymax-ymin)/0.02+1;
	likelihood_contour = new TH2F("likelihood_contour","likelihood_contour",xxbin,xmin,xmax,yybin,ymin,ymax);
	for(int i=1;i<xxbin;i++){
		for(int j=1;j<yybin;j++){
			double likelihood_x = chi2->GetXaxis()->FindBin(xmin+0.02*(i-1));
			double likelihood_y = chi2->GetYaxis()->FindBin(ymin+0.02*(j-1));
			double bincontent = chi2->GetBinContent(likelihood_x,likelihood_y);
			likelihood_contour->SetBinContent(i,j,bincontent);
		}
	}
	likelihood_contour->Write("",TObject::kOverwrite);
//	likelihood_contour->Fit("contfcn");

	contfcn->SetTitle("contour fit function;#lambda_{#theta};#lambda_{#phi}");
	contfcn->Write("",TObject::kOverwrite);


	TF2* Fouri = new TF2("Fouri","[0]+[1]*sin(0.5*pi*x)+[2]*cos(0.5*pi*x)+[3]*sin(0.5*pi*y)+[4]*cos(0.5*pi*y)+[5]*sin(0.5*2*pi*x)+[6]*sin(0.5*pi*x)*sin(0.5*pi*y)+[7]*sin(0.5*pi*x)*cos(0.5*pi*y)+[8]*cos(0.5*pi*x)*sin(0.5*pi*y)+[9]*cos(0.5*pi*x)*cos(0.5*pi*y)+[10]*cos(0.5*2*pi*x)+[11]*sin(0.5*2*pi*y)+[12]*cos(0.5*2*pi*y)+[13]*sin(0.5*3*pi*x)+[14]*cos(0.5*3*pi*x)+[15]*sin(0.5*3*pi*y)+[16]*cos(0.5*3*pi*y)+[17]*sin(0.5*2*pi*x)*sin(0.5*pi*y)+[18]*sin(0.5*2*pi*x)*cos(0.5*pi*y)+[19]*cos(0.5*2*pi*x)*sin(0.5*pi*y)+[20]*cos(0.5*2*pi*x)*cos(0.5*pi*y)+[21]*sin(0.5*pi*x)*sin(0.5*2*pi*y)+[22]*sin(0.5*pi*x)*cos(0.5*2*pi*y)+[23]*cos(0.5*pi*x)*sin(0.5*2*pi*x)+[24]*cos(0.5*pi*x)*cos(0.5*2*pi*y)+[25]*sin(0.5*3*pi*y)+[26]*cos(0.5*3*pi*y)+[27]*sin(0.5*4*pi*x)+[28]*cos(0.5*4*pi*x)+[29]*sin(0.5*4*pi*x)+[30]*cos(0.5*4*pi*y)+[31]*sin(0.5*3*pi*x)*sin(0.5*pi*y)+[32]*sin(0.5*3*pi*x)*cos(0.5*pi*y)+[33]*cos(0.5*3*pi*x)*sin(0.5*pi*y)+[34]*cos(0.5*3*pi*x)*cos(0.5*pi*y)+[35]*sin(0.5*2*pi*x)*sin(0.5*2*pi*y)+[36]*sin(0.5*2*pi*x)*cos(0.5*2*pi*y)+[37]*cos(0.5*2*pi*x)*sin(0.5*2*pi*y)+[38]*cos(0.5*2*pi*x)*cos(0.5*2*pi*y)+[39]*sin(0.5*pi*x)*sin(0.5*3*pi*y)+[40]*sin(0.5*pi*x)*cos(0.5*3*pi*y)+[41]*cos(0.5*pi*x)*sin(0.5*3*pi*y)+[42]*cos(0.5*pi*x)*cos(0.5*3*pi*y)",-1,1,-1,1);
	Fouri->SetTitle("Fourier series");
	cout<<"fourier function fit"<<endl;
	chi2->Fit("Fouri","0");
	Fouri->Write("",TObject::kOverwrite);

	TBox *minpoint;

	if(contfcn!=NULL){
		contfcn->GetMinimumXY(xx,yy);
		lambda_theta[file][trig][pt][frame] = xx;
		lambda_phi[file][trig][pt][frame] = yy;

		double zlevel = contfcn->GetMinimum()+0.5;
		contfcn->SetContour(1,&zlevel);		

		contfcn->Draw("CONT Z LIST");
		canvas[file][trig][pt][frame]->Update();
		TObjArray *smallconts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
		TList* smallcontlevel = NULL;
		TGraph* smallgc = NULL;
		smallcontlevel = (TList*)smallconts->At(0);	
		smallgc = (TGraph*)smallcontlevel->First();
		if(smallgc==NULL) {
			cout<<"contour of minL + 0.5 is NULL"<<endl;
			return;
		}
//		mlefcn->Draw("COLZ");
//		gc2->Draw("same");
		minpoint = new TBox(xx-0.02,yy-0.02,xx+0.02,yy+0.02);
		minpoint->SetFillColor(kRed);
		minpoint->SetLineColor(kRed);
		minpoint->Draw("same");
		xvector.clear();
		yvector.clear();
		while(1){
		for(int i=0;i<smallgc->GetN();i++){
			smallgc->GetPoint(i,xcont,ycont);
			xvector.push_back(xcont);
			yvector.push_back(ycont);
		}
		xmin = getminimum(xvector);
		xmax = getmaximum(xvector);
		ymin = getminimum(yvector);
		ymax = getmaximum(yvector);	
		if(xx<xmin || xx>xmax || yy<ymin || yy>ymax){
			smallgc = (TGraph*)smallcontlevel->After(smallgc);
			xvector.clear();
			yvector.clear();
		}
		else break;
		}
		smallgc->SetName("small_cont_min_half");
		smallgc->SetTitle("contour of minL+0.5");
		smallgc->SetLineWidth(1);
		smallgc->Write("",TObject::kOverwrite);	
		smallgc->Draw("same");

		lambda_theta_err[file][trig][pt][frame] = 0.5*(xmax-xmin);
		lambda_phi_err[file][trig][pt][frame] = 0.5*(ymax-ymin);
}
/*
	if(xmin>xx || xmax<xx || ymin>yy || ymax<yy){

		cout<<"update values"<<"xmin="<<xmin<<"xmax="<<xmax<<"ymin="<<ymin<<"ymax="<<ymax<<endl;

		xvector.clear();
		yvector.clear();
		contlevel = NULL;
		gc = NULL;

		contlevel = (TList*)conts->At(1);	
		gc = (TGraph*)contlevel->First();
		gc->SetName("cont_min_half");
		gc->SetTitle("contour of minL+0.5");
		gc->Write("",TObject::kOverwrite);	

		for(int i=0;i<gc->GetN();i++){
			gc->GetPoint(i,xcont,ycont);
			xvector.push_back(xcont);
			yvector.push_back(ycont);
		}
		xmin = getminimum(xvector);
		xmax = getmaximum(xvector);
		ymin = getminimum(yvector);
		ymax = getmaximum(yvector);	
		}
*/
	canvas[file][trig][pt][frame]->cd(3);
	contfcn->GetXaxis()->SetLimits(-1,1);
	contfcn->GetYaxis()->SetLimits(-1,1);
//	contfcn->GetYaxis()->SetRangeUser(-1,1);
	contfcn->Draw("colz2");
	smallgc->Draw("same");
	TGraphErrors* lambda_parameter = new TGraphErrors(1,&lambda_theta[file][trig][pt][frame],&lambda_phi[file][trig][pt][frame],&lambda_theta_err[file][trig][pt][frame],&lambda_phi_err[file][trig][pt][frame]);
	lambda_parameter->SetName("lambda");
	lambda_parameter->SetTitle("lambda parameter;#lambda_{#theta};#lambda_{#phi}");
	lambda_parameter->GetXaxis()->SetLimits(-1,1);
	lambda_parameter->GetYaxis()->SetRangeUser(-1,1);
	lambda_parameter->Write("",TObject::kOverwrite);
	lambda_parameter->Draw("same");

	canvas[file][trig][pt][frame]->cd(4);
	contfcn->Draw("colz");
	smallgc->Draw("same");
	lambda_parameter->Draw("same");

	canvas[file][trig][pt][frame]->cd(5);
	TH2F* difference = new TH2F("difference","difference;#lambda_{#theta};#lambda_{#phi}",xxbin,xmin,xmax,yybin,ymin,ymax);
	for(int i=1;i<xxbin;i++){
		for(int j=1;j<yybin;j++){
			double chi2_x = difference->GetXaxis()->GetBinCenter(i);
			double chi2_y = difference->GetYaxis()->GetBinCenter(j);
			double difference_xy = chi2->GetBinContent(chi2->GetXaxis()->FindBin(chi2_x),chi2->GetYaxis()->FindBin(chi2_y))-contfcn->Eval(chi2_x,chi2_y);
			difference->SetBinContent(i,j,difference_xy);
		}
	}		
	//difference->SetMaximum(5);
	//difference->SetMinimum(-5);
//	difference->GetXaxis()->SetLimits(-1,1);
//	difference->GetYaxis()->SetLimits(-1,1);
	difference->Draw("colz2");	
	//		mlefcn->Draw("CONT2 Z colz2 same");
	smallgc->Draw("same");
	minpoint->Draw("same");
	difference->Write("",TObject::kOverwrite);

/*
	canvas[file][trig][pt][frame]->cd(4);
	TH2F* ratio = new TH2F("ratio","ratio;#lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1);
	for(int i=1;i<101;i++){
		for(int j=1;j<101;j++){
			double chi2_x = difference->GetXaxis()->GetBinCenter(i);
			double chi2_y = difference->GetYaxis()->GetBinCenter(j);
			double ratio_xy = contfcn->Eval(chi2_x,chi2_y)/chi2->GetBinContent(i,j);
			ratio->SetBinContent(i,j,ratio_xy);
		}
	}		
	ratio->Draw("colz2");	
	//		mlefcn->Draw("CONT2 Z colz2 same");
	gc2->Draw("same");
	minpoint->Draw("same");
	ratio->Write("",TObject::kOverwrite);		
*/
	cout<<"lambda_theta ="<<lambda_theta[file][trig][pt][frame]<<"+/-"<<lambda_theta_err[file][trig][pt][frame]<<";"<<endl;
	cout<<"lambda_phi ="<<lambda_phi[file][trig][pt][frame]<<"+/-"<<lambda_phi_err[file][trig][pt][frame]<<endl;

	if(contfcn!=NULL)canvas[file][trig][pt][frame]->SaveAs(Form("~/polresults/20160707/pdf/sys_%d/mle_%d_%d_%d_%d.pdf",file,file,trig,pt,frame));
	inputfile->Close();	
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
		for(int frame=0;frame<NFRAME;frame++){ 
			for(int phase=0;phase<NPHASE+1;phase++){
				for(int pt=0;pt<NPT;pt++) {
					//					x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.+0.2*trig; // shift a little bit
					x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.;
					x_err[pt] = (PtEdge[pt+1]-PtEdge[pt])/2.;
					theta[pt] = lambda_theta[file][trig][pt][frame];
					thetaerr[pt] = lambda_theta_err[file][trig][pt][frame];
					phi[pt] = lambda_phi[file][trig][pt][frame];
					phierr[pt] = lambda_phi_err[file][trig][pt][frame];

					if(phase==0){
						y[pt] = lambda_theta[file][trig][pt][frame];
						y_err[pt] = lambda_theta_err[file][trig][pt][frame];
						//						y_sys[pt] = systematic_error[trig][pt][frame][phase];
					}
					if(phase==1){
						y[pt] = lambda_phi[file][trig][pt][frame];
						y_err[pt] = lambda_phi_err[file][trig][pt][frame];
						//						y_sys[pt] = systematic_error[trig][pt][frame][phase];
					}
					if(phase==2){//calculate lambda_invariant and its statistic error
						y[pt] = (theta[pt]+3*phi[pt])/(1-phi[pt]);
						y_err[pt] = TMath::Sqrt(theta[pt]*theta[pt]/((1-phi[pt])*(1-phi[pt]))+phi[pt]*phi[pt]*TMath::Power((3+theta[pt]/((1-phi[pt])*(1-phi[pt]))),2)+2*(3+theta[pt])/TMath::Power((1-phi[pt]),3)*covariant[trig][pt][frame]);
						//						y_sys[pt] = y[pt]-(lambda_theta[0][trig][pt][frame]+3*lambda_phi[0][trig][pt][frame])/(1-lambda_phi[0][trig][pt][frame]);
						//							y_sys[pt] = y[pt]-(fitparameters[0][trig][pt][frame][0]+3*fitparameters[0][trig][pt][frame][2])/(1-fitparameters[0][trig][pt][frame][2]);
						//cout<<"invariant ==================="<<y[pt]<<"trig"<<trig<<"   "<<"frame"<<frame<<endl;
					}
				}
				lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
				lambda_parameters[trig][frame][phase]->SetName(Form("lambdas_%d_%d_%d",trig,frame,phase));
				lambda_parameters_sys[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_sys);
				lambda_parameters_sys[trig][frame][phase]->SetName(Form("lambdas_sys_%d_%d_%d",trig,frame,phase));
				if(sys3==100) {
					//			lambda_parameters[trig][frame][phase]->Write();		
					//			lambda_parameters_sys[trig][frame][phase]->Write();		
				}
			}
		}
	}
	drawlambdas(sys3,trig,0);
}

void drawlambdas(int sys,int trig, int drawoptions = 0){
	TCanvas* lambdacanvas; 
	lambdas = new TFile(Form("~/polresults/20160707/lambdas/lambdas_trig%d.root",trig),"recreate");
	lambdas->cd();
	if(drawoptions==1){
		lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
		lambdacanvas->Divide(3,2);
		for(int phase=0;phase<NPHASE+1;phase++){
			for(int frame=0;frame<NFRAME;frame++){ 
				lambdacanvas->cd(frame*3+phase+1);
				legend = new TLegend(0.7,0.7,0.89,0.89);
				legend->SetBorderSize(0);
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
				lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s frame; J/#psi p_{T};%s",phaseName[phase].Data(),frameName[frame].Data(),phaseName[phase].Data()));
				lambda_parameters[trig][frame][phase]->Draw("ap");
				lambda_parameters[trig][frame][phase]->Write();
				compare_lambda(trig,frame,phase);
				legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
				legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
				legend->Draw("same");
				lambda_parameters[trig][frame][phase]->Draw("psame");
				lambda_parameters_sys[trig][frame][phase]->Draw("psame[]");
			}
		}
		lambdacanvas->SaveAs(Form("~/polresultspdsf/20160707/figures/sys_%d/lambdas_comparison_all.pdf",sys));
		return;
	}

	lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
	lambdacanvas->Divide(3,2);
	for(int phase=0;phase<NPHASE+1;phase++){
		for(int frame=0;frame<NFRAME;frame++){ 
			lambdacanvas->cd(frame*3+phase+1);
			legend = new TLegend(0.7,0.7,0.89,0.89);
			legend->SetBorderSize(0);
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
			lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s frame; J/#psi p_{T};%s",phaseName[phase].Data(),frameName[frame].Data(),phaseName[phase].Data()));
			lambda_parameters[trig][frame][phase]->Draw("ap");
			lambda_parameters[trig][frame][phase]->Write();
			//				if(phase!=2)
			compare_lambda(trig,frame,phase);
			legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
			legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
			legend->Draw("same");
			lambda_parameters[trig][frame][phase]->Draw("psame");
			//				lambda_parameters_sys[trig][frame][phase]->Draw("samep[]");
		}
	}
	lambdacanvas->SaveAs(Form("~/polresultspdsf/20160707/figures/sys_%d/lambdas_comparison_%d.pdf",sys,trig));
}

void plotaxissetting(int trig,int frame,int phase){
	lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s; J/#psi p_{T};#lambda_{#theta}",trigName[trig].Data(),frameName[frame].Data()));
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

void compare_lambda(int trig,int frame, int phase){
	if(phase==0) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
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

void minimization(int ifile,int itrig,int ipt,int iframe){
	//	TF2* fourier = new TF2("fourier","1+[0]*y*y+[1]*(1-y*y)cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);
	//	TH2F* chi2 = (TH2F*)inputfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));	

	float deltax = 0.02,deltay=2*TMath::Pi()/100;

	float Hessian[2][2];
	float invHessian[2][2];
	float gradient[2];

	lambda[0][0] = lambda[0][1] = lambda[1][0] = lambda[1][1] = 0.;

	while(1){
		int xbin,ybin;
		double x=0.,y=0.;

		gradient[0] = (mlefcn->Eval(x+deltax,y)-mlefcn->Eval(x-deltax,y))/(2*deltax);
		gradient[1] = (mlefcn->Eval(x,y+deltay)-mlefcn->Eval(x,y-deltay))/(2*deltay);

		Hessian[0][0] = (mlefcn->Eval(x+deltax,y)-2*mlefcn->Eval(x,y)+mlefcn->Eval(x-deltax,y))/(deltax*deltax);
		Hessian[1][1] = (mlefcn->Eval(x,y+deltay)-2*mlefcn->Eval(x,y)+mlefcn->Eval(x,y-deltay))/(deltay*deltay);
		Hessian[0][1] = Hessian[1][0] = (mlefcn->Eval(x+deltax,y+deltay) + mlefcn->Eval(x,y) - mlefcn->Eval(x+deltax,y) - mlefcn->Eval(x,y+deltay))/(deltax*deltay);

		invHessian[0][0] = Hessian[1][1]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[0][1] = -1*Hessian[0][1]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[1][0] = -1*Hessian[1][0]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[1][1] = Hessian[0][0]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		//		cout<<"invHessian[0][0] = "<<invHessian[0][0]<<"lambda[0][1] = "<<invHessian[0][1]<<"lambda[1][0]"<<invHessian[1][0]<<"lambda[1][1] = "<<invHessian[1][1]<<endl;

		lambda[0][1] = lambda[0][0] - (invHessian[0][0]*gradient[0]+invHessian[0][1]*gradient[1]);
		lambda[1][1] = lambda[1][0] - (invHessian[1][0]*gradient[0]+invHessian[1][1]*gradient[1]);

		//		cout<<"lambda[0][0] = "<<lambda[0][0]<<"lambda[0][1] = "<<lambda[0][1]<<"lambda[1][0]"<<lambda[1][0]<<"lambda[1][1] = "<<lambda[1][1]<<endl;
		//		cout<<"min x bin = "<< chi2->GetXaxis()->FindBin(lambda[0][0])<<"min y bin ="<<chi2->GetYaxis()->FindBin(lambda[1][0])<<endl;
		if(mlefcn->Eval(lambda[0][1],lambda[1][1])-mlefcn->Eval(lambda[0][0],lambda[1][0])<1e-3) {
			cout<<"lambda[0][1] = "<<lambda[0][1]<<"lambda[1][1] = "<<lambda[1][1]<<endl;
			//		cout<<"min x bin = "<< chi2->GetXaxis()->FindBin(lambda[0][0])<<"min y bin ="<<chi2->GetYaxis()->FindBin(lambda[1][0])<<endl;
			break;
		}
		else{
			lambda[0][0] = lambda[0][1];
			lambda[1][0] = lambda[1][1];
		}
	}
}




