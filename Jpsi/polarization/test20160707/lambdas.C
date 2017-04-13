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
TCanvas *fitlambda[NFILE][NTRIG][NPT][NFRAME];
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
//TH2F* chi2_smooth;
TH2F* result_mean;
TH2F* result_sigma;

double lambda[2][2];
TGraphErrors* combinedlambda[3][3][2];

void lambdas(){
	gStyle->SetOptStat(0);
	outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	

	int file=0; //marker for loop 
	std::vector<double> sysvector_theta_hx,sysvector_theta_cs,sysvector_phi_hx,sysvector_phi_cs;

	//	for(;file<1;file++){ //marker for files number
	for(int trig=0;trig<NTRIG;trig++){
		for(int pt=0;pt<NPT;pt++){
			for(int frame=0;frame<NFRAME;frame++){
				scan(file,trig,pt,frame);
				//					for(int phase=0;phase<2;phase++){
				//						systematic_error[trig][pt][frame][phase] = 0;				
				//					}
			}
			/*				
							if(file>=1 && file<12){
							sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
							sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);// revise the frame 
							sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
							sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
							}
							if(file==11){
							ptsmearing[trig][pt][0][0] = getaverage(sysvector_theta_hx);
							ptsmearing[trig][pt][1][0] = getaverage(sysvector_theta_cs);
							ptsmearing[trig][pt][0][1] = getaverage(sysvector_phi_hx);// marker getmaximum()
							ptsmearing[trig][pt][1][1] = getaverage(sysvector_phi_cs);// marker getmaximum()
							sysvector_theta_hx.clear();
							sysvector_theta_cs.clear();
							sysvector_phi_hx.clear();
							sysvector_phi_cs.clear();
							}

							if(file>=12 && file<20){
							sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
							sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
							sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
							sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
							}
							if(file==19){
							weight[trig][pt][0][0] = getaverage(sysvector_theta_hx);
							weight[trig][pt][1][0] = getaverage(sysvector_theta_cs);
							weight[trig][pt][0][1] = getaverage(sysvector_phi_hx);
							weight[trig][pt][1][1] = getaverage(sysvector_phi_cs);
							sysvector_theta_hx.clear();
							sysvector_theta_cs.clear();
							sysvector_phi_hx.clear();
							sysvector_phi_cs.clear();
							}

							if(file>=20 && file<24){
							sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
							sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
							sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
							sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
							}
							if(file==23){
							nhitsfit[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
							nhitsfit[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
							nhitsfit[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
							nhitsfit[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
							sysvector_theta_hx.clear();
							sysvector_theta_cs.clear();
							sysvector_phi_hx.clear();
							sysvector_phi_cs.clear();
							}

							if(file>=24 && file<29){
							sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
							sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
							sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
							sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
							}
							if(file==28){
							nsigma[trig][pt][0][0] = getaverage(sysvector_theta_hx);
							nsigma[trig][pt][1][0] = getaverage(sysvector_theta_cs);
							nsigma[trig][pt][0][1] = getaverage(sysvector_phi_hx);
							nsigma[trig][pt][1][1] = getaverage(sysvector_phi_cs);
							sysvector_theta_hx.clear();
							sysvector_theta_cs.clear();
							sysvector_phi_hx.clear();
							sysvector_phi_cs.clear();
							}

							if(file>=29 && file<31){
							sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
	sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
	sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
	sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
}
if(file==30){
	dsmadc[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	dsmadc[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	dsmadc[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	dsmadc[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	sysvector_theta_hx.clear();
	sysvector_theta_cs.clear();
	sysvector_phi_hx.clear();
	sysvector_phi_cs.clear();
}

if(file>=31 && file<35){
	sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
	sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
	sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
	sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
}
if(file==34){
	beta[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
	beta[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
	beta[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
	beta[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
	sysvector_theta_hx.clear();
	sysvector_theta_cs.clear();
	sysvector_phi_hx.clear();
	sysvector_phi_cs.clear();
}
if(file==35){
	sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
	sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
	sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
	sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);

	poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	poe[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	poe[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	poe[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	sysvector_theta_hx.clear();
	sysvector_theta_cs.clear();
	sysvector_phi_hx.clear();
	sysvector_phi_cs.clear();
}

if(file>=36 && file<40){
	sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
	sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
	sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
	sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
}
if(file==39){
	pol[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	pol[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	pol[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	pol[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	sysvector_theta_hx.clear();
	sysvector_theta_cs.clear();
	sysvector_phi_hx.clear();
	sysvector_phi_cs.clear();
}
for(int phase=0;phase<2;phase++){
	//					systematic_error[trig][pt][frame][phase] = 0;
	//					systematic_error[trig][pt][frame][phase] = TMath::Sqrt(fabs(ptsmearing[trig][pt][frame][phase])*fabs(ptsmearing[trig][pt][frame][phase])+fabs(weight[trig][pt][frame][phase])*fabs(weight[trig][pt][frame][phase])+fabs(nhitsfit[trig][pt][frame][phase])*fabs(nhitsfit[trig][pt][frame][phase])+fabs(nsigma[trig][pt][frame][phase])*fabs(nsigma[trig][pt][frame][phase])+fabs(dsmadc[trig][pt][frame][phase])*fabs(dsmadc[trig][pt][frame][phase])+fabs(beta[trig][pt][frame][phase])*fabs(beta[trig][pt][frame][phase])+fabs(poe[trig][pt][frame][phase])*fabs(poe[trig][pt][frame][phase])+fabs(pol[trig][pt][frame][phase])*fabs(pol[trig][pt][frame][phase]));
}
*/
}	
lambdaparameters(file,trig);
}
combinelambdas(file);
//	}
}


void scan(int file = 0, int trig = 0, int pt = 0, int frame = 0){
	gStyle->SetPadRightMargin(0.2);	

	TFile* inputfile = new TFile(Form("~/polresultspdsf/20160707/rootcombined/lambda_file%d_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");
	//	TFile* inputfile = new TFile(Form("~/polresultspdsf/results_backup_20161130_2/rootcombined/lambda_file%d_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");
	double x,y,xerr,yerr;

	if((TGraphErrors*)inputfile->Get("lambda")!=0x0){
		TGraphErrors* get_lambda_parameter = (TGraphErrors*)inputfile->Get("lambda");
		get_lambda_parameter->GetPoint(0,x,y);
		lambda_theta[file][trig][pt][frame]=x;
		lambda_phi[file][trig][pt][frame]=y;
		lambda_theta_err[file][trig][pt][frame]=get_lambda_parameter->GetErrorX(0);
		lambda_phi_err[file][trig][pt][frame]=get_lambda_parameter->GetErrorY(0);
	}
	else{
		lambda_theta[file][trig][pt][frame] = -100;
		lambda_theta_err[file][trig][pt][frame] = 0.;
		lambda_phi[file][trig][pt][frame] = -100;
		lambda_phi_err[file][trig][pt][frame] = 0.;	
		return;
	}



	/*
	   chi2 = (TH2F*)inputfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));
	   TF2* mlefit = (TF2*)inputfile->Get("mlefcn");

	   if(chi2==0x0 || chi2->GetMean()!=chi2->GetMean() || chi2->GetMean()==0 ){
	   lambda_theta[file][trig][pt][frame] = -100;
	   lambda_theta_err[file][trig][pt][frame] = 0.;
	   lambda_phi[file][trig][pt][frame] = -100;
	   lambda_phi_err[file][trig][pt][frame] = 0.;	
	   return;
	   }

	   result_sigma = new TH2F("result_sigma","#sigma_{#chi^{2}};#lambda_{#theta};#lambda_{#phi}",CHIXBIN,-1,1,CHIYBIN,-1,1);
	   std::vector<double> lambda_theta_error,lambda_phi_error;

	   double xx,yy,zz;
	   mlefit->GetMinimumXY(xx,yy);
	   double minimum = mlefcn->GetMinimum();

	//	double minimum = chi2->GetBinContent(xx,yy);
	//	lambda_theta[file][trig][pt][frame] = chi2->GetXaxis()->GetBinCenter(xx);
	//	lambda_phi[file][trig][pt][frame] = chi2->GetYaxis()->GetBinCenter(yy);

	lambda_theta[file][trig][pt][frame] = xx;
	lambda_phi[file][trig][pt][frame] = yy;

	double contour[1];
	contour[0] = minimum +0.5;
	mlefit->SetContour(1,contour);


	canvas[file][trig][pt][frame] = new TCanvas(Form("mle_%d_%d_%d_%d",file,trig,pt,frame),Form("mle_%d_%d_%d_%d",file,trig,pt,frame),2600,600);
	canvas[file][trig][pt][frame]->Divide(4,1);
	canvas[file][trig][pt][frame]->cd(1);
	//	chi2->Draw("colz");
	//		chi2_smooth->Draw("colz");
	chi2->Draw("surf2");

	canvas[file][trig][pt][frame]->cd(2);
	//	result_sigma->SetTitle("contour minL+0.5");
	//	result_sigma->Draw("colz");
	mlefit->Draw("cont2 z");
	canvas[file][trig][pt][frame]->cd(3);

	//		result_sigma->ProjectionX("result_sigma_px")->Fit("gaus","Q");
	//	result_sigma_px->Fit("gaus","Q");
	//		lambda_theta_err[file][trig][pt][frame] = 2*gaus->GetParameter(2);
	//		result_sigma_px->SetTitle(Form("2#sigma_{#lambda_{#theta}} = %.2f",lambda_theta_err[file][trig][pt][frame]));
	canvas[file][trig][pt][frame]->cd(4);
	//		result_sigma->ProjectionY("result_sigma_py")->Fit("gaus","Q");
	//		lambda_phi_err[file][trig][pt][frame] = 2*gaus->GetParameter(2);
	//		result_sigma_py->SetTitle(Form("2#sigma_{#lambda_{#phi}} = %.2f",lambda_phi_err[file][trig][pt][frame]));
	cout<<"lambda_theta="<<lambda_theta[file][trig][pt][frame]<<"          "<<"lambda_phi="<<lambda_phi[file][trig][pt][frame]<<endl;
	//		canvas[file][trig][pt][frame]->SaveAs(Form("~/polresultspdsf/20160707/pdf/sys_%d/mle_%d_%d_%d_%d.pdf",file,file,trig,pt,frame));

	 */

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
	lambdas = new TFile(Form("~/polresultspdsf/20160707/lambdas/lambdas_trig%d.root",trig),"recreate");
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

void combinelambdas(int sys){
	TFile* inputfile[3];
	TGraphErrors* lam;
	TCanvas* combined = new TCanvas("combined","combined",2400,1200);
	combined->Divide(3,2);
	TLegend* leg[3][2];

	for(int phase=0;phase<3;phase++){
		for(int frame=0;frame<2;frame++){
			combined->cd(phase+frame*3+1);
			leg[phase][frame] = new TLegend(0.6,0.7,0.8,0.9);
			for(int trig=0;trig<3;trig++){
				inputfile[trig] = new TFile(Form("~/polresultspdsf/20160707/lambdas/lambdas_trig%d.root",trig),"read");
				combinedlambda[trig][phase][frame]=(TGraphErrors*)inputfile[trig]->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));
				if(phase==2)combinedlambda[trig][phase][frame]->GetYaxis()->SetRangeUser(-3,3);
				if(trig==0)combinedlambda[trig][phase][frame]->Draw("ap");
				else combinedlambda[trig][phase][frame]->Draw("psame");
				
				if(trig==0 && frame==0){
					combinedlambda[trig][phase][frame]->RemovePoint(5);
					combinedlambda[trig][phase][frame]->RemovePoint(4);
					combinedlambda[trig][phase][frame]->RemovePoint(0);
				}

				if(trig==1 && frame==0){
					combinedlambda[trig][phase][frame]->RemovePoint(5);
					combinedlambda[trig][phase][frame]->RemovePoint(4);
				}
				if(trig==2 && frame==0){
					combinedlambda[trig][phase][frame]->RemovePoint(2);
				}
				if(trig==0 && frame==1){
					combinedlambda[trig][phase][frame]->RemovePoint(5);
					combinedlambda[trig][phase][frame]->RemovePoint(1);
					combinedlambda[trig][phase][frame]->RemovePoint(0);
				}
				if(trig==1 && frame==1){
					combinedlambda[trig][phase][frame]->RemovePoint(5);
					combinedlambda[trig][phase][frame]->RemovePoint(2);
				}
				if(trig==2 && frame==1){
					combinedlambda[trig][phase][frame]->RemovePoint(2);
				}	
				leg[phase][frame]->AddEntry(combinedlambda[trig][phase][frame],Form("HT%d",trig),"lep");
			}
		leg[phase][frame]->Draw("same");
		}
	}
	combined->SaveAs(Form("~/polresultspdsf/20160707/figures/sys_%d/combinedlambdas.pdf",sys));
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

		cout<<"chi2 ="<<chi2<<endl;

		xbin = chi2->GetXaxis()->FindBin(lambda[0][0]);
		ybin = chi2->GetYaxis()->FindBin(lambda[1][0]);

		gradient[0] = (chi2->GetBinContent(xbin+1,ybin)-chi2->GetBinContent(xbin-1,ybin))/(2*deltax);
		gradient[1] = (chi2->GetBinContent(xbin,ybin+1)-chi2->GetBinContent(xbin,ybin-1))/(2*deltay);

		Hessian[0][0] = (chi2->GetBinContent(xbin+1,ybin)-2*chi2->GetBinContent(xbin,ybin)+chi2->GetBinContent(xbin-1,ybin))/(deltax*deltax);
		Hessian[1][1] = (chi2->GetBinContent(xbin,ybin+1)-2*chi2->GetBinContent(xbin,ybin)+chi2->GetBinContent(xbin,ybin-1))/(deltay*deltay);
		Hessian[0][1] = Hessian[1][0] = (chi2->GetBinContent(xbin+1,ybin+1) + chi2->GetBinContent(xbin,ybin) - chi2->GetBinContent(xbin+1,ybin) - chi2->GetBinContent(xbin,ybin+1))/(deltax*deltay);

		invHessian[0][0] = Hessian[1][1]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[0][1] = -1*Hessian[0][1]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[1][0] = -1*Hessian[1][0]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		invHessian[1][1] = Hessian[0][0]/(Hessian[0][0]*Hessian[1][1] - Hessian[0][1]*Hessian[1][0]);
		//		cout<<"invHessian[0][0] = "<<invHessian[0][0]<<"lambda[0][1] = "<<invHessian[0][1]<<"lambda[1][0]"<<invHessian[1][0]<<"lambda[1][1] = "<<invHessian[1][1]<<endl;

		lambda[0][1] = lambda[0][0] - (invHessian[0][0]*gradient[0]+invHessian[0][1]*gradient[1]);
		lambda[1][1] = lambda[1][0] - (invHessian[1][0]*gradient[0]+invHessian[1][1]*gradient[1]);

		//		cout<<"lambda[0][0] = "<<lambda[0][0]<<"lambda[0][1] = "<<lambda[0][1]<<"lambda[1][0]"<<lambda[1][0]<<"lambda[1][1] = "<<lambda[1][1]<<endl;
		//		cout<<"min x bin = "<< chi2->GetXaxis()->FindBin(lambda[0][0])<<"min y bin ="<<chi2->GetYaxis()->FindBin(lambda[1][0])<<endl;
		if(chi2->GetXaxis()->FindBin(lambda[0][0])==chi2->GetXaxis()->FindBin(lambda[0][1]) && chi2->GetYaxis()->FindBin(lambda[1][0])==chi2->GetYaxis()->FindBin(lambda[1][1])) {
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




