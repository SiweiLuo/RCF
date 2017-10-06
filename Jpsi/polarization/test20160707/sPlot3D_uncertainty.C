#include <algorithm>
#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TPDF.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <fstream>

#define NTRIG 4
#define NPT 6
#define NFRAME 2
#define NPHASE 3
#define NPLOT 6
#define NFILE 40
#define NREBIN 4

TFile* rootfile;
TFile* combinedfile;
TH2F* rawdata;
TH2F* eff;
TH2F* chi2;
TH2F* truth;
TH2F* bestfit;
TH2F* check1Dfit;

TH2F* template_theta;
TH2F* template_phi;
TH2F* template_theta_eff;
TH2F* template_phi_eff;

TH2F* truthcomparison;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];
TGraphErrors* poissonlambda[NTRIG][NPLOT];
TGraphErrors* splotlambda[NTRIG][NPLOT];

TString trigName[NTRIG] = {"MB","HT0","HT1","HT2"};
TString frameName[NFRAME] = {"HX","CS"};
TString phaseName[NPHASE+1] = {"#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}","#lambda_{inv}"};

TLegend* leg;

TFile *efficiencyfile[NFILE];
TH1F *eff1D[NFILE][NTRIG][NPT][NFRAME][3]; //pass  , total , ratio;
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][2]; // pass , total ;

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};
TFile* oldfile;
TFile* poissonfile;
TFile* splotfile;
TFile* fit1dfile;
TGraphErrors* lambda_parameters[NTRIG][NFRAME][4];// theta , phi , thetaphi, invariant
//TGraphErrors* lambda[NTRIG][NPT][NFRAME];
Double_t lambda[NTRIG][NPT][NFRAME][4][30];
Double_t lambda_err[NTRIG][NPT][NFRAME][4][30];
Double_t ltheta[NTRIG][NPT][NFRAME];
Double_t lphi[NTRIG][NPT][NFRAME];

Double_t uncertainty_theta_hx;
Double_t uncertainty_theta_cs;
Double_t uncertainty_phi_hx;
Double_t uncertainty_phi_cs;
Double_t uncertainty_thetaphi_hx;
Double_t uncertainty_thetaphi_cs;


std::vector<double> sysvector_theta_hx,sysvector_theta_cs,sysvector_phi_hx,sysvector_phi_cs,sysvector_thetaphi_hx,sysvector_thetaphi_cs,sysvector_inv_hx,sysvector_inv_cs;

Double_t variation[9][NTRIG][6][2][4];

TFile* outputfile;


void sPlot3D_uncertainty(int trig=2,int pt=3,int frame=0, int phase=0){
	//	trig=3;
	//	phase=3;
	//	pt=3;
	gStyle->SetOptStat("eMR");

	TString dir;

	Double_t par[NTRIG][NPT][NFRAME][NPHASE+1];
	Double_t parerr[NTRIG][NPT][NFRAME][NPHASE+1];
	Double_t x;
	Double_t x_err;
	Double_t y;
	Double_t y_err;
	outputfile = new TFile(Form("~/polresults/20160707/functional/splot_3D_20170928_2eid_inv30/outtest_%d_%d_%d_%d.root",trig,pt,frame,phase),"recreate");
	dir = Form("~/polresults/20160707/functional/splot_3D_20171005_2eid_inv30/splot_3D_%d_%d_%d_0_20171005.root",trig,pt,frame);	
	cout<<"raeding file =>"<<dir<<endl;
	rootfile = new TFile(dir,"read");
	if(rootfile==0) continue;
	if((TH3F*)rootfile->Get("hlambda")!=0x0){
		if(phase!=3)par[trig][pt][frame][phase] = ((TH3F*)rootfile->Get("hlambda"))->GetMean(phase+1);
		else par[trig][pt][frame][phase] = ((TH1F*)rootfile->Get("hlambda_inv"))->GetMean();
		cout<<" par[trig][pt][frame][phase] === "<<par[trig][pt][frame][phase]<<endl;
		parerr[trig][pt][frame][phase] = uncertainty(trig,pt,frame,phase);
	}
	else continue;
	x = 	0.5*(PtEdge[pt]+PtEdge[pt+1]); 					
	x_err = 0.5*(PtEdge[pt+1]-PtEdge[pt]);

	y =     par[trig][pt][frame][phase];
	y_err = parerr[trig][pt][frame][phase];
	cout<<"pt ===="<<pt<<" y ==="<< par[trig][pt][frame][phase]<<"    "<<parerr[trig][pt][frame][phase]<<endl;

	lambda_parameters[trig][frame][phase] = new TGraphErrors(1,&x,&y,&x_err,&y_err);
	lambda_parameters[trig][frame][phase]->SetName(Form("lambda_%d_%d_%d_%d",trig,pt,frame,phase));
	lambda_parameters[trig][frame][phase]->SetTitle(Form("lambda_%d_%d_%d_%d;p_{T} GeV/c",trig,pt,frame,phase));
	lambda_parameters[trig][frame][phase]->GetXaxis()->SetLimits(0,12);
	lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
	TFile* outputfile2 = new TFile("outtest_20171005.root","update");
	outputfile2->cd();
	lambda_parameters[trig][frame][phase]->Write("",TObject::kOverwrite);
	outputfile2->Close();
}

Double_t uncertainty(int trig=2,int pt=3,int frame=0,int phase=0){

	double sysarray[21];
	gStyle->SetPadRightMargin(0.2);	

	TDatime* time = new TDatime();
	TFile* lambdafile;
	TCanvas* canvas = new TCanvas("syscanvas","uncertainty",4000,4000);
	canvas->Divide(2,2);
	//	TFile* reference = new TFile(Form("~/polresults/20160707/functional/splot_3D_20170706/splot_3D_%d_%d_%d_0_20170706.root",trig,pt,frame),"read");
	//	TFile* reference = new TFile(Form("~/polresults/20160707/functional/splot_3D_20170928_2eid_inv30/splot_3D_%d_%d_%d_0_20170928.root",trig,pt,frame),"read");
	TFile* reference = new TFile(Form("~/polresults/20160707/functional/splot_3D_20171005_2eid_inv30/splot_3D_%d_%d_%d_0_20171005.root",trig,pt,frame),"read");

	for(int file=0;file<=20;file++){
		//		lambdafile = new TFile(Form("~/polresults/20160707/functional/splot_3D_20170501/splot_3D_%d_%d_%d_%d_20170501.root",trig,pt,frame,file),"read");
		lambdafile = new TFile(Form("~/polresults/20160707/functional/splot_3D_20171005_2eid_inv30/splot_3D_%d_%d_%d_%d_20171005.root",trig,pt,frame,file),"update");
		if((TH3F*)lambdafile->Get("hlambda")==0x0) continue;
		if(phase!=3)lambda[trig][pt][frame][phase][file] = ((TH3F*)lambdafile->Get("hlambda"))->GetMean(phase+1);
		else lambda[trig][pt][frame][phase][file] = ((TH1F*)lambdafile->Get("hlambda_inv"))->GetMean();
		canvas->cd(1);
		((TH2F*)lambdafile->Get("datahist"))->ProjectionX("data_px");
		data_px->SetTitle("data;#phi");
		data_px->SetMarkerStyle(30);
		data_px->Draw("same");
		((TH2F*)reference->Get("datahist"))->ProjectionX("ref_data_px");
		ref_data_px->SetMarkerStyle(29);
		ref_data_px->Draw("same");

		TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.9);
		leg1->AddEntry(data_px,"systematic","lep");
		leg1->AddEntry(ref_data_px,"reference","lep");
		leg1->Draw("same");	

		canvas->cd(2);
		((TH2F*)lambdafile->Get("datahist"))->ProjectionY("data_py");
		data_py->SetTitle("data;cos#theta");
		data_py->SetMarkerStyle(30);
		data_py->Draw();
		((TH2F*)reference->Get("datahist"))->ProjectionY("ref_data_py");
		ref_data_py->SetMarkerStyle(29);
		ref_data_py->Draw("same");
		canvas->cd(3);
		((TH2F*)lambdafile->Get("effhist"))->ProjectionX("eff_px");
		eff_px->SetTitle("efficiency;#phi");
		eff_px->SetMarkerStyle(30);
		eff_px->Draw();
		((TH2F*)reference->Get("effhist"))->ProjectionX("ref_eff_px");
		ref_eff_px->SetMarkerStyle(29);
		ref_eff_px->Draw("same");
		canvas->cd(4);
		((TH2F*)lambdafile->Get("effhist"))->ProjectionY("eff_py");
		eff_py->SetTitle("efficiency;cos#theta");
		eff_py->SetMarkerStyle(30);
		eff_py->Draw();
		((TH2F*)reference->Get("effhist"))->ProjectionY("ref_eff_py");
		ref_eff_py->SetMarkerStyle(29);
		ref_eff_py->Draw("same");
		canvas->Write("",TObject::kOverwrite);
		canvas->SaveAs(Form("~/polresults/20160707/functional/splot_3D_20171005_2eid_inv30/pdf/uncertainty_%d_%d_%d_%d_20171005.pdf",trig,pt,frame,file));

		cout<<"     sys ====> "<<file<<"     delta lambda ="<<lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]<<endl;	
		sysarray[file] = lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0];
		cout<<" sysarray = "<<sysarray[file]<<endl;
	
	
	
	}

	ofstream unfile;
	unfile.open("sPlot3D_uncertainty.txt",std::ofstream::out | std::ofstream::app);
//	if(unfile.is_open()){
		for(int i=0;i<20;i++){
			unfile<<sysarray[i]<<"    ";
		}
			unfile<<'\n';
	//	unfile<<'\n'<<endl;
//	}
		unfile.close();	

	Double_t ptsmearing[NTRIG][6][2][3];
	int file=1;
	int part=0;
	/*
	   while(file<3){
	   if(frame==0 && phase==0)sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	   if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
	   if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	   if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	   if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	   if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));

	   fill_array(trig,pt,frame,phase,file);
	   file++;
	   }
	   component(part,trig,pt,frame,phase,0);
	   ptsmearing[trig][pt][0][0] = getaverage(sysvector_theta_hx);// marker getmaximum()
	   ptsmearing[trig][pt][1][0] = getaverage(sysvector_theta_cs);// marker getmaximum()
	   ptsmearing[trig][pt][0][1] = getaverage(sysvector_phi_hx);// marker getmaximum()
	   ptsmearing[trig][pt][1][1] = getaverage(sysvector_phi_cs);// marker getmaximum()
	   ptsmearing[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	   ptsmearing[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	   clear_array();
	   */

	Double_t weight[NTRIG][6][2][3];
	while(file>=1 && file<3){
		if(frame==0 && phase==0){
			sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
			//			cout<" fabs(lambda[trig][pt][frame][phase][file] = "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		}	
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		fill_array(trig,pt,frame,phase,file);
		file++;
	}
	component(1,trig,pt,frame,phase,1);
	/*	weight[trig][pt][0][0] = getaverage(sysvector_theta_hx);
		weight[trig][pt][1][0] = getaverage(sysvector_theta_cs);
		weight[trig][pt][0][1] = getaverage(sysvector_phi_hx);
		weight[trig][pt][1][1] = getaverage(sysvector_phi_cs);
		weight[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
		weight[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
		*/
	clear_array();

	Double_t nhitsfit[NTRIG][6][2][3];
	for(int file=5;file<9;file++){
		//	while(file>=12 && file<16){
		if(frame==0 && phase==0)sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		cout<<" part ====2 "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		fill_array(trig,pt,frame,phase,file);
	}
	component(2,trig,pt,frame,phase,1);
	if(frame==0 && phase==0)nhitsfit[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
	if(frame==1 && phase==0)nhitsfit[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
	if(frame==0 && phase==1)nhitsfit[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
	if(frame==1 && phase==1)nhitsfit[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
	if(frame==0 && phase==2)nhitsfit[trig][pt][0][2] = getmaximum(sysvector_thetaphi_hx);
	if(frame==1 && phase==2)nhitsfit[trig][pt][1][2] = getmaximum(sysvector_thetaphi_cs);

	clear_array();

	Double_t nsigma[NTRIG][6][2][3];
	for(int file=16;file<20;file++){
		if(frame==0 && phase==0)sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		fill_array(trig,pt,frame,phase,file);
	}
	component(3,trig,pt,frame,phase,0);
	nsigma[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	nsigma[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	nsigma[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	nsigma[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	nsigma[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	nsigma[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	clear_array();

	Double_t dedx[NTRIG][6][2][3];
	for(int file=3;file<4;file++){
		if(frame==0 && phase==0){
			sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
			//		cout<" fabs(lambda[trig][pt][frame][phase][file] = "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		}
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));

	}
	component(4,trig,pt,frame,phase,1);

	dedx[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	dedx[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	dedx[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	dedx[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	dedx[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	dedx[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	clear_array();

	Double_t Dca[NTRIG][6][2][3];
	for(int file=4;file<5;file++){
		if(frame==0 && phase==0){
			sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
			//		cout<" fabs(lambda[trig][pt][frame][phase][file] = "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		}
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	}
	component(5,trig,pt,frame,phase,1);
	Dca[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	Dca[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	Dca[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	Dca[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	Dca[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	Dca[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	clear_array();

	Double_t dsmadc[NTRIG][6][2][3];
	for(int file=9;file<11;file++){
		if(frame==0 && phase==0){
			sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
			//		cout<" fabs(lambda[trig][pt][frame][phase][file] = "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		}
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	}
	component(6,trig,pt,frame,phase,1);
	dsmadc[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	dsmadc[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	dsmadc[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	dsmadc[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	dsmadc[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	dsmadc[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	clear_array();

	Double_t beta[NTRIG][6][2][3];
	for(int file=11;file<15;file++){
		if(frame==0 && phase==0){
			sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
			//		cout<" fabs(lambda[trig][pt][frame][phase][file] = "<<fabs(lambda[trig][pt][frame][phase][file])<<endl;
		}
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	}
	component(7,trig,pt,frame,phase,1);
	if(frame==0 && phase==0)beta[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
	if(frame==1 && phase==0)beta[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
	if(frame==0 && phase==1)beta[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
	if(frame==1 && phase==1)beta[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
	if(frame==0 && phase==2)beta[trig][pt][0][2] = getmaximum(sysvector_thetaphi_hx);
	if(frame==1 && phase==2)beta[trig][pt][1][2] = getmaximum(sysvector_thetaphi_cs);
	clear_array();

	Double_t poe[NTRIG][6][2][3];
	for(int file=15;file<16;file++){
		if(frame==0 && phase==0)sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
	}
	component(8,trig,pt,frame,phase,1);
	poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
	poe[trig][pt][1][0] = getaverage(sysvector_theta_cs);
	poe[trig][pt][0][1] = getaverage(sysvector_phi_hx);
	poe[trig][pt][1][1] = getaverage(sysvector_phi_cs);
	poe[trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
	poe[trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
	clear_array();	

	/*	cout<<"uncertainty ============"<<fabs(ptsmearing[trig][pt][frame][phase])+fabs(weight[trig][pt][frame][phase])+fabs(nhitsfit[trig][pt][frame][phase])+fabs(nsigma[trig][pt][frame][phase])+fabs(dsmadc[trig][pt][frame][phase])+fabs(beta[trig][pt][frame][phase])+fabs(poe[trig][pt][frame][phase])+fabs(dedx[trig][pt][frame][phase])+fabs(Dca[trig][pt][frame][phase])<<endl;
		cout<<"ptsmaring ===="<<ptsmearing[trig][pt][frame][phase]<<endl;
		cout<<"weight ======="<<weight[trig][pt][frame][phase]<<endl;
		cout<<"nhitsfit ====="<<nhitsfit[trig][pt][frame][phase]<<endl;
		cout<<"nsigma========"<<nsigma[trig][pt][frame][phase]<<endl;
		cout<<"dsmadc========"<<dsmadc[trig][pt][frame][phase]<<endl;
		cout<<"beta=========="<<beta[trig][pt][frame][phase]<<endl;
		cout<<"poe==========="<<poe[trig][pt][frame][phase]<<endl;
		cout<<"dedx=========="<<dedx[trig][pt][frame][phase]<<endl;
		cout<<"dca =========="<<Dca[trig][pt][frame][phase]<<endl;
		*/
	//	return fabs(ptsmearing[trig][pt][frame][phase])+fabs(weight[trig][pt][frame][phase])+fabs(nhitsfit[trig][pt][frame][phase])+fabs(nsigma[trig][pt][frame][phase])+fabs(dsmadc[trig][pt][frame][phase])+fabs(beta[trig][pt][frame][phase])+fabs(poe[trig][pt][frame][phase])+fabs(dedx[trig][pt][frame][phase])+fabs(Dca[trig][pt][frame][phase]);
	Double_t result = fabs(variation[0][trig][pt][frame][phase])+fabs(variation[1][trig][pt][frame][phase])+fabs(variation[2][trig][pt][frame][phase])+fabs(variation[3][trig][pt][frame][phase])+fabs(variation[4][trig][pt][frame][phase])+fabs(variation[5][trig][pt][frame][phase])+fabs(variation[6][trig][pt][frame][phase])+fabs(variation[7][trig][pt][frame][phase])+fabs(variation[8][trig][pt][frame][phase]);

	outputfile->cd();
	Float_t systematic,central;
	TTree* tree;
	if((TTree*)outputfile->Get("tree")==0x0) tree = new TTree("tree","tree");
	else tree = (TTree*)outputfile->Get("tree");
	//	tree->Branch("systematic",&systematic,"systematic/F");
	//	TBranch *branch = tree->Branch(Form("systematic_%d_%d_%d_%d",trig,pt,frame,phase),&systematic,"systematic/F");
	TBranch *branch = tree->Branch("systematic",&systematic,"systematic/F");
	//	TBranch *branch2 = tree->Branch("central",&central,"central/F");
	for(int i=0;i<9;i++){
		tree->GetEntry(i);
		if(i==0){
			systematic = result;
		}
		else {
			systematic = variation[i-1][trig][pt][frame][phase];
			cout<<" variation ====="<<variation[i-1][trig][pt][frame][phase]<<"  i-1 === "<<i-1<<endl;
		}
		tree->Fill();
	}
	tree->Write("",TObject::kOverwrite);
	cout<<" result ==="<<result<<endl;
	return result;
	}

	void compare_lambda(int trig,int frame, int phase){
		if(phase!=2){
			comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)fit1dfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));
			comparisonlambda[trig][phase+3*frame]->SetName(Form("lambdas1dfit_%d_%d_%d",trig,frame,phase));
			poissonlambda[trig][phase+3*frame] = (TGraphErrors*)poissonfile->Get(Form("lambda_parameters_%d_%d_%d",trig,frame,phase));
			splotlambda[trig][phase+3*frame] = (TGraphErrors*)splotfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));

			if(trig==0){
				comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(24);
				comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
				poissonlambda[trig][phase+3*frame]->SetMarkerStyle(29); 
				poissonlambda[trig][phase+3*frame]->SetMarkerSize(1.2); 
				splotlambda[trig][phase+3*frame]->SetMarkerStyle(30);
				splotlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
			}
			if(trig==1){
				comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(26);
				comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
				poissonlambda[trig][phase+3*frame]->SetMarkerStyle(29); 
				poissonlambda[trig][phase+3*frame]->SetMarkerSize(1.2); 
				splotlambda[trig][phase+3*frame]->SetMarkerStyle(30);
				splotlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
			}
			if(trig==2){
				comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(29);
				comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
				poissonlambda[trig][phase+3*frame]->SetMarkerStyle(27); 
				poissonlambda[trig][phase+3*frame]->SetMarkerSize(1.2); 
				splotlambda[trig][phase+3*frame]->SetMarkerStyle(30);
				splotlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
			}
			comparisonlambda[trig][phase+3*frame]->Draw("samep");
			poissonlambda[trig][phase+3*frame]->Draw("samep"); 
			splotlambda[trig][phase+3*frame]->Draw("samep");
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

	void fill_array(int trig,int pt,int frame,int phase,int file){
		if(frame==0 && phase==0)sysvector_theta_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==0)sysvector_theta_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
		if(frame==0 && phase==1)sysvector_phi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==1)sysvector_phi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==2)sysvector_thetaphi_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==2)sysvector_thetaphi_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==0 && phase==3)sysvector_inv_hx.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));
		if(frame==1 && phase==3)sysvector_inv_cs.push_back(fabs(lambda[trig][pt][frame][phase][file]-lambda[trig][pt][frame][phase][0]));	
	}

	void component(int part,int trig,int pt,int frame,int phase,int mode){
		if(mode==0){
			if(frame==0 && phase==0)variation[part][trig][pt][0][0] = getaverage(sysvector_theta_hx);
			if(frame==1 && phase==0)variation[part][trig][pt][1][0] = getaverage(sysvector_theta_cs);
			if(frame==0 && phase==1)variation[part][trig][pt][0][1] = getaverage(sysvector_phi_hx);
			if(frame==1 && phase==1)variation[part][trig][pt][1][1] = getaverage(sysvector_phi_cs);
			if(frame==0 && phase==2)variation[part][trig][pt][0][2] = getaverage(sysvector_thetaphi_hx);
			if(frame==1 && phase==2)variation[part][trig][pt][1][2] = getaverage(sysvector_thetaphi_cs);
			if(frame==0 && phase==3)variation[part][trig][pt][0][3] = getaverage(sysvector_inv_hx);
			if(frame==1 && phase==3)variation[part][trig][pt][1][3] = getaverage(sysvector_inv_cs);
		}
		else{
			if(frame==0 && phase==0)variation[part][trig][pt][0][0] = getmaximum(sysvector_theta_hx);
			if(frame==1 && phase==0)variation[part][trig][pt][1][0] = getmaximum(sysvector_theta_cs);
			if(frame==0 && phase==1)variation[part][trig][pt][0][1] = getmaximum(sysvector_phi_hx);
			if(frame==1 && phase==1)variation[part][trig][pt][1][1] = getmaximum(sysvector_phi_cs);
			if(frame==0 && phase==2)variation[part][trig][pt][0][2] = getmaximum(sysvector_thetaphi_hx);
			if(frame==1 && phase==2)variation[part][trig][pt][1][2] = getmaximum(sysvector_thetaphi_cs);
			if(frame==0 && phase==3)variation[part][trig][pt][0][3] = getmaximum(sysvector_inv_hx);
			if(frame==1 && phase==3)variation[part][trig][pt][1][3] = getmaximum(sysvector_inv_cs);
		}
	}

	void clear_array(){
		sysvector_theta_hx.clear();
		sysvector_theta_cs.clear();
		sysvector_phi_hx.clear();
		sysvector_phi_cs.clear();
		sysvector_thetaphi_hx.clear();
		sysvector_thetaphi_cs.clear();
		sysvector_inv_hx.clear();
		sysvector_inv_cs.clear();
	}
