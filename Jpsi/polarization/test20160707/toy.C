#include "TH2F.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TError.h"
#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>

#define NEXP 1000
#define	NEVENT 1000
#define XNBIN 5
#define YNBIN 5
#define CHIXBIN 100
#define CHIYBIN 100
#define CHI2DIFF 1.

TH2F* sample1;
TH2F* sample2;
TCanvas* canvas;
TF2* sigma1[CHIXBIN][CHIXBIN];
TH2F* random1[CHIXBIN][CHIXBIN];
TF2* sigma2;
TH2F* random2;
TH2F* chi2;
TH2F* chi2plus;
TH2F* result_mean;
TH2F* result_error1;
TH2F* result_error2;
TF2* fEff;
TF1* fitx;
TF1* fity;

TFile* effrootfile;
TH2F* effhist;


void pull(){
	TRandom3 rndeff;
	chi2 = new TH2F("chi2","chi2",CHIXBIN,-1,1,CHIYBIN,-1,1);
	for(int xnbin=1;xnbin<CHIXBIN+1;xnbin++){
		for(int ynbin=1;ynbin<CHIYBIN+1;ynbin++){
			double lambda_theta = chi2->GetXaxis()->GetBinCenter(xnbin);
			double lambda_phi = chi2->GetYaxis()->GetBinCenter(ynbin);
			sigma1[xnbin-1][ynbin-1] = new TF2(Form("sigma_%d_%d",xnbin,ynbin),"1+[0]*x*x+[1]*(1-x*x)*TMath::Cos(2*y)",-1,1,-TMath::Pi(),TMath::Pi());
			sigma1[xnbin-1][ynbin-1]->SetParameters(lambda_theta,lambda_phi);

			random1[xnbin-1][ynbin-1] = new TH2F(Form("random_%d_%d",xnbin,ynbin),"random",XNBIN,-1,1,YNBIN,-TMath::Pi(),TMath::Pi());
			random1[xnbin-1][ynbin-1]->Sumw2();
			for(int i=0;i<100*NEVENT;i++){
				double costheta,phi;
				sigma1[xnbin-1][ynbin-1]->GetRandom2(costheta,phi);
				double eff=fEff->Eval(costheta,phi);
				random1[xnbin-1][ynbin-1]->Fill(costheta,phi,eff);
			}
			random1[xnbin-1][ynbin-1]->Scale(1./random1[xnbin-1][ynbin-1]->Integral());
		}
	}
	chi2->Clear();
}

void truth(int seed=0){
	double lambda_theta=-0.2,lambda_phi=0.2,lambda_theta_phi=0.2;

	TRandom3 rnd;
	rnd.SetSeed(seed);
	sample1 = new TH2F("sample1","truth;cos#theta;#phi",XNBIN,-1,1,YNBIN,-TMath::Pi(),TMath::Pi());
	sample1->Sumw2();
	sample2 = new TH2F("sample2","truth;cos#theta;#phi",XNBIN,-1,1,YNBIN,-TMath::Pi(),TMath::Pi());
	sample2->Sumw2();
	TF2* input = new TF2("input","1+[0]*x*x+[1]*(1-x*x)*cos(2*y)+[2]*cos(y)*sin(2*acos(x))",-1,1,-TMath::Pi(),TMath::Pi());
	input->SetParameters(lambda_theta,lambda_phi,lambda_theta_phi);

	effrootfile = new TFile("~/polresultspdsf/20160707/rootfiles/file0_trig0_pt0_frame0_x31_y31.root","read");
	effhist = (TH2F*)effrootfile->Get("eff_2dhist_0_0_0_0");

	canvas->cd(4);
	input->Draw("COLZ");
	for(int i=0;i<NEVENT;i++){
		double costheta,phi;
		input->GetRandom2(costheta,phi);
		sample1->Fill(costheta,phi);
//		double eff=fEff->Eval(costheta,phi);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(costheta),effhist->GetYaxis()->FindBin(phi));
//		if(rnd.Uniform()<eff) sample2->Fill(costheta,phi);
		sample2->Fill(costheta,phi,eff);
	}
	sample2->Scale(1./sample2->Integral());

	canvas->cd(1);
	sample2->SetTitle("raw data");
	sample2->Draw("colz");
	
	effhist->Clear();
	effrootfile->Close();
}

void calchi2(int seed=0){
	chi2 = new TH2F("chi2","chi2",CHIXBIN,-1,1,CHIYBIN,-1,1);
	chi2plus = new TH2F("chi2plus","#chi^{2}+1",CHIXBIN,-1,1,CHIYBIN,-1,1);
	chi2->Sumw2();
	chi2plus->Sumw2();
	for(int xnbin=1;xnbin<CHIXBIN+1;xnbin++){
		for(int ynbin=1;ynbin<CHIYBIN+1;ynbin++){
			double lambda_theta = chi2->GetXaxis()->GetBinCenter(xnbin);
			double lambda_phi = chi2->GetYaxis()->GetBinCenter(ynbin);
			double chisquare = 0.;
			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					chisquare += TMath::Power((random1[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin))/(sample2->GetBinError(xbin,ybin)),2);
				}
			}
			chi2->Fill(lambda_theta,lambda_phi,chisquare);
		}
	}
	canvas->cd(3);
	chi2->SetTitle("#chi^{2}");
	chi2->Draw("colz");
	int xx=1,yy=1,zz;
	chi2->GetMinimumBin(xx,yy,zz);
	double minimum = chi2->GetBinContent(xx,yy);
	cout<<"minimum========="<<minimum<<endl;
	result_mean->Fill(chi2->GetXaxis()->GetBinCenter(xx),chi2->GetYaxis()->GetBinCenter(yy));
	std::vector<double> xvector,yvector,xminvector,yminvector;	

	for(int xnbin=1;xnbin<=CHIXBIN+1;xnbin++){
		for(int ynbin=1;ynbin<=CHIYBIN+1;ynbin++){
			if(chi2->GetBinContent(xnbin,ynbin)<=minimum+1){
				chi2plus->Fill(chi2->GetXaxis()->GetBinCenter(xnbin),chi2->GetYaxis()->GetBinCenter(ynbin),chi2->GetBinContent(xnbin,ynbin));
				xvector.push_back(chi2->GetXaxis()->GetBinCenter(xnbin));
				yvector.push_back(chi2->GetYaxis()->GetBinCenter(ynbin));
				if(xnbin==xx)yminvector.push_back(chi2->GetYaxis()->GetBinCenter(ynbin));
				if(ynbin==yy)xminvector.push_back(chi2->GetXaxis()->GetBinCenter(xnbin));
			}		
		}
		xvector.clear();
		yvector.clear();
	}

	//cout<<getmaximum(xvector)<<endl;
	//cout<<getminimum(yvector)<<endl;

	//	result_error1->Fill(0.5*(getmaximum(xminvector)-getminimum(xminvector)),0.5*(getmaximum(yminvector)-getminimum(yminvector)));
	//result_error2->Fill(0.5*(getmaximum(xvector)-getminimum(xvector)),0.5*(getmaximum(yvector)-getminimum(yvector)));	
	if(xx>=1&& yy>=1){
		canvas->cd(2);
		random1[xx-1][yy-1]->Draw("colz");

		canvas->cd(5);
		sigma1[xx-1][yy-1]->Draw("colz");
	}

	TFile* outputfile = new TFile("outputfile.root","recreate");
	outputfile->cd();
	chi2->Write();
	chi2plus->Write();
	if(!chi2plus){
		result_error2 = (TH2F*)chi2plus->Clone();
		result_error2->SetName("result_error2");
		result_error2->Write();
	}
	outputfile->Close();
	sample1->Clear();
	sample2->Clear();
	chi2->Clear();
	chi2plus->Clear();
}

void toy(){
	//gStyle->SetOptStat(false);
	gStyle->SetPadRightMargin(0.2);
	canvas = new TCanvas("canvas","canvas",1200,1200);
	canvas->Divide(3,4);

	fEff=new TF2("fEff","(0.5-0.5*x)*cos(2*y)*cos(2*y)",-1,1,-TMath::Pi(),TMath::Pi());
	fitx=new TF1("fitx","pol3",-1,1);
	fity=new TF1("fity","pol3",-1,1);

	result_mean=new TH2F("result_mean","",CHIXBIN,-1,1,CHIYBIN,-1,1);
	result_error1=new TH2F("result_error1","",CHIXBIN,0,1,CHIYBIN,0,1);
	result_error2=new TH2F("result_error2","",CHIXBIN,0,1,CHIYBIN,0,1);

	pull();
	for(int i=0;i<NEXP;i++) {
		truth(i);
		calchi2();
		if(i==NEXP-1){
			canvas->cd(7);	
			result_error2->Draw();
		}
	}

	canvas->cd(6);
	result_mean->SetTitle("Mean Values");
	result_mean->Draw("colz");
	canvas->cd(10);
	result_error1->SetTitle("Uncertainties");
	result_error1->Draw("colz");
	canvas->cd(11);
	result_error1->ProjectionX()->Draw();
	canvas->cd(12);
	result_error1->ProjectionY()->Draw();

	canvas->Print("fig.pdf");
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



