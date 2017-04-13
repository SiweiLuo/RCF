#include <algorithm>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

#define NFILE 40
#define NTRIG 3
#define NPT 6
#define NPHASE 2 // theta, phi
#define NFRAME 2
#define NPLOT 6
#define NREBIN 4
#define	nsample 200
#define XNBIN 10
#define YNBIN 10
#define CHIXBIN 10
#define CHIYBIN 10

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
TF2* sigma[100][100];
TH2F* random2[100][100];
TH2F* chi2;
TF2* fEff;

void Minimization(int uncertainty = 0){
	//	gStyle->SetPadRightMargin(0.2);
	outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	
	gStyle->SetOptStat(false);
	correcteddata(0);
	for(gfile=0,gtrig=0;gtrig<1;gtrig++){
		for(gpt=2;gpt<3;gpt++){
			for(gframe=0;gframe<1;gframe++){
				NumericalMinimization();
			}
		}
	}
	//simultaneousfit(0);
	//systematic(uncertainty);
	//	lambdaparameters(uncertainty);
}

Double_t angular(Double_t *x, Double_t *par){
	Double_t result = par[3]*(1+par[0]*x[1]*x[1]+par[1]*(1-x[1]*x[1])*TMath::Cos(2*x[0]));
	return result;
}

void calchi2(int ifile,int itrig,int ipt,int iframe){
	chi2 = new TH2F("chi2","chi2",CHIXBIN,-1,1,CHIYBIN,-1,1);
	chi2->Sumw2();
	TRandom3 rndeff;
	std::vector<double> chi2vector;

	sample2 = (TH2F*)rawdata2D[ifile][itrig][ipt][iframe][0]->Clone();
	sample2->SetName(Form("raw%d_%d_%d_%d",ifile,itrig,ipt,iframe));
	sample2->Scale(1./sample2->Integral());
	TCanvas* sampleness  = new TCanvas();
	sampleness->Divide(2,2);
	sampleness->cd(1);
	sample2->Draw("colz");
	cout<<"==========>"<<sample2->GetXaxis()->GetTitle()<<endl;

	for(int xnbin=1;xnbin<CHIXBIN+1;xnbin++){
		for(int ynbin=1;ynbin<CHIYBIN+1;ynbin++){
			double lambda_theta_toy = chi2->GetXaxis()->GetBinCenter(xnbin);
			double lambda_phi_toy = chi2->GetYaxis()->GetBinCenter(ynbin);
			double chisquare = 0.;
			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					if(sample2->GetBinContent(xbin,ybin)>=1e-3)chisquare += TMath::Power((random2[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin))/(sample2->GetBinError(xbin,ybin)),2);
					//                	chisquare +=TMath::Power((random2[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)-sample2->GetBinContent(xbin,ybin)),2);
				}
			}
			cout<<"==================? chisquare"<<chisquare<<endl;
			chi2vector.push_back(chisquare);	
			chi2->Fill(lambda_theta_toy,lambda_phi_toy,chisquare);
		}
	}
	//	getminimum(chi2vector);

	int xx=1,yy=1,zz;
	chi2->GetMinimumBin(xx,yy,zz);
	//		result->Fill(chi2->GetXaxis()->GetBinCenter(xx),chi2->GetYaxis()->GetBinCenter(yy));

	lambda_theta[ifile][itrig][ipt][iframe] = chi2->GetYaxis()->GetBinCenter(xx);
	lambda_phi[ifile][itrig][ipt][iframe] = chi2->GetXaxis()->GetBinCenter(yy);

	cout<<"=====================>"<<lambda_theta[ifile][itrig][ipt][iframe]<<endl;

	if(xx>=1 && yy>=1){
		sampleness->cd(2); 
		random2[xx-1][yy-1]->Draw("colz");
		//        sigma[xx-1][yy-1]->Draw("colz");
	}

	sampleness->cd(3);
	chi2->Draw("colz");

	sampleness->SaveAs(Form("figures/samples2_%d_%d_%d_%d.pdf",ifile,itrig,ipt,iframe));
	chi2->Clear();
}

void correcteddata(int sys2 = 0){
	int selectedfile,file;
	if(sys2>=0 && sys2<=NFILE) file=sys2,selectedfile=sys2+1;

	efficiencyfile[file] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/OutFile_cent_0_9_%d.root",file),"read"); 
	for(int trig=0;trig<NTRIG;trig++){//marker
		if((sys2>=20 && sys2<=23) || (sys2>=29 && sys2<=30) || sys2==35) rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,sys2),"read");//marker
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
					rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][0]->Get("hJpsiCosThetaPhiPtCS"));
					rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][0]->Get("hJpsiCosThetaPhiPtCSBG"));
					rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_cs");
					rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_cs");
				}
				rawdata3D[file][trig][pt][frame][0]->Draw();

				canvas[trig][pt][frame]->cd(2);
				rawdata3D[file][trig][pt][frame][1]->Draw();

				canvas[trig][pt][frame]->cd(3);
				rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike-like");
				if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_hx");
				else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_cs");
				rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);
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
				rawdata2D[file][trig][pt][frame][0]->SetName(Form("rawdata_%s_pT%d_%d",trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1]));
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

				if(pt==0) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf(",trigName[trig].Data(),frame));			
				else if(pt==NPT-1) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf)",trigName[trig].Data(),frame));			
				else canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf",trigName[trig].Data(),frame));			
			}
		}
	}
}

void compare_lambda(int trig,int frame, int variable){
	if(variable==0)	comparisonlambda[trig][variable+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_0_1eID",trigName[trig].Data()));
	if(variable==1) comparisonlambda[trig][variable+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_phi_0_1eID",trigName[trig].Data()));
	if(trig==0){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kRed);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kRed);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(24);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==1){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kBlue);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kBlue);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(26);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==2){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kBlack);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kBlack);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(27);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	if(variable!=2) comparisonlambda[trig][variable+3*frame]->Draw("samep");
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
			//			else continue;
			ymin = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(-TMath::Pi());
			ymax = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(TMath::Pi());
			zmin = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt]);//marker
			zmax = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt+1])-1;

			//		cout<<"string================="<<Form("hHt%dJpsiCosThetaPt%d",trig,pt)<<endl;	
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

			//	h1[1]->Draw("p");
			//	defaulthist[0]->Draw("samep");

			TGraphAsymmErrors* test = new TGraphAsymmErrors(h1[0],h1[1],"N");
			//	test->Print("all");	
			//	h1[0]->Draw("p");
			oldhist[trig][pt]->SetTitle("efficiency comparison");
			oldhist[trig][pt]->GetXaxis()->SetTitle("cos#theta");
			oldhist[trig][pt]->Draw("p");
			//			oldhist[trig][pt]->Print("all");
			test->Draw("samep");
			//			test->Print("all");

			//	oldhist[0]->Print("all");
			//			canvas->cd(2);
			//			defaulthist[1]->Draw("p");
			//			h1[1]->Draw("samep");

			canvas2->cd(7+pt);
			defaultrawdata[0] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPt"))->ProjectionX("rawdata",zmin,zmax);
			defaultrawdata[0]->SetTitle("Unlike Sign Comparison");
			//	defaultrawdata[0]->Draw("p");	
			h2[pt][0] = (TH1F*)rawdata3D[0][0][2][0][0]->ProjectionX("h4",ymin,ymax,zmin,zmax);
			h2[pt][0]->RebinX(4);
			h2[pt][0]->SetMarkerStyle(24);
			//	h2[0]->Draw("same");

			//canvas->cd(5);
			defaultrawdata[1] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPtBG"))->ProjectionX("rawdatabkg",1,8);
			defaultrawdata[1]->SetTitle("Like Sign Comparison");
			//	defaultrawdata[1]->Draw("p");	
			h2[pt][1] = (TH1F*)rawdata3D[0][0][2][0][1]->ProjectionX("h5",1,41,1,8);
			h2[pt][1]->RebinX(4);
			h2[pt][1]->SetMarkerStyle(24);
			//	h2[1]->Draw("same");
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

Double_t fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	Double_t chisquare;

	int u1 = (int)((par[1]+1)*50);
	int u2 = (int)((par[2]+1)*50);

	if(u1<0 || u2<0 || u1>=99 || u2>=99) return 1e9*(par[1]*par[1]+par[2]*par[2]);

	TF2* sigma1 = new TF2("sigma1","[0]*(1+[1]*y*y+[2]*(1-y*y)*TMath::Cos(2*x))",-TMath::Pi(),TMath::Pi(),-1,1);	
	sigma1->SetParameters(par[0],par[1],par[2]);
	random2[u1][u2] = new TH2F(Form("theta%dphi%d_%d_%d_%d_%d",u1,u2,gfile,gtrig,gpt,gframe),"random2",YNBIN,-TMath::Pi(),TMath::Pi(),XNBIN,-1,1);
	random2[u1][u2]->Sumw2();

	TRandom3 rndeff;
	if(eff2D[gfile][gtrig][gpt][gframe][2]==0x0 || rawdata2D[gfile][gtrig][gpt][gframe][0]==0x0) return;
	TH2F* effhist = (TH2F*)eff2D[gfile][gtrig][gpt][gframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",gfile,gtrig,gpt,gframe));
	double theta,phi;
	for(int i=0;i<10000;i++){
		sigma1->GetRandom2(phi,theta);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi),effhist->GetYaxis()->FindBin(theta));
		if(rndeff.Uniform()<eff) random2[u1][u2]->Fill(phi,theta);	
	}
	random2[u1][u2]->Scale(1./random2[u1][u2]->Integral());

	Int_t THETANBIN=10,PHINBIN=10;

	sample2 = (TH2F*)rawdata2D[gfile][gtrig][gpt][gframe][0]->Clone();
	sample2->SetName(Form("raw%d_%d_%d_%d",gfile,gtrig,gpt,gframe));
	sample2->Scale(1./sample2->Integral());

	for(int thetabin=1;thetabin<THETANBIN+1;thetabin++){
		for(int phibin=1;phibin<PHINBIN+1;phibin++){
			if(sample2->GetBinContent(thetabin,phibin)>=1e-6 && sample2->GetBinError(thetabin,phibin)>=1e-6) chisquare += TMath::Power((random2[u1][u2]->GetBinContent(thetabin,phibin)-sample2->GetBinContent(thetabin,phibin))/(sample2->GetBinError(thetabin,phibin)),2);	
		}
	}
	sigma1->Clear();

	f = chisquare;

	return chisquare;
}

Double_t func1(Double_t *x,Double_t *par)
{
	Double_t value1=par[0]*(1+par[1]*x[0]*x[0]);
	return value1;
}

Double_t func2(Double_t *y,Double_t *par)
{
	Double_t value2=par[3]*(1+2*par[2]*TMath::Cos(2*y[0])/(3+par[1]));
	return value2;
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

void lambdaparameters(int sys3){// plot and write lambda parameters 

	int selectedfile,file;
	if(sys3>=0 && sys3<=NFILE) file=sys3,selectedfile=sys3+1;
	else if(sys3==100) {
		file=0,selectedfile=NFILE;
		lambdas = new TFile("lambdas.root","recreate");
		lambdas->cd();
	}
	else return;

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT],y_sys[NPT],theta[NPT],thetaerr[NPT],phi[NPT],phierr[NPT];

	for(;file<selectedfile;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){// marker add cs frame 
				for(int phase=0;phase<NPHASE+1;phase++){
					for(int pt=0;pt<NPT;pt++) {
						x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.;
						x_err[pt] = (PtEdge[pt+1]-PtEdge[pt])/2.;
						theta[pt] = fitparameters[file][trig][pt][frame][0];
						thetaerr[pt] = fitparameters[file][trig][pt][frame][1];
						phi[pt] = fitparameters[file][trig][pt][frame][2];
						phierr[pt] = fitparameters[file][trig][pt][frame][3];
						if(phase!=2){
							if(fitdata[trig][pt][frame][phase]->GetEntries()!=0){
								y[pt] = fitparameters[file][trig][pt][frame][phase*2];
								y_err[pt] = fitparameters[file][trig][pt][frame][phase*2+1];
								y_sys[pt] = systematic_error[trig][pt][frame][phase];
							}
							else{
								y[pt] = -100.;
								y_err[pt] = 0.;
								y_sys[pt] = 0.;
							}
						}
						else{//calculate lambda_invariant and its statistic error
							y[pt] = (theta[pt]+3*phi[pt])/(1-phi[pt]);
							y_err[pt] = TMath::Sqrt(theta[pt]*theta[pt]/((1-phi[pt])*(1-phi[pt]))+phi[pt]*phi[pt]*TMath::Power((3+theta[pt]/((1-phi[pt])*(1-phi[pt]))),2)+2*(3+theta[pt])/TMath::Power((1-phi[pt]),3)*covariant[trig][pt][frame]);
							y_sys[pt] = y[pt]-(fitparameters[0][trig][pt][frame][0]+3*fitparameters[0][trig][pt][frame][2])/(1-fitparameters[0][trig][pt][frame][2]);
						}
					}
					lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
					lambda_parameters_sys[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_sys);
					if(sys3==100) {
						lambda_parameters[trig][frame][phase]->Write();		
						//			lambda_parameters_sys[trig][frame][phase]->Write();		
					}
				}
			}
		}
	}

	TCanvas* lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
	lambdacanvas->Divide(3,2);

	for(int frame=0;frame<NFRAME;frame++){// marker add cs frame 
		for(int phase=0;phase<NPHASE+1;phase++){
			for(int trig=0;trig<NTRIG;trig++){
				lambdacanvas->cd(frame*3+phase+1);
				legend = new TLegend(0.7,0.7,0.89,0.89);
				lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
				if(trig==0){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
				}
				if(trig==1){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
				}
				plotsetting();
				legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
				if(trig==0)	lambda_parameters[trig][frame][phase]->Draw("ap");
				else lambda_parameters[trig][frame][phase]->Draw("samep");
				if(phase!=2)compare_lambda(trig,frame,phase);
				//				lambda_parameters_sys[trig][frame][phase]->Draw("samep[]");
			}
		}
	}
	lambdacanvas->SaveAs("figures/lambdas.pdf");
}

int NumericalMinimization()
{
	static Double_t vstart[4];
	static Double_t step[4];

	vstart[0] = 0.3;
	vstart[1] = 0.1;
	vstart[2] = 0.3;
	vstart[3] = 0.1;

	step[0] = 0.01;
	step[1] = 0.01;
	step[2] = 0.01;
	step[3] = 0.01;

	arglist[0] = 1;
	TMinuit *gMinuit = new TMinuit(4);
	gMinuit->SetFCN(fcn);
	gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

	gMinuit->mnparm(0,"norm1",vstart[0],step[0],0,0,ierflg);
	gMinuit->mnparm(1,"lambda1",vstart[1],step[1],0,0,ierflg);
	gMinuit->mnparm(2,"lambda2",vstart[2],step[2],0,0,ierflg);

	arglist[0] = 500.;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	return 0;
}

void plotsetting(){
	lambda_parameters[0][0][0]->SetTitle("#lambda_{#theta} in HX frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][1][0]->SetTitle("#lambda_{#theta} in CS frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][0][1]->SetTitle("#lambda_{#phi} in HX frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][1][1]->SetTitle("#lambda_{#phi} in CS frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][0][2]->SetTitle("#lambda_{inv} in HX frame;J/#psi p_{T};#lambda_{inv}");
	lambda_parameters[0][1][2]->SetTitle("#lambda_{inv} in CS frame;J/#psi p_{T};#lambda_{inv}");
}

void pull(int ifile,int itrig,int ipt,int iframe){
	TRandom3 rndeff;
	chi2 = new TH2F("chi2","chi2",CHIXBIN,-1,1,CHIYBIN,-1,1);
	TH2F* effhist = (TH2F*)eff2D[ifile][itrig][ipt][iframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",ifile,itrig,ipt,iframe));
	for(int xnbin=1;xnbin<CHIXBIN+1;xnbin++){
		for(int ynbin=1;ynbin<CHIYBIN+1;ynbin++){
			double lambda_theta_toy = chi2->GetXaxis()->GetBinCenter(xnbin);
			double lambda_phi_toy = chi2->GetYaxis()->GetBinCenter(ynbin);
			sigma[xnbin-1][ynbin-1] = new TF2(Form("sigma_%d_%d",xnbin,ynbin),"1+[0]*y*y+[1]*(1-y*y)*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);
			sigma[xnbin-1][ynbin-1]->SetParameters(lambda_theta_toy,lambda_phi_toy);
			random2[xnbin-1][ynbin-1] = new TH2F(Form("theta%dphi%d_%d_%d_%d_%d",xnbin,ynbin,ifile,itrig,ipt,iframe),"random2",YNBIN,-TMath::Pi(),TMath::Pi(),XNBIN,-1,1);
			random2[xnbin-1][ynbin-1]->Sumw2();
			for(int i=0;i<1000*nsample;i++){
				double costheta_toy,phi_toy;
				sigma[xnbin-1][ynbin-1]->GetRandom2(phi_toy,costheta_toy);
				//                double eff=fEff->Eval(costheta_toy,phi_toy);
				double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi_toy),effhist->GetYaxis()->FindBin(costheta_toy));// marker
				if(rndeff.Uniform()<eff) random2[xnbin-1][ynbin-1]->Fill(phi_toy,costheta_toy);
			}
			random2[xnbin-1][ynbin-1]->Scale(1./random2[xnbin-1][ynbin-1]->Integral());
		}
	}

	TCanvas* randomsss = new TCanvas();
	randomsss->cd();
	random2[5][5]->Draw("colz");
	//	effhist->Draw("colz");
	randomsss->SaveAs(Form("figures/rand_%d_%d_%d_%d.pdf",ifile,itrig,ipt,iframe));
	chi2->Clear();
}

void simultaneousfit(int sys1 = 0){// marker add frame 
	TCanvas* chi2resultplot[NFILE][NTRIG][NPT][NFRAME];
	// 	= new TCanvas("chi2resultplot","chi2resultplot",100,-1,1,100,-1,1);
	int file,selectedfile;
	if(sys1>=0 && sys1<=NFILE) file=sys1,selectedfile=sys1+1;
	else if(sys1==100) file=0,selectedfile=NFILE;
	for(;file<selectedfile;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){
				for(int pt=0;pt<NPT;pt++){
					chi2resultplot[file][trig][pt][frame] = new TCanvas(Form("chi2result_%d_%d_%d_%d",file,trig,pt,frame));
					chi2result[file][trig][pt][frame] = new TH2F(Form("chi2result_%d_%d_%d_%d",file,trig,pt,frame),"#chi^2; #lambda_{#theta};#lambda_{#phi}",100,-1,1,100,-1,1);
					pull(file,trig,pt,frame);
					calchi2(file,trig,pt,frame);
					chi2result[file][trig][pt][frame]->Fill(lambda_theta[file][trig][pt][frame],lambda_phi[file][trig][pt][frame]);
					chi2resultplot[file][trig][pt][frame]->cd();
					chi2result[file][trig][pt][frame]->Draw("colz");
					chi2resultplot[file][trig][pt][frame]->SaveAs(Form("figures/chi2resultplot_%d_%d_%d_%d.pdf",file,trig,pt,frame));
				}
			}
		}
	}
}

void systematic(int sys = 0){

	correcteddata(0);
	//crosscheck();
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
	systematic[0]->SaveAs("figures/systematic1.pdf");	

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

double RosenBrock(const double *xx)
{
	const Double_t p1 = xx[0];//lambda_theta
	const Double_t p2 = xx[1];//lambda_phi

	cout<<"p1=========="<<xx[0]<<endl;
	cout<<"p2=========="<<xx[1]<<endl;

	int u1 = (int)((p1+1)*50);
	int u2 = (int)((p2+1)*50);

	cout<<"u1=============="<<u1<<endl;
	cout<<"u2=============="<<u2<<endl;

	if(u1<0 || u2<0 || u1>=99 || u2>=99) return 1e9*(p1*p1+p2*p2);

	TF2* sigma1 = new TF2("sigma1","1+[0]*y*y+[1]*(1-y*y)*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);	
	sigma1->SetParameters(p1,p2);
	random2[u1][u2] = new TH2F(Form("theta%dphi%d_%d_%d_%d_%d",u1,u2,gfile,gtrig,gpt,gframe),"random2",YNBIN,-TMath::Pi(),TMath::Pi(),XNBIN,-1,1);
	random2[u1][u2]->Sumw2();

	TRandom3 rndeff;
	if(eff2D[gfile][gtrig][gpt][gframe][2]==0x0 || rawdata2D[gfile][gtrig][gpt][gframe][0]==0x0) return;
	TH2F* effhist = (TH2F*)eff2D[gfile][gtrig][gpt][gframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",gfile,gtrig,gpt,gframe));
	double theta,phi;
	for(int i=0;i<10000;i++){
		sigma1->GetRandom2(phi,theta);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi),effhist->GetYaxis()->FindBin(theta));
		if(rndeff.Uniform()<eff) random2[u1][u2]->Fill(phi,theta);	
	}
	random2[u1][u2]->Scale(1./random2[u1][u2]->Integral());

	Int_t THETANBIN=10,PHINBIN=10;
	double chisquare=0;

	sample2 = (TH2F*)rawdata2D[gfile][gtrig][gpt][gframe][0]->Clone();
	sample2->SetName(Form("raw%d_%d_%d_%d",gfile,gtrig,gpt,gframe));
	sample2->Scale(1./sample2->Integral());

	for(int thetabin=1;thetabin<THETANBIN+1;thetabin++){
		for(int phibin=1;phibin<PHINBIN+1;phibin++){
			if(sample2->GetBinContent(thetabin,phibin)>=1e-3) chisquare += TMath::Power((random2[u1][u2]->GetBinContent(thetabin,phibin)-sample2->GetBinContent(thetabin,phibin))/(sample2->GetBinError(thetabin,phibin)),2);	
		}
	}
	sigma1->Clear();
	return chisquare;
}

void simultaneousfit1(int sys1 = 0){// marker add frame 
	static Double_t vstart[NTRIG][NPT][5];
	static Double_t step[NTRIG][NPT][5];

	int file,selectedfile;
	if(sys1>=0 && sys1<=NFILE) file=sys1,selectedfile=sys1+1;
	else if(sys1==100) file=0,selectedfile=NFILE;
	for(;file<selectedfile;file++){ //marker
		TCanvas *canvas1[NTRIG];
		for(int trig=0;trig<NTRIG;trig++){
			canvas1[trig]	= new TCanvas(Form("canvas1%d",trig),Form("canvas1%d",trig),1200,800);
			canvas1[trig]->Divide(NPT,4);
			for(int frame=0;frame<NFRAME;frame++){
				for(int pt=0;pt<NPT;pt++){
					//	int frame = 0;

					fitdata[trig][pt][frame][0] = (TH1F*)costheta[file][trig][pt][frame][1]->Clone();
					fitdata[trig][pt][frame][0]->SetName("costheta_corrected_data");
					fitdata[trig][pt][frame][0]->Scale(1./fitdata[trig][pt][frame][0]->Integral());

					fitdata[trig][pt][frame][1] = (TH1F*)phi[file][trig][pt][frame][1]->Clone();
					fitdata[trig][pt][frame][1]->SetName("phi_corrected_data");
					fitdata[trig][pt][frame][1]->Scale(1./fitdata[trig][pt][frame][1]->Integral());

					for(int nbin=1;nbin<40/NREBIN+1;nbin++){
						z1[nbin-1] = fitdata[trig][pt][frame][0]->GetBinContent(nbin);
						errorz1[nbin-1] = fitdata[trig][pt][frame][0]->GetBinError(nbin);
						x[nbin-1] = fitdata[trig][pt][frame][0]->GetBinCenter(nbin);
						//								if(fabs(z2[nbin-1])>1e-6)empty[trig][pt][frame][eID-1][1]--;
					}

					for(int nbin=1;nbin<40/NREBIN+1;nbin++){
						z2[nbin-1] = fitdata[trig][pt][frame][1]->GetBinContent(nbin);
						errorz2[nbin-1] = fitdata[trig][pt][frame][1]->GetBinError(nbin);
						y[nbin-1] = fitdata[trig][pt][frame][1]->GetBinCenter(nbin);
						//								if(fabs(z2[nbin-1])>1e-6)empty[trig][pt][frame][eID-1][1]--;
					}

					vstart[trig][pt][0]=0.1;
					vstart[trig][pt][1]=0.1;
					vstart[trig][pt][2]=0.1;
					vstart[trig][pt][3]=-0.1;
					vstart[trig][pt][4]=0.001;

					step[trig][pt][0]=0.001;
					step[trig][pt][1]=0.001;
					step[trig][pt][2]=0.001;
					step[trig][pt][3]=0.001;
					step[trig][pt][4]=0.001;

					arglist[0] = 1;
					TMinuit *gMinuit = new TMinuit(5);
					gMinuit->SetFCN(fcn);
					gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

					gMinuit->mnparm(0,"norm1",vstart[trig][pt][0],step[trig][pt][0],0,0,ierflg);
					gMinuit->mnparm(1,"lambda1",vstart[trig][pt][1],step[trig][pt][1],0,0,ierflg);
					fitflag[trig][pt][frame][1]=ierflg;
					gMinuit->mnparm(2,"lambda2",vstart[trig][pt][2],step[trig][pt][2],0,0,ierflg);
					fitflag[trig][pt][frame][2]=ierflg;
					gMinuit->mnparm(3,"norm2",vstart[trig][pt][3],step[trig][pt][3],0,0,ierflg);

					arglist[0] = 500.;
					arglist[1] = 1.;
					gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
					Double_t amin,edm,errdef;
					Int_t nvpar,nparx,icstat;
					gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

					for(int var=0;var<4;var++)gMinuit->GetParameter(var,fParamVal[file][trig][pt][frame][var],fParamErr[file][trig][pt][frame][var]);
					lambda_theta[file][trig][pt][0] = fParamVal[file][trig][pt][frame][1];
					lambda_theta_err[file][trig][pt][0] = fParamErr[file][trig][pt][frame][1];
					lambda_phi[file][trig][pt][0] = fParamVal[file][trig][pt][frame][2];
					lambda_phi_err[file][trig][pt][0] = fParamErr[file][trig][pt][frame][2];

					gMinuit->mnemat(&matrix[0][0],4);
					covariant[trig][pt][frame] = matrix[1][2];

					fitparameters[file][trig][pt][frame][0] = fParamVal[file][trig][pt][frame][1];
					fitparameters[file][trig][pt][frame][1] = fParamErr[file][trig][pt][frame][1];
					fitparameters[file][trig][pt][frame][2] = fParamVal[file][trig][pt][frame][2];
					fitparameters[file][trig][pt][frame][3] = fParamErr[file][trig][pt][frame][2];

					fit1[pt][frame][0] = new TF1(Form("fit1_%d_%d",pt,0),func1,-1,1,2);
					fit2[pt][frame][0] = new TF1(Form("fit2_%d_%d",pt,0),func2,-TMath::Pi(),TMath::Pi(),4);
					fit1[pt][frame][1] = new TF1(Form("fit1_%d_%d",pt,1),func1,-1,1,2);
					fit2[pt][frame][1] = new TF1(Form("fit2_%d_%d",pt,1),func2,-TMath::Pi(),TMath::Pi(),4);

					fit1[pt][frame][0]->SetParameters(fParamVal[0][trig][pt][frame][0],fParamVal[0][trig][pt][frame][1]);
					fit2[pt][frame][0]->SetParameters(fParamVal[0][trig][pt][frame][0],fParamVal[0][trig][pt][frame][1],fParamVal[0][trig][pt][frame][2],fParamVal[0][trig][pt][frame][3]);
					canvas1[trig]->cd(pt+1);
					if(fitdata[trig][pt][frame][0]->GetEntries()!=0){
						fitdata[trig][pt][frame][0]->Draw();
						fit1[pt][frame][0]->Draw("same");	
					}
					canvas1[trig]->cd(pt+NPT+1);
					if(fitdata[trig][pt][frame][1]->GetEntries()!=0){
						fitdata[trig][pt][frame][1]->Draw();
						fit2[pt][frame][0]->Draw("same");
					}
					canvas1[trig]->SaveAs(Form("figures/simultaneousfit%d.pdf",trig));
				}
			}
		}
	}
}

