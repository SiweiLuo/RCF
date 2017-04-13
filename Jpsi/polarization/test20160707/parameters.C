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

#define NFILE 40
#define NTRIG 3
#define NPT 6
#define NPHASE 2 
#define NFRAME 2
#define NREBIN 4
#define XNBIN 10 // 40
#define YNBIN 10 // 40
#define CHIXBIN 100
#define CHIYBIN 100

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht0","ht1","ht2"};

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME];
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][2]; // pass , total ;
TH2F *chi2result[NFILE][NTRIG][NPT][NFRAME];

TFile *efficiencyfile[NFILE];
TFile *rawdatafile[NFILE][NTRIG];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker

TH2F* chi2[NFILE][NTRIG][NPT][NFRAME];
TFile* histograms;
TFile* histograms_copy;
TFile* chi2histograms[NFILE][NTRIG][NPT][NFRAME];

void parameters(int ifile = 0,int itrig = 0, int xsect = 11, int ysect = 11){
	histograms = new TFile("histograms20161102.root","read");
	histograms_copy = new TFile("templates20161102.root","read");
	correcteddata(ifile);
	for(int ipt=0;ipt<6;ipt++){
		for(int iframe=0;iframe<2;iframe++){
			chisquare(ifile,itrig,ipt,iframe,xsect,ysect);
		}
	}
}

void chisquare(int ifile = 0,int itrig = 0,int ipt = 0,int iframe = 0,int xsect = 11,int ysect = 11){
	//	if(rawdata2D[ifile][itrig][ipt][iframe]->GetMean()!=rawdata2D[ifile][itrig][ipt][iframe]->GetMean() || eff2D[ifile][itrig][ipt][iframe][2]->GetMean()!=eff2D[ifile][itrig][ipt][iframe][2]->GetMean()) return;
	TH2F* templates[100][100];
	chi2histograms[ifile][itrig][ipt][iframe] = new TFile(Form("~/polresults/20160707/rootfiles/chi2histograms_%d/chi2histograms_%d_%d_%d_%d_%d_%d.root",ifile,ifile,itrig,ipt,iframe,xsect,ysect),"recreate");
	chi2[ifile][itrig][ipt][iframe] = new TH2F(Form("chi2_%d_%d_%d_%d",ifile,itrig,ipt,iframe),"#chi^{2};#lambda_{#theta};#lambda_{#phi}",CHIXBIN,-1,1,CHIYBIN,-1,1);
	chi2[ifile][itrig][ipt][iframe]->Sumw2();

	TH2F* rawdata = rawdata2D[ifile][itrig][ipt][iframe]->Clone();
	rawdata->Sumw2();
	//	rawdata->Scale(1./rawdata->Integral());

	TH2F* eff = (TH2F*)eff2D[ifile][itrig][ipt][iframe][2]->Clone();
	eff->SetName(Form("eff_%d_%d_%d_%d",ifile,itrig,ipt,iframe));
	eff->Sumw2();
	eff->SetTitle("efficiency;#phi;cos#theta");
	chi2histograms[ifile][itrig][ipt][iframe]->cd();
	TH2F* minchi2;
	TH2F* TRUTH;

	double minX2=1e4;	

	for(int xnbin=xsect;xnbin<xsect+10;xnbin++){
		for(int ynbin=ysect;ynbin<ysect+10;ynbin++){
			double lambda1,lambda2,X2=0;
			lambda1 = chi2[ifile][itrig][ipt][iframe]->GetXaxis()->GetBinCenter(xnbin);
			lambda2 = chi2[ifile][itrig][ipt][iframe]->GetYaxis()->GetBinCenter(ynbin);
			templates[xnbin-1][ynbin-1] = (TH2F*)histograms->Get(Form("theta_%d_phi_%d",xnbin-1,ynbin-1));
			templates[xnbin-1][ynbin-1]->Sumw2();
			templates[xnbin-1][ynbin-1]->SetName(Form("template_%d_%d",xnbin-1,ynbin-1));
			templates[xnbin-1][ynbin-1]->Multiply(eff);
			templates[xnbin-1][ynbin-1]->Scale(rawdata->Integral()/templates[xnbin-1][ynbin-1]->Integral());

			int empty = 0;
			double sum_error = 0.;
			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					if(rawdata->GetBinError(xbin,ybin)>1e-4) {
						sum_error += rawdata->GetBinError(xbin,ybin);
					}
					if(rawdata->GetBinError(xbin,ybin)<1e-4) {
						empty++;
					}
				}
			}

			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					//					if(rawdata->GetBinError(xbin,ybin)==0) rawdata->SetBinError(xbin,ybin,1e6);
					//					if(rawdata->GetBinError(xbin,ybin)==0 || rawdata->GetBinContent(xbin,ybin)==0) continue;
					if(rawdata->GetBinError(xbin,ybin)<1e-4) rawdata->SetBinError(xbin,ybin,sum_error/empty);
					//					X2 += TMath::Power((templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)-rawdata->GetBinContent(xbin,ybin))/rawdata->GetBinError(xbin,ybin),2);
					X2 += templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)*templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)/(rawdata->GetBinError(xbin,ybin)*rawdata->GetBinError(xbin,ybin)+templates[xnbin-1][ynbin-1]->GetBinError(xbin,ybin)*templates[xnbin-1][ynbin-1]->GetBinError(xbin,ybin));
				}
			}
			chi2[ifile][itrig][ipt][iframe]->Fill(lambda1,lambda2,X2);
			if(X2<=minX2){
				minX2=X2;
				int xx=xnbin,yy=ynbin;
				minchi2 = (TH2F*)templates[xnbin-1][ynbin-1]->Clone();
				minchi2->Sumw2();
				minchi2->SetName(Form("template_file%d_trg%d_pt%d_frame%d_x%d_y%d",ifile,itrig,ipt,iframe,xx,yy));
				minchi2->SetTitle(Form("#lambda_{#theta}=%.2f #lambda_{#phi}=%.2f, #chi^{2}=%.2f",lambda1,lambda2,X2));
				TRUTH = (TH2F*)histograms_copy->Get(Form("theta_%d_phi_%d",xx-1,yy-1));
			}
		}
	}
	rawdata->Write();
	eff->Write();
	chi2[ifile][itrig][ipt][iframe]->Write();
	if(minchi2!=0x0){
		minchi2->Write();
		TRUTH->Write();
	}
}

void correcteddata(int file=0){
	int selectedfile,file;
	int rebinning = 1;

	if(file>=0 && file<=NFILE) selectedfile=file+1;
	else if(file==100) file=0,selectedfile=NFILE;
	else return;
	
	for(;file<selectedfile;file++){
		efficiencyfile[file] = new TFile(Form("rootfile/OutFile_sys%d.root",file),"read"); 
		for(int trig=0;trig<NTRIG;trig++){
			if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),file),"read");//marker
			else rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),0),"read");//marker
			for(int pt=0;pt<NPT;pt++){
				for(int frame=0;frame<NFRAME;frame++){
					//					canvas[trig][pt][frame] = new TCanvas(Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),3200,1600);
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
					//					canvas[trig][pt][frame]->Divide(4,3);
					//					canvas[trig][pt][frame]->cd(1);
					//					rawdata3D[file][trig][pt][frame][0]->Draw();

					//					canvas[trig][pt][frame]->cd(2);
					//					rawdata3D[file][trig][pt][frame][1]->Draw();

					//					canvas[trig][pt][frame]->cd(3);
					rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike_like");
					if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_hx");
					else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_cs");
					rawdata3D[file][trig][pt][frame][2]->Sumw2();
					rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);
					//					rawdata3D[file][trig][pt][frame][2]->Draw();
					int max = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
					int min = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	
					if(frame==0){
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->Sumw2();
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPt1");
						eff3D[file][trig][pt][frame][1]->Sumw2();	
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
					}
					else{	
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->Sumw2();
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPtCS1");
						eff3D[file][trig][pt][frame][1]->Sumw2();
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
					}

					rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);
					//					rawdata3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					//					rawdata3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));

					//					canvas[trig][pt][frame]->cd(4);
					rawdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
					rawdata2D[file][trig][pt][frame]->Sumw2();
					//					rawdata2D[file][trig][pt][frame]->SetName(Form("rawdata_%d_%s_pT%d_frame%d",file,trigName[trig].Data(),pt,frame));
					rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
					if(rebinning==1){
						rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);
						rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);	
					}
					//					rawdata2D[file][trig][pt][frame]->Draw("colz");

					//					canvas[trig][pt][frame]->cd(5);
					eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
					eff2D[file][trig][pt][frame][0]->Sumw2();
					eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
					if(rebinning==1){
						eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);
						eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);
					}

					eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
					eff2D[file][trig][pt][frame][1]->Sumw2();
					eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
					if(rebinning==1){
						eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);
						eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);
					}
					//					rebining(file,trig,pt,frame);

					eff2D[file][trig][pt][frame][2] = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
					eff2D[file][trig][pt][frame][2]->Sumw2();
					eff2D[file][trig][pt][frame][2]->SetName(Form("eff_ratio_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][2]->Divide(eff2D[file][trig][pt][frame][1]);
					//					eff2D[file][trig][pt][frame][2]->Draw("colz");					
				}
			}
		}
	}
}

