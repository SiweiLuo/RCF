#include "TFile.h"
#include "TCanvas.h"
#include "TPDF.h"

#define NTRIG 3
#define NPT 6
#define NFRAME 2
#define NPHASE 2
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

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString frameName[NFRAME] = {"HX","CS"};
TString phaseName[NPHASE+1] = {"#lambda_{#theta}","#lambda_{#phi}","#lambda_{inv}"};

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
TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda[NTRIG][NPT][NFRAME];
Double_t ltheta[NTRIG][NPT][NFRAME];
Double_t lphi[NTRIG][NPT][NFRAME];

void uncertainty(int file=0,int rebin=1){
	gStyle->SetPadRightMargin(0.2);	

	TDatime* time = new TDatime();
	TFile* lambdafile;
	TF2* likelihood[3][6][2];
	TGraph* contour[3][6][2];
	TCanvas* canvas = new TCanvas("canvas","canvas",4000,4000);
	canvas->Divide(6,6);

	for(int trig=0;trig<3;trig++){
		for(int pt=0;pt<6;pt++){
			for(int frame=0;frame<2;frame++){
				//				lambdafile = new TFile(Form("~/polresults/20160707/functional/sys_%d/functional_%d_%d_%d_%d.root",file,file,trig,pt,frame),"read");
				//				lambdafile = new TFile(Form("~/polresults/20160707/poisson/sys_%d/functional_%d_%d_%d_%d.root",file,file,trig,pt,frame),"read");
				lambdafile = new TFile(Form("~/polresults/20160707/poisson/sys_%d/functional_%d_%d_%d_%d.root",file,file,trig,pt,frame),"read");
				canvas->cd(pt+6*frame+12*trig+1);
				likelihood[trig][pt][frame] = (TF2*)lambdafile->Get("likelihood");
				contour[trig][pt][frame] = (TGraph*)lambdafile->Get("cont_min");
				lambda[trig][pt][frame] = (TGraphErrors*)lambdafile->Get("lambda");
				if(likelihood[trig][pt][frame]){
					likelihood[trig][pt][frame]->SetName(Form("likelihood_%d_%d_%d_%d",file,trig,pt,frame));
					likelihood[trig][pt][frame]->Draw("colz");
				}
				if(likelihood[trig][pt][frame] && contour[trig][pt][frame] && lambda[trig][pt][frame]){
					contour[trig][pt][frame]->SetName(Form("contour_%d_%d_%d_%d",file,trig,pt,frame));
					contour[trig][pt][frame]->Draw("same");
					lambda[trig][pt][frame]->SetName(Form("lambda_%d_%d_%d_%d",file,trig,pt,frame));
					lambda[trig][pt][frame]->Draw("same");
				}
			}
		}
	}
	canvas->SaveAs(Form("~/polresults/20160707/poisson/pdf_0/fit%d.pdf",file));

	oldfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");
	poissonfile = new TFile("~/polresults/20160707/poisson/lambda_0/lambdas.root","read");
	splotfile = new TFile("~/polresults/20160707/splot/sys_0/lambdasvspt.root","read");
	fit1dfile = new TFile("~/polresults/20160707/fit1D/sys_0/lambdasvspt.root","read");
	TCanvas* lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas");
	lambdacanvas->Divide(3,2);
	TLegend* legend;

	TFile* lambdaroot = new TFile(Form("~/polresults/20160707/functional/sys_%d/lambdas.root",file),"recreate");
	double x[NPT]=NULL,y[NPT]=NULL,x_err[NPT]=NULL,y_err[NPT]=NULL;

	for(int trig=0;trig<NTRIG;trig++){
		for(int frame=0;frame<NFRAME;frame++){					
			for(int phase=0;phase<NPHASE+1;phase++){	
				lambdacanvas->cd(frame*3+phase+1);

				for(int pt=0;pt<NPT;pt++){
					x[pt] = 	0.5*(PtEdge[pt]+PtEdge[pt+1]); 					
					x_err[pt] = 0.5*(PtEdge[pt+1]-PtEdge[pt]);

					lambdafile = new TFile(Form("~/polresults/20160707/functional/functional_0_%d_%d_%d.root",trig,pt,frame),"read");
					lambda[trig][pt][frame] = (TGraphErrors*)lambdafile->Get("lambda");
					if(lambda[trig][pt][frame]==0x0){
						y[pt] = -100;
						y_err[pt] = 0.;
					}
					else{
						lambda[trig][pt][frame]->GetPoint(0,ltheta[trig][pt][frame],lphi[trig][pt][frame]);
						cout<<"ltheta = "<<ltheta[trig][pt][frame]<<"lphi = "<<lphi[trig][pt][frame]<<endl;
						//					x[pt] = ltheta[trig][pt][frame];
						if(phase==0)		y[pt] = ltheta[trig][pt][frame];
						else if(phase==1)	y[pt] = lphi[trig][pt][frame];
						else 				y[pt] = (ltheta[trig][pt][frame]+3*lphi[trig][pt][frame])/(1-lphi[trig][pt][frame]); 

						//x_err[pt] = lambda[trig][pt][frame]->GetErrorX(0);
						if(phase==0)		y_err[pt] = lambda[trig][pt][frame]->GetErrorX(0);
						else if(phase==1) 	y_err[pt] = lambda[trig][pt][frame]->GetErrorY(0);
					}
				}

//				cout<<Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame)<<endl;	
//				comparisonlambda[trig][3*frame] = (TGraphErrors*)oldfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
//				comparisonlambda[trig][1+3*frame] = (TGraphErrors*)oldfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));

				lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
				lambda_parameters[trig][frame][phase]->SetName(Form("lambdas_%d_%d_%d",trig,phase,frame));
				lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
				lambda_parameters[trig][frame][phase]->GetXaxis()->SetLimits(0,14);

				if(trig==0){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
					//					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
					//					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(20);
				}
				else if(trig==1){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
					//					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
					//					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(22);
				}
				else lambda_parameters[trig][frame][phase]->SetMarkerStyle(23);

				lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s frame; J/#psi p_{T};%s",phaseName[phase].Data(),frameName[frame].Data(),phaseName[phase].Data()));
				//				if(trig==2)lambda_parameters[trig][frame][phase]->Draw("ap");
				lambda_parameters[trig][frame][phase]->Draw("ap");
				compare_lambda(trig,frame,phase);
				if(phase==2){
					legend = new TLegend(0.55,0.6,0.89,0.89);
					legend->SetBorderSize(0);
					legend->SetFillColor(-1);
					legend->AddEntry(lambda_parameters[trig][0][0],Form("mle %s",trigName[trig].Data()),"p");	
					//					if(phase!=2)
					legend->AddEntry(poissonlambda[trig][3*frame],Form("poisson",trigName[trig].Data()),"p");	
					legend->AddEntry(comparisonlambda[trig][3*frame],"1D fit","p");	
					legend->AddEntry(splotlambda[trig][3*frame],Form("splot",trigName[trig].Data()),"p");	
					legend->Draw("same");
				}

				//				legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("1Dfit %s",trigName[trig].Data()),"p");
				lambda_parameters[trig][frame][phase]->Draw("psame");
				lambdaroot->cd();
				lambda_parameters[trig][frame][phase]->Write();
				Double_t ptx,lambdatheta,lambdaphi;
				for(int pt=0;pt<6;pt++){
					if(phase!=2)comparisonlambda[trig][phase+3*frame]->GetPoint(pt,ptx,lambdatheta);
//					comparisonlambda[trig][1+3*frame]->GetPoint(pt,ptx,lambdaphi);
				}
			}
		}
		lambdacanvas->SaveAs(Form("~/polresults/20160707/poisson/pdf_0/lambdacanvas_%d_%d.pdf",trig,time->GetDate()));
	}
}

void compare_lambda(int trig,int frame, int phase){
	/*	if(phase==0) {
		comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)oldfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
		poissonlambda[trig][phase+3*frame] = (TGraphErrors*)poissonfile->Get(Form("lambda_parameters_%d_%d_%d",trig,frame,phase));
		}
		if(phase==1) {
		comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)oldfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));
		poissonlambda[trig][phase+3*frame] = (TGraphErrors*)poissonfile->Get(Form("lambda_parameters_%d_%d_%d",trig,frame,phase));

		}
		if(phase==2) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)oldfile->Get(Form("plot_%s_inv_%d_1eID",trigName[trig].Data(),frame));
		*/

	if(phase!=2){
		//	if(phase==0){
		comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)fit1dfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));
		comparisonlambda[trig][phase+3*frame]->SetName(Form("lambdas1dfit_%d_%d_%d",trig,frame,phase));
		poissonlambda[trig][phase+3*frame] = (TGraphErrors*)poissonfile->Get(Form("lambda_parameters_%d_%d_%d",trig,frame,phase));
		splotlambda[trig][phase+3*frame] = (TGraphErrors*)splotfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));

		//		cout<<"    "<<comparisonlambda[trig][phase+3*frame]<<"       "<<poissonlambda[trig][phase+3*frame]<<"     "<<splotlambda[trig][phase+3*frame]<<endl;


		//	}
		/*
		   if(phase==1){
		   comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)fit1dfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));
		   poissonlambda[trig][phase+3*frame] = (TGraphErrors*)poissonfile->Get(Form("lambda_parameters_%d_%d_%d",trig,frame,phase));
		   splotlambda[trig][phase+3*frame] = (TGraphErrors*)splotfile->Get(Form("lambdas_%d_%d_%d",trig,frame,phase));
		   }
		   */
		if(trig==0){
			//			comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kRed);
			//			comparisonlambda[trig][phase+3*frame]->SetLineColor(kRed);
			comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(24);
			comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
			poissonlambda[trig][phase+3*frame]->SetMarkerStyle(29); 
			poissonlambda[trig][phase+3*frame]->SetMarkerSize(1.2); 
			splotlambda[trig][phase+3*frame]->SetMarkerStyle(30);
			splotlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
		}
		if(trig==1){
			//			comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlue);
			//			comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlue);
			comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(26);
			comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
			poissonlambda[trig][phase+3*frame]->SetMarkerStyle(29); 
			poissonlambda[trig][phase+3*frame]->SetMarkerSize(1.2); 
			splotlambda[trig][phase+3*frame]->SetMarkerStyle(30);
			splotlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
		}
		if(trig==2){
			//			comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlack);
			//			comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlack);
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





