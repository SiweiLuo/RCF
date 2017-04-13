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

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TFile* rootfile;

TGraphErrors* lambda[NTRIG][NPT][NFRAME][4];//method 
TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE][4];//method
Double_t ltheta[NTRIG][NPT][NFRAME];
Double_t lphi[NTRIG][NPT][NFRAME];

void getlambdas(){

	TString dir,smethod;
	int file=0;
	double x[NPT]=NULL,y[NPT]=NULL,x_err[NPT]=NULL,y_err[NPT]=NULL;

	TFile* outputfile;

	for(int method=0;method<2;method++){
		if(method==0)	outputfile = new TFile("~/polresults/20160707/splot/sys_0/lambdasvspt.root","recreate");
		else 			outputfile = new TFile("~/polresults/20160707/fit1D/sys_0/lambdasvspt.root","recreate");
		for(int trig=0;trig<3;trig++){
			for(int frame=0;frame<2;frame++){
				for(int phase=0;phase<2;phase++){
					for(int pt=0;pt<6;pt++){
						if(method==0) dir = Form("~/polresults/20160707/splot/sys_0/functional_%d_%d_%d_%d.root",file,trig,pt,frame);
						if(method==1) dir = Form("~/polresults/20160707/fit1D/sys_0/functional_%d_%d_%d_%d.root",file,trig,pt,frame);
						rootfile = new TFile(dir,"read");
						lambda[trig][pt][frame][method] = (TGraphErrors*)rootfile->Get("lambda");	

						x[pt] = 	0.5*(PtEdge[pt]+PtEdge[pt+1])+0.2*(method+1); 					
						x_err[pt] = 0.5*(PtEdge[pt+1]-PtEdge[pt]);

						if(lambda[trig][pt][frame][method]==0x0){
							y[pt] = -100;
							y_err[pt] = 0.;
						}
						else{
							lambda[trig][pt][frame][method]->GetPoint(0,ltheta[trig][pt][frame],lphi[trig][pt][frame]);
							if(phase==0)		y[pt] = ltheta[trig][pt][frame];
							else if(phase==1)	y[pt] = lphi[trig][pt][frame];
							else 				y[pt] = (ltheta[trig][pt][frame]+3*lphi[trig][pt][frame])/(1-lphi[trig][pt][frame]); 

							if(phase==0)		y_err[pt] = lambda[trig][pt][frame][method]->GetErrorX(0);
							else if(phase==1) 	y_err[pt] = lambda[trig][pt][frame][method]->GetErrorY(0);
							else ; //marker
						}
					}
					outputfile->cd();
					lambda_parameters[trig][frame][phase][method] = new TGraphErrors(NPT,x,y,x_err,y_err);
					lambda_parameters[trig][frame][phase][method]->SetName(Form("lambdas_%d_%d_%d",trig,frame,phase));

					lambda_parameters[trig][frame][phase][method]->GetYaxis()->SetRangeUser(-2,2);
					lambda_parameters[trig][frame][phase][method]->Write();
				}
			}
		}
	}
}



