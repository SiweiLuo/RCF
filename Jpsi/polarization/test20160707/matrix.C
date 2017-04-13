#include <algorithm>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TMinuit.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFormula.h"

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
TF2* sigmab; // background pdf
TF2* sigmafcn;// generate mc data 
TH2F* effhist;
TH2F* datahist;
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

double matrix[10][10][2][2],invmatrix[10][10][2][2];


Double_t likelihoodfcn(Double_t *x, Double_t *par)
{
	Double_t result=0.;

	int nbinxeff = effhist->GetNbinsX();
	int nbinyeff = effhist->GetNbinsY();

	double costheta=0.,phi=0.,eff=0.;
	double si = 0.,xi = 0.;
	int Npositron = 0;

	// normalization of PDFs of signal and background

	sigma = new TF2("sigma","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",-TMath::Pi(),TMath::Pi(),-1,1);
	sigmab = new TF2("sigmab","[0]*(1+[1]*y**2+[2]*(1-y**2)*cos(2*x))",-TMath::Pi(),TMath::Pi(),-1,1);

	double scale = 0.,scaleb = 0.;

	sigma->SetParameter(0,1);
	sigma->SetParameter(1,par[0]);
	sigma->SetParameter(2,par[1]);
	
	sigmab->SetParameter(0,1);
	sigmab->SetParameter(1,par[2]);
	sigmab->SetParameter(2,par[3]);

	cout<<"par[0] = "<<par[0]<<"    "<<par[1]<<"    "<<par[2]<<"    "<<par[3]<<"    "<<sigma->GetParameter(1)<<"     "<<sigma->GetParameter(2)<<endl;

	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){

			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			si = sigma->Eval(phi,costheta);
			xi = sigmab->Eval(phi,costheta);
			eff = effhist->GetBinContent(i,j);
			
			if(si*eff<1e-16) continue;
			scale += si*eff;
			scaleb += xi*eff;
		}
	}

	sigma->SetParameter(1,par[0]);
	sigma->SetParameter(2,par[1]);
	sigma->SetParameter(0,1./scale);

	sigmab->SetParameter(1,par[2]);
	sigmab->SetParameter(2,par[3]);
	sigmab->SetParameter(0,1./scaleb);	


	double normal1=0.,normal2=0.;
	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			si = sigma->Eval(phi,costheta);
			xi = sigmab->Eval(phi,costheta);
			eff = effhist->GetBinContent(i,j);
			normal1 += si*eff;
			normal2 += xi*eff;
		}
	}
	cout<<"normalization ======"<<normal1<<"         "<<normal2<<endl;

//	cout<<"sigma integral = "<<sigma->Integral(-TMath::Pi(),TMath::Pi(),-1,1)<<"  ============"<<sigmab->Integral(-TMath::Pi(),TMath::Pi(),-1,1)<<endl;
// calculate invertible matrix V^-1 in sFit technique
	
	double sum = 0.,sum1 = 0.;
	double sWeight[10][10]=0.;

	invV(0,0) = invV(0,1) =invV(1,0) =invV(1,1) = 0;
	V(0,0) = V(0,1) =V(1,0) =V(1,1) = 0;

	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){
			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			si = sigma->Eval(phi,costheta);
			xi = sigmab->Eval(phi,costheta);

			eff = effhist->GetBinContent(i,j);
			Npositron = datahist->GetBinContent(i,j);

//			if(Ns*si*eff+Nb*xi*eff<=1e-8) continue;
			if(Ns*si*eff+Nb*fb[i-1][j-1]<=1e-8) continue;

//			invV(0,0) += tothist->GetBinContent(i,j)*(si*eff)**2/(Ns*si*eff+Nb*xi*eff)**2;
//			invV(0,1) += tothist->GetBinContent(i,j)*si*xi*eff**2/(Ns*si*eff+Nb*xi*eff)**2;
//			invV(1,0) = invV(0,1);
//			invV(1,1) += tothist->GetBinContent(i,j)*(xi*eff)**2/(Ns*si*eff+Nb*xi*eff)**2;

//			invV(0,0) += tothist->GetBinContent(i,j)*fs[i-1][j-1]**2/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;
//			invV(0,1) += tothist->GetBinContent(i,j)*fs[i-1][j-1]*fb[i-1][j-1]/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;
//			invV(1,0) = invV(0,1);
//			invV(1,1) += tothist->GetBinContent(i,j)*fb[i-1][j-1]**2/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;

			invV(0,0) += tothist->GetBinContent(i,j)*(si*eff)**2/(Ns*si*eff+Nb*fb[i-1][j-1])**2;
			invV(0,1) += tothist->GetBinContent(i,j)*(si*eff)*fb[i-1][j-1]/(Ns*si*eff+Nb*fb[i-1][j-1])**2;
			invV(1,0) = invV(0,1);
			invV(1,1) += tothist->GetBinContent(i,j)*fb[i-1][j-1]**2/(Ns*si*eff+Nb*fb[i-1][j-1])**2;

			cout<<"invV = "<<invV(0,0)<<"     "<<invV(0,1)<<"     "<<invV(1,0)<<"      "<<invV(1,1)<<endl;
		}
	}
	
	cout<<"Ns ====="<<Ns<<"     "<<"  Nb ======"<<Nb<<"            "<<tothist->Integral()<<endl;
	cout<<"V^-1 matrix = "<<invV(0,0)<<"     "<<invV(0,1)<<endl;
	cout<<"              "<<invV(1,0)<<"     "<<invV(1,1)<<endl;
//	V = invV.Invert();

	V(0,0) = 1/invV(0,0);
	V(0,1) = V(1,0) = 1/invV(0,1);
	V(1,1) = 1/invV(1,1);

	cout<<"V matrix    = "<<V(0,0)<<"      "<<V(0,1)<<endl;
	cout<<"              "<<V(1,0)<<"      "<<V(1,1)<<endl;
	cout<<" tothist bin = "<<tothist->Integral()<<"==="<<V(0,0)+2*V(0,1)+V(1,1)<<endl;
	
	//calculate sWeight in sFit technique

	double normweight=0.;

	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){

			phi = effhist->GetXaxis()->GetBinCenter(i);
			costheta = effhist->GetYaxis()->GetBinCenter(j);
			si = sigma->Eval(phi,costheta);
			xi = sigmab->Eval(phi,costheta);

			eff = effhist->GetBinContent(i,j);
			Npositron = datahist->GetBinContent(i,j);

//			invV(0,0) = tothist->GetBinContent(i,j)*(si*eff)**2/(Ns*si*eff+Nb*xi*eff)**2;
//			invV(1,1) = tothist->GetBinContent(i,j)*(xi*eff)**2/(Ns*si*eff+Nb*xi*eff)**2;
//			invV(0,1) = tothist->GetBinContent(i,j)*(si*eff*xi*eff)/(Ns*si*eff+Nb*xi*eff)**2;
//			invV(1,0) = invV(0,1);

//			V =invV.Invert();

			sWeight[i-1][j-1]=0.;

//			if(Ns*si*eff+Nb*xi*eff<=1e-8) continue;
			if(Ns*si*eff+Nb*fb[i-1][j-1]<=1e-8) continue;

//			sWeight[i-1][j-1] = (V(0,0)*si+V(0,1)*xi)*eff/(Ns*si*eff+Nb*xi*eff);
			sWeight[i-1][j-1] = (V(0,0)*si*eff+V(0,1)*fb[i-1][j-1])/(Ns*si*eff+Nb*fb[i-1][j-1]);
//			sWeight[i-1][j-1] = (V(0,0)*fs[i-1][j-1]+V(0,1)*fb[i-1][j-1])/(tothist->Integral());

			if(sWeight[i-1][j-1]!=sWeight[i-1][j-1]) continue;
			if(si*eff<1e-16) continue;
//			if(tothist->GetBinContent(i,j)==0x0) continue;

			cout<<"sWeight ==========="<<sWeight[i-1][j-1]<<endl;
//			normweight += tothist->GetBinContent(i,j)*sWeight[i-1][j-1];
			normweight += sWeight[i-1][j-1];
			sum += si*eff;
			sum1 += xi*eff;
			result += -1*sWeight[i-1][j-1]*(tothist->GetBinContent(i,j))*TMath::Log(si*eff); //sFit technique
		}
	}
	cout<<"sum ===="<<sum<<endl;
	cout<<"normweight ="<<normweight<<endl;
	return result;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

	f = likelihoodfcn(0,par);
	cout<<"f =================="<<f<<endl;
}

void matrix(int file,int trig,int pt,int frame){

	correcteddata(file,trig,pt,frame,1);

	TFile* savefile = new TFile(Form("~/polresults/20160707/sfit/sys_%d/functional_%d_%d_%d_%d.root",file,file,trig,pt,frame),"recreate");
	savefile->cd();
	datahist->Write();	

	int nbinxeff = effhist->GetNbinsX();
	int nbinyeff = effhist->GetNbinsY();

	double sumsquare;
	Ns = datahist->Integral();
	Nb = bkghist->Integral();

	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){
			fs[i-1][j-1] = datahist->GetBinContent(i,j)/datahist->Integral();
			fb[i-1][j-1] = bkghist->GetBinContent(i,j)/bkghist->Integral();
		}
	}

// minimization of likelihood function 
	TMinuit *gMinuit = new TMinuit(2);
	gMinuit->SetFCN(fcn);

Double_t arglist[10];
Int_t ierflg =0;
arglist[0]=1;
gMinuit->mnexcm("SET ERR", arglist,1,ierflg);

static Double_t vstart[4] = {-0.,-0.05,0.,0.05};
static Double_t step[4] = {0.01,0.01,0.01,0.01};
gMinuit->mnparm(0,"a1",vstart[0],step[0],-1,1,ierflg);
gMinuit->mnparm(1,"a2",vstart[1],step[1],-1,1,ierflg);
//gMinuit->mnparm(2,"a3",vstart[2],step[2],-1,1,ierflg);
//gMinuit->mnparm(3,"a4",vstart[3],step[3],-1,1,ierflg);

arglist[0] =500;
arglist[1]=1.;
//gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

Double_t amin,edm,errdef;
Int_t nvpar,nparx,icstat;
gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	
	/*
	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){
			
//			fs[i-1][j-1] = datahist->GetBinContent(i,j)/datahist->Integral()/(0.04*TMath::Pi());
//			fb[i-1][j-1] = bkghist->GetBinContent(i,j)/bkghist->Integral()/(0.04*TMath::Pi());

			fs[i-1][j-1] = datahist->GetBinContent(i,j)/datahist->Integral();
			fb[i-1][j-1] = bkghist->GetBinContent(i,j)/bkghist->Integral();

			cout<<"fs[i-1][j-1] = "<<fs[i-1][j-1]<<"     "<<fb[i-1][j-1]<<endl;

			cout<<"      "<< Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1]<<"      "<<sumsquare<<endl;
			sumsquare = Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1]; 
			sumsquare = sumsquare*sumsquare;
			cout<<"sumsquare = "<<sumsquare<<"==========="<<tothist->GetBinContent(i,j)**2<<endl;
		
			if(fs[i-1][j-1]==0 && fb[i-1][j-1]==0) continue;

//			invV(0,0) += tothist->GetBinContent(i,j)*fs[i-1][j-1]**2/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;
//			invV(0,1) += tothist->GetBinContent(i,j)*fs[i-1][j-1]*fb[i-1][j-1]/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;
//			invV(1,0) = invV(0,1) ;
//			invV(1,1) += tothist->GetBinContent(i,j)*fb[i-1][j-1]**2/(Ns*fs[i-1][j-1]+Nb*fb[i-1][j-1])**2;

			invV(0,0) = tothist->GetBinContent(i,j)*fs[i-1][j-1]**2/sumsquare;
			invV(0,1) = tothist->GetBinContent(i,j)*fs[i-1][j-1]*fb[i-1][j-1]/sumsquare;
			invV(1,0) = invV(0,1);
			invV(1,1) = tothist->GetBinContent(i,j)*fb[i-1][j-1]**2/sumsquare;

//			invV(0,0) = tothist->GetBinContent(i,j)*(datahist->GetBinContent(i,j)/tothist->GetBinContent(i,j))**2/sumsquare;
//			invV(0,1) = datahist->GetBinContent(i,j)*bkghist->GetBinContent(i,j)/tothist->GetBinContent(i,j)/sumsquare;
//			invV(1,0) = invV(0,1);
//			invV(1,1) = tothist->GetBinContent(i,j)*(bkghist->GetBinContent(i,j)/tothist->GetBinContent(i,j))**2/sumsquare;

//			V =invV;

//			V = invV.Invert();	

			double det = fabs(invV(0,0)*invV(1,1) - invV(0,1)**2);
			V(0,0) = invV(1,1)/det;
			V(0,1) = -1*invV(1,0)/det;
			V(1,0) = V(0,1);
			V(1,1) = invV(0,0)/det;

			cout<<"after invert V =========="<<V(0,0)+2*V(0,1)+V(1,1)<<" ============="<<tothist->GetBinContent(i,j)<<endl;
			cout<<"tothist bin content = "<<tothist->GetBinContent(i,j)<<"======="<<fs[i-1][j-1]<<"    "<<fb[i-1][j-1]<<endl;
	
			cout<<"inV matrix = "<<invV(0,0)<<"    "<<invV(0,1)<<"     "<<invV(1,1)<<endl;
			cout<<"matrix = "<<V(0,0)<<"    "<<V(0,1)<<"     "<<V(1,1)<<endl;

			matrix[i-1][j-1][0][0] = V(0,0);
			matrix[i-1][j-1][0][1] = V(0,1);
			matrix[i-1][j-1][1][0] = V(1,0);
			matrix[i-1][j-1][1][1] = V(1,1);
		}
	}
*/
	
/*
	for(int i=1;i<=nbinxeff;i++){
		for(int j=1;j<=nbinyeff;j++){
			double phi = effhist->GetXaxis()->GetBinCenter(i);
			double costheta = effhist->GetYaxis()->GetBinCenter(j);

	//		fs[i-1][j-1] = sfcn->Eval(phi,costheta)/sfcn->Integral(-TMath::Pi(),TMath::Pi(),-1,1);
	//		fb[i-1][j-1] = bfcn->Eval(phi,costheta)/bfcn->Integral(-TMath::Pi(),TMath::Pi(),-1,1);

//			fs[i-1][j-1] = datahist->GetBinContent(i,j)/datahist->Integral();
//			fb[i-1][j-1] = bkghist->GetBinContent(i,j)/bkghist->Integral();

//			Ns = datahist->Integral();
//			Nb = bkghist->Integral();

				}
	}
*/






/*
	TF2* likelihood;

	double ltheta,lphi,lthetaerr,lphierr;

	for(int iteration = 0;iteration<2;iteration++){
		if(iteration==0)likelihood = new TF2("likelihood",likelihoodfcn,-1,1,-1,1,0);
		else likelihood = new TF2("likelihood",likelihoodfcn,ltheta-0.5,ltheta+0.5,lphi-0.5,lphi+0.5,0);

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
		//		cout<<"normalization = "<<normalization<<endl;
		sigma->SetParameter(0,par0ini/normalization);
		//		cout<<"parameter 0 = "<<sigma->GetParameter(0)<<endl;

		//		cout<<"integral ="<<sigma->Integral(-TMath::Pi(),TMath::Pi(),-1,1)<<endl;	
		sigma->Write("",TObject::kOverwrite);

		double deltacontour = 0.5;
		double contour = likelihood->Eval(ltheta,lphi) + deltacontour;	

		likelihood->SetContour(1,&contour);
		likelihood->Draw("CONT Z LIST");
		c->Update();

		TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

		TList* contlevel = NULL;
		TGraph* gc2 = NULL;
		contlevel = (TList*)conts->At(0);	
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


		TGraphErrors* lambda = new TGraphErrors(1,&ltheta,&lphi,&lthetaerr,&lphierr);
		lambda->SetName("lambda");
		lambda->SetTitle("lambda parameter;#lambda_{#theta};#lambda_{#phi}");

		//	c->SaveAs("functional.pdf");
		//	effhist->Write();
		//	sigmafcn->Write();
		lambda->Write();
	}
*/

}

void rawMCdata(int file,int trig,int pt,int frame){

	gRandom = new TRandom3();
	gRandom->SetSeed(0);

	int rebin=1;
	rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),0),"read");
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

	int max = (rawdata3D[file][trig][pt][frame][2]->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
	int min = (rawdata3D[file][trig][pt][frame][2]->GetZaxis())->FindBin(PtEdge[pt]+0.001);	

	rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);

	rawdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
	//					rawdata2D[file][trig][pt][frame]->SetName(Form("rawdata_%d_%s_pT%d_frame%d",file,trigName[trig].Data(),pt,frame));
	rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
	rawdata2D[file][trig][pt][frame]->Sumw2();
	if(rebin==1){
		rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
		rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
	}

	double sample = rawdata2D[file][trig][pt][frame]->Integral();
	cout<<"sample = "<<sample<<endl;

	efficiencyfile = new TFile("rootfile/OutFile_sys0.root","read"); 
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

	eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
	eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
	eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
	eff3D[file][trig][pt][frame][0]->Sumw2();
	eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));
	eff3D[file][trig][pt][frame][1]->Sumw2();

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

	rawdatafile1[file][trig] = new TFile("../zbackupfiles/histograms20161025.root","read");//marker
	sigmafcn = (TF2*)rawdatafile1[file][trig]->Get("sigma_50_50");
	sigmafcn->SetName(Form("sigma_%d_%d_%d_%d",file,trig,pt,frame));

	double phi,costheta;			

	//	TFile* rawmcfile = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/raw_%d_%d_%d_%d.root",file,file,trig,pt,frame),"recreate");
	//	rawmcfile->cd();

	rawdata2D[file][trig][pt][frame] = new TH2F(Form("raw_%d_%d_%d_%d",file,trig,pt,frame),Form("raw_%d_%d_%d_%d",file,trig,pt,frame),40,-TMath::Pi(),TMath::Pi(),40,-1,1);
	rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
	rawdata2D[file][trig][pt][frame]->Sumw2();

	for(int i=0;i<sample;){	
		sigmafcn->GetRandom2(phi,costheta);	
		int xbin = eff2D[file][trig][pt][frame][2]->GetXaxis()->FindBin(phi);
		int ybin = eff2D[file][trig][pt][frame][2]->GetYaxis()->FindBin(costheta);
		if(gRandom->Uniform()<eff2D[file][trig][pt][frame][2]->GetBinContent(xbin,ybin)){
			rawdata2D[file][trig][pt][frame]->Fill(phi,costheta);		
			i++;
		}
	}

	//	rawdata2D[file][trig][pt][frame]->Scale(sample/10000);
	//					rawdata2D[file][trig][pt][frame]->Scale(sample);

	if(rebin==1){
		rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
		rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
	}
	//	rawdata2D[file][trig][pt][frame]->Write();	
	datahist = (TH2F*)rawdata2D[file][trig][pt][frame]->Clone("datahist");
	effhist = (TH2F*)eff2D[file][trig][pt][frame][2]->Clone("effhist");
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
	}
}
