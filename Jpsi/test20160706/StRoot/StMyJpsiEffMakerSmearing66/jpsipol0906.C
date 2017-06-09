#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TString.h"
//#include "StRoot/StMyElectronMaker/StMyElectron.h"
//#include "StRoot/StMyElectronMaker/StMyElectronEvent.h"
//#include "StRoot/StMyElectronMaker/StMyElectronMaker.h"
#include "/afs/rhic.bnl.gov/star/packages/SL12d_embed/StRoot/StarClassLibrary/SystemOfUnits.h"
#include <iostream>
//#define Pi 3.14159265
//#define twoPi 6.28318530
#define EMASS 0.000511
#define nHitsFitCut 20

#include "cuts.h"


Double_t InvPt(Double_t *x,Double_t *par)
{
	return par[0]*TMath::Power(TMath::Exp(par[1]*x[0])+x[0]/par[2],par[3]);
}

Double_t InvPt_test(Double_t *x,Double_t *par)
{
	return par[0]*TMath::Power(1+TMath::Power(x[0]*par[1],2),-par[2]);
}
Double_t sigmaFit(Double_t *x,Double_t *par)
{
	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];

}
class StMyElectronEvent;
StMyElectronEvent *mElectronEvent;
class StMyElectron;
StMyElectron *mElectron;
StMyElectron *mElectron2;



TH1F *hMcVertexZ = new TH1F("mcVertexZ","mcVertexZ;Vz^{mc} (cm)",400,-200,200);
TH2F *hMcVertexXY = new TH2F("mcVertexXY","mcVertexXY;Vx^{mc} (cm);Vy^{mc} (cm)", 40, -2, 2, 40, -2, 2);
TH1F *hRcVertexZ = new TH1F("rcVertexZ","rcVertexZ;Vz^{rc} (cm)",400,-200,200);
TH1F *hVertexZdiff = new TH1F("vertexZdiff", "vertexZdiff;(Vz^{rc}-Vz^{mc} (cm)", 100, -5, 5);
TH1F *hRefMult = new TH1F ("refMult","refMult;Reference Multiplicity",1000,0,1000);
TH1F *hRefMultCut = new TH1F ("refMultCut","refMultCut;Reference Multiplicity after Cut",1000,0,1000);
TH2F *hMCdeltaEtavsdeltaPhi = new TH2F("hMCdeltaEtavsdeltaPhi", "D_Eta vs D_Phi; #Delta #eta; #Delta #phi", 800, -4, 4, 360, -TMath::Pi(), TMath::Pi());
TH2F *hRCdeltaEtavsdeltaPhi = new TH2F("hRCdeltaEtavsdeltaPhi", "D_Eta vs D_Phi; #Delta #eta; #Delta #phi", 800, -4, 4, 360, -TMath::Pi(), TMath::Pi());

TH1F *hNJpsi = new TH1F("hNJpsi","#Jpsi;#Jpsi;#Event", 10, 0, 10);
TH1D *hMcJpsiPt   = new TH1D("mcJpsiPt",  "mcJpsiPt;p_{T}^{mc} (GeV/c)",  300, 0, 30.);
TH1D *hMcJpsiPt_Or   = new TH1D("mcJpsiPt_Or",  "mcJpsiPt;p_{T}^{mc} (GeV/c)",  300, 0, 30.);

TH1F *hMcJpsiY  = new TH1F("mcJpsiY", "mcJpsiY; Y^{mc};", 100, -5,  5.);
TH2D *hMcJpsiPtY = new  TH2D("mcJpsiPtY","mcJpsiPtY;Y^{mc};p_{T}^{mc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hMcJpsiPtY_Or = new  TH2D("mcJpsiPtY_Or","mcJpsiPtY;Y^{mc};p_{T}^{mc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH1F *hMcJpsiPhi  = new TH1F("mcJpsiPhi", "mcJpsiPhi,#phi^{mc}", 360, -TMath::Pi(), TMath::Pi());
TH2D *hMcJpsiPtM = new  TH2D("mcJpsiPtM","mcJpsiPtM;M_{inv}^{mc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{mc} (GeV/c)", 300, 0, 30, 1000,0,5);
TH2D *hMcJpsiPtM_1 = new  TH2D("mcJpsiPtM_1","mcJpsiPtM;M_{inv}^{mc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{mc} (GeV/c)", 300, 0, 30, 1000,0,5);

TH1D *hMCElectronPt = new TH1D("mcElectronPt","input electron pt",300,0,30);

TH2D *hmcPtvsrcPt = new TH2D("mcPtvsrcPt","mcPt vs rcPt;p_{T}^{rc};p_{T}^{mc}",3500,0,35,4000,0,40);
TH2D *hmcPtvsrcPt_Cut = new TH2D("mcPtvsrcPt_Cut","mcPt vs rcPt;p_{T}^{rc};p_{T}^{mc}",3500,0,35,4000,0,40);


TH1D *hRcJpsiPt   = new TH1D("rcJpsiPt",  "rcJpsiPt;p_{T}^{rc} (GeV/c)",  300, 0, 30.);

TH1D *hRcJpsiPt_Or   = new TH1D("rcJpsiPt_Or",  "rcJpsiPt;p_{T}^{rc} (GeV/c)",  300, 0, 30.);
TH1F *hRcJpsiY  = new TH1F("rcJpsiY", "rcJpsiY;y^{rc};", 100, -5,  5.);
TH2D *hRcJpsiPtY = new TH2D("rcJpsiPtY", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hRcJpsiPtY_Or = new TH2D("rcJpsiPtY_Or", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hRcJpsiPtY_B = new TH2D("rcJpsiPtY_B", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH2D *hRcJpsiPtY_BOr = new TH2D("rcJpsiPtY_BOr", "rcJpsiPtY;y^{rc};p_{T}^{rc} (GeV/c)", 40, -2, 2, 300, 0, 30);
TH1F *hRcJpsiPhi  = new TH1F("rcJpsiPhi", "rcJpsiPhi;#phi^{rc}", 360, -TMath::Pi(), TMath::Pi());
TH2D *hRcJpsiPtM = new TH2D("rcJpsiPtM", "rcJpsiPtM;M_{inv}^{rc}(e^{+}e^{-}) [GeV/c^{2}];p_{T}^{rc} (GeV/c)", 300, 0, 30, 1000,0,5);
TH2D *hRcJpsiPtDiff = new TH2D("rcJpsiPtDiff", "rcJpsiPtDiff;p_{T}^{mc}(GeV/c);p_{T}^{rc}/p_{T}^{mc}", 300, 0, 30,1000, 0, 2);
TH2D *hRcJpsiPtDiff_rc = new TH2D("rcJpsiPtDiff_rc", "rcJpsiPtDiff;p_{T}^{rc}(GeV/c);p_{T}^{rc}/p_{T}^{mc}", 300, 0, 30,1000, 0, 2);

TH3F *hRcJpsiYDiff = new TH3F("rcJpsiYDiff", "rcJpsiYDiff;p_{T}^{mc}(GeV/c);y^{mc};y^{rc}-y^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);
TH3F *hRcJpsiPhiDiff = new TH3F("rcJpsiPhiDiff", "rcJpsiPhiDiff;p_{T}^{mc}(GeV/c);y^{mc};#phi^{rc}-#phi^{mc}", 250,0,25,40, -1, 1, 1000, -0.1, 0.1);

TH2D *hJpsiMc = new TH2D("hJpsiMc","hJpsiMc;y^{mc};p_{T}^{mc} (GeV/c)",80,-2,2,300,0,30);

TH2D *hHt0JpsiTrg = new TH2D("hHt0JpsiTrg","hHt0JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiAdc0 = new TH2D("hHt0JpsiAdc0","hHt0JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiPE = new TH2D("hHt0JpsiPE","hHt0JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiNSMD = new TH2D("hHt0JpsiNSMD","hHt0JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt0JpsiDist = new TH2D("hHt0JpsiDist","hHt0JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);


TH2D *hHt2JpsiTrg = new TH2D("hHt2JpsiTrg","hHt2JpsiTrg;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiAdc0 = new TH2D("hHt2JpsiAdc0","hHt2JpsiAdc0;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiPE = new TH2D("hHt2JpsiPE","hHt2JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiPE_Or = new TH2D("hHt2JpsiPE_Or","hHt2JpsiPE;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH2D *hHt2JpsiNSMD = new TH2D("hHt2JpsiNSMD","hHt2JpsiNSMD;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);
TH2D *hHt2JpsiDist = new TH2D("hHt2JpsiDist","hHt2JpsiDist;y^{rc};p_{T}^{rc} (GeV/c)",80,-2,2,300,0,30);

TH3F *hHt0JpsiMassPE = new TH3F("hHt0JpsiMassPE","hHt0JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt0JpsiMassDist = new TH3F("hHt0JpsiMassDist","hHt0JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassPE_1 = new TH3F("hHt2JpsiMassPE_1","hHt2JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassPE = new TH3F("hHt2JpsiMassPE","hHt2JpsiMassPE;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);
TH3F *hHt2JpsiMassDist = new TH3F("hHt2JpsiMassDist","hHt2JpsiMassDist;y^{rc};p_{T}^{rc} (GeV/c);Mass (Gev/c^{2})",80,-2,2,300,0,30,400,0,4);

TH3F *hJpsi3DMc = new TH3F("hJpsi3DMc","hJpsi3DMc;y^{mc};Vz(cm);p_{T}^{mc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt0Jpsi3DAdc0 = new TH3F("hHt0Jpsi3DAdc0","hHt0Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DPE = new TH3F("hHt0Jpsi3DPE","hHt0Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DNSMD = new TH3F("hHt0Jpsi3DNSMD","hHt0Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt0Jpsi3DDist = new TH3F("hHt0Jpsi3DDist","hHt0Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH3F *hHt2Jpsi3DAdc0 = new TH3F("hHt2Jpsi3DAdc0","hHt2Jpsi3DAdc0;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DPE = new TH3F("hHt2Jpsi3DPE","hHt2Jpsi3DPE;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DNSMD = new TH3F("hHt2Jpsi3DNSMD","hHt2Jpsi3DNSMD;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);
TH3F *hHt2Jpsi3DDist = new TH3F("hHt2Jpsi3DDist","hHt2Jpsi3DDist;y^{rc};Vz(cm);p_{T}^{rc} (GeV/c)",80,-2,2,200,-200,200,300,0,30);

TH2D *hCommonhitsvsRCPt = new TH2D("hCommonhitsvsRCPt","commonhits vs RC pT;tpc commonHits;RC p_{T} (GeV/c)",50,0,50,300,0,300);
TH2D *hCommonhitsvsMCPt = new TH2D("hCommonhitsvsMCPt","commonhits vs MC pT;tpc commonHits;MC p_{T} (GeV/c)",50,0,50,300,0,300);

TH2D *hdeltaEtavsMCPt = new TH2D("hdeltaEtavsMCPt","delta Eta vs MCpT; #Delta #eta; MC p_{T} GeV/c", 800,-20,20,300,0,30);
TH2D *hdeltaEtavsRCPt = new TH2D("hdeltaEtavsRCPt","delta Eta vs RCpT; #Delta #eta; RC p_{T} GeV/c", 800,-20,20,300,0,30);
TH2D *hdeltaPhivsMCPt = new TH2D("hdeltaPhivsMCPt","delta phi vs MCpT; #Delta #phi; MC p_{T} GeV/c", 360, -TMath::Pi(), TMath::Pi(), 300, 0, 30);
TH2D *hdeltaPhivsRCPt = new TH2D("hdeltaPhivsRCPt","delta phi vs MCpT; #Delta #phi; RC p_{T} GeV/c", 360, -TMath::Pi(), TMath::Pi(), 300, 0, 30);
TH2D *hdeltaRvsRCPt = new TH2D("hdeltaRvsRCPt","delta R vs RCpT; #Delta R; RC p_{T} GeV/c", 800, -4, 4, 300, 0, 30);
TH2D *hdeltaRvsMCPt = new TH2D("hdeltaRvsMCPt","delta R vs MCpT; #Delta R; MC p_{T} GeV/c", 800, -4, 4, 300, 0, 30);

TH2F *hHt2Adc0vsPt = new TH2F("hHt2Adc0vsPt","Ht2Adc0vsPt;p_{T}^{mc} (GeV/c);Adc0 (HT2)",300,0,30,1000,0,1000);
TH2F *hHt2Adc0vsrcPt = new TH2F("hHt2Adc0vsrcPt","Ht2Adc0vsrcPt;p_{T}^{rc} (GeV/c);Adc0 (HT2)",300,0,30,1000,0,1000);
TH2F *hHt2Adc0vsPt_1 = new TH2F("hHt2Adc0vsPt_1","Ht2Adc0vsPt_1;p_{T}^{mc} (GeV/c);Adc0 (HT2)",250,0,25,1000,0,1000);
TH3F *hHt2Adc0vsPtY = new TH3F("hHt2Adc0vsPtY","Ht2Adc0vsPtY;p_{T}^{mc} (GeV/c);#phi;Adc0 (HT0)",250,0,25,40,-1,1,1000,0,1000);

TH1F *hMcJpsiPol = new TH1F("hMcJpsiPol","Theta",300,0,TMath::Pi());
TH2F *hMcJpsiPolThetaPt = new TH2F("hMcJpsiPolThetaPt","Pt vs Cos(#theta)",300,-1,1,300,0,10);
TH1F *hHtJpsiPol = new TH1F("hHtJpsiPol","Theta",300,0,TMath::Pi());
TH2F *hHtJpsiPolThetaPt = new TH2F("hHtJpsiPolThetaPt","Pt vs Cos(#theta)",300,-1,1,300,0,10);

TH1F *hNreal = new TH1F("hNrea","nreal",300,0,100);

void jpsipol0906(const char* fileInDir="./out_mb__100_20151202/", const char* OutFile="jpsiEff_200", const int mCentCut1=0, const int mCentCut2=9){
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StMyElectronMaker");
	const Int_t nFilesMax = 10000;

	TF1 *myGaus = new TF1("myGaus","gaus",-10,10);
	myGaus->SetParameters(1,0,0.9);
	TF1 *myGaus_1 = new TF1("myGaus_1","gaus",-10,10);
	myGaus_1->SetParameters(1,0,0.9);

	mElectronEvent = new StMyElectronEvent();
	TChain chMc("mcT");
	chMc.SetBranchAddress("mcE",&mElectronEvent);

	const Double_t mDsmAdcCut[4] = {11, 15, 18, 25};
	const float pi = 3.1415926;
	void* dir = gSystem->OpenDirectory(gSystem->ExpandPathName(fileInDir));
	int nruns=0;
	char *file_name;
	TString Tname;
	char file_list[nFilesMax][256];
	do {
		file_name = (char*)gSystem->GetDirEntry(dir);
		Tname=file_name;
		if(file_name && Tname.Contains("myminimc.root")&&Tname.Contains("emb")) {
			sprintf(file_list[nruns],"%s/%s",fileInDir,file_name);
			chMc.Add(file_list[nruns]);
			cout << " read in " << file_list[nruns] << endl;
			nruns++;
		}
	} while (file_name && nruns<=nFilesMax);

	int iret = 0;
	int nb=0;
	cout << chMc.GetEntries() << " events in chain" << endl;
	int nevents = chMc.GetEntries();
	nevents = 100;
	//	cout<<"nevents ========"<< nevents <<endl;
	//nevents = 6000;
	//+--------------------------------+
	//| initialize wegiht function     |
	//+--------------------------------+
	TF1 *function_W = new TF1("function_W",InvPt,0,30,4);
	// function_W->SetParameters(0.4,-0.4042,4.053,-6);
	function_W->SetParameters(0.4,-0.4796,4.229,-7.54);
	//+---------------------------------+
	//|   loop the events               |
	//+---------------------------------+
	TF1 *function_sigma = new TF1("function_sigma",sigmaFit,0,30,3);
	function_sigma->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);


	for (int i = 0; i < nevents; i++) {
		nb +=chMc.GetEvent(i);
		chMc.GetEvent(i);
		if(i%10==0) cout<<"event# "<<i<<endl;
		int n;
		if (!mElectronEvent) continue;
		if (mElectronEvent->eventID()<=0) continue;

		hMcVertexZ->Fill(mElectronEvent->mcVertexZ());
		hMcVertexXY->Fill(mElectronEvent->mcVertexX(), mElectronEvent->mcVertexY());
		hRcVertexZ->Fill(mElectronEvent->vertexZ());
		hVertexZdiff->Fill(mElectronEvent->vertexZ() - mElectronEvent->mcVertexZ());

		Double_t vz = mElectronEvent->vertexZ();
		int Mult = mElectronEvent->refMultPos() + mElectronEvent->refMultNeg();
		if(vz<mVzCut[0]||vz>mVzCut[1]) continue;
		hRefMult->Fill(Mult);
		//cout<<"vertexZ = "<<mElectronEvent->vertexZ()<<endl;
		int cent = getCentrality(Mult);
		if(cent<mCentCut1||cent>mCentCut2) continue;
		hRefMultCut->Fill(Mult);

		TLorentzVector jpsiMcRaw(0.,0.,0.,0.),jpsiMc(0.,0.,0.,0.), ePosMc(0.,0.,0.,0.), eNegMc(0.,0.,0.,0.);
		TLorentzVector jpsiRcRaw(0.,0.,0.,0.),jpsiRc(0.,0.,0.,0.), ePosRc(0.,0.,0.,0.), eNegRc(0.,0.,0.,0.);
		TLorentzVector jpsiMc_tem(0.,0.,0.,0.), ePosMc_tem(0.,0.,0.,0.), eNegMc_tem(0.,0.,0.,0.);

		Int_t nJpsi = 0;
//		cout<<"nReal ========="<< mElectronEvent->nReal()<<endl;
		hNreal->Fill(mElectronEvent->nReal());
		for(int j=0;j<mElectronEvent->nReal();j++){
			mElectron = (StMyElectron) mElectronEvent->real()->UncheckedAt(j);
			//cout<<j<<"  "<<mElectron->geantId<<"  "<<mElectron->pGeantId<<endl;
			if(mElectron->pGeantId!=160) continue;//not from Jpsi electron
			//if(mElectron->pGeantId>0) continue;//not original electron
			if(mElectron->mcId<0) continue;
			hCommonhitsvsMCPt->Fill(mElectron->tpcCommonHits,mElectron->mcPt);
			hCommonhitsvsRCPt->Fill(mElectron->tpcCommonHits,mElectron->pt);
			bool tag = kFALSE;
			for(int k=0;k<mElectronEvent->nReal();k++) {
				mElectron2 = (StMyElectron) mElectronEvent->real()->UncheckedAt(k);
				//cout<<k<<"  "<<mElectron2->geantId<<"  "<<mElectron2->pGeantId<<endl;
				if(mElectron2->pGeantId!=160) continue;
				if(mElectron2->mcId<0) continue;
				if(mElectron2->mcId==mElectron->mcId) continue;
				//cout<<deta<<"  "<<dphi<<"  "<<tag<<endl;
				if(mElectron->geantId==2&&mElectron2->geantId==3) {
					ePosMc_tem.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					eNegMc_tem.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
				} else if(mElectron->geantId==3&&mElectron2->geantId==2) {
					eNegMc_tem.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					ePosMc_tem.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
				} else {
					//	  cout<<"1!!!! something wrong with the J/psi decay !!!!!!"<<endl;
					continue;
				};
				jpsiMc_tem = ePosMc_tem + eNegMc_tem;
				Double_t deta = mElectron->mcY - mElectron2->mcY;
				Double_t dphi = mElectron->mcPhi - mElectron2->mcPhi;
				while(dphi>2*TMath::Pi()) dphi -= 2.*TMath::Pi();
				while(dphi<0) dphi += 2.*TMath::Pi();
				while(dphi>TMath::Pi()) dphi = dphi - 2*TMath::Pi();
				Double_t dReta = mElectron->eta - mElectron2->eta;
				Double_t dRphi = mElectron->phi - mElectron2->phi;
				while(dRphi>2*TMath::Pi()) dRphi -= 2.*TMath::Pi();
				while(dRphi<0) dRphi += 2.*TMath::Pi();
				while(dRphi>TMath::Pi()) dRphi = dRphi - 2*TMath::Pi();
				if(mElectron->pId!=mElectron2->pId)
				{
					if(jpsiMc_tem.Rapidity()<1&&jpsiMc_tem.Rapidity()>-1)
					{
						hMCdeltaEtavsdeltaPhi->Fill(deta, dphi);	
						hRCdeltaEtavsdeltaPhi->Fill(dReta, dRphi);	
					}
				}
				if(TMath::Abs(deta)<0.1&&TMath::Abs(dphi)<0.5&&mElectron2->pId!=mElectron->pId) tag = kTRUE;
			}
			//cout<<tag<<endl;
			if(tag) continue;
			Double_t pt_tem = mElectron->mcPt;
			hMCElectronPt->Fill(pt_tem);
			Int_t count1 =0;
			for(int k=j+1; k<mElectronEvent->nReal();k++) {
				mElectron2 =  (StMyElectron) mElectronEvent->real()->UncheckedAt(k);
				//cout<<"xxx"<<k<<"  "<<mElectron2->geantId<<"  "<<mElectron2->pGeantId<<endl;
				if(mElectron2->pGeantId!=160) continue;
				if(mElectron2->mcId<0) continue;
				if(mElectron2->mcId==mElectron->mcId) continue;
				if(mElectron2->pId!=mElectron->pId) continue;//from same parent
				if(mElectron->geantId==2&&mElectron2->geantId==3) {
					ePosMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					eNegMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
//					ePosRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
//					eNegRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
				} else if(mElectron->geantId==3&&mElectron2->geantId==2) {
					eNegMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					ePosMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
//					eNegRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
//					ePosRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
				} else {
					cout<<"!!!! something wrong with the J/psi decay !!!!!!"<<endl;
					continue;
				}
				nJpsi++;
				jpsiMc = ePosMc + eNegMc;
				//	cout<<jpsiMc.M()<<endl;

				//  Double_t weight1 = function_W->Eval(jpsiMc.Pt());
				Double_t weight1 = 1;
				
				if(mElectron->mcId>=0 && mElectron2->mcId>=0) {
					hMcJpsiPt_Or->Fill(jpsiMc.Pt());
					hMcJpsiPt->Fill(jpsiMc.Pt(),weight1);
					hMcJpsiY->Fill(jpsiMc.Rapidity());
					hMcJpsiPtY->Fill(jpsiMc.Rapidity(), jpsiMc.Pt(),weight1);
					hMcJpsiPtY_Or->Fill(jpsiMc.Rapidity(), jpsiMc.Pt());
					hMcJpsiPhi->Fill(jpsiMc.Phi());
					hMcJpsiPtM->Fill(jpsiMc.Pt(), jpsiMc.M());
					hMcJpsiPtM_1->Fill(jpsiMc.Pt(), jpsiMc.M(),weight1);
					
					TLorentzVector ePosMcRest = ePosMc;
					ePosMcRest.Boost(-jpsiMc.Px()/jpsiMc.E(),-jpsiMc.Py()/jpsiMc.E(),-jpsiMc.Pz()/jpsiMc.E());
					Double_t dtheta = jpsiMc.Angle(ePosMcRest.Vect());
					Double_t costheta = TMath::Cos(dtheta);
					hMcJpsiPol->Fill(dtheta);
					hMcJpsiPolThetaPt->Fill(costheta,jpsiMc.Pt());
				
					if(mElectron->id>=0)
					{
						if(fabs(mElectron->eta)<1)
						{
							if(mElectron->dsmAdc0>18)
							{
								hHt2Adc0vsPt->Fill(mElectron->mcPt, mElectron->adc0,weight1);
								hHt2Adc0vsrcPt->Fill(mElectron->pt, mElectron->adc0,weight1);
								hHt2Adc0vsPtY->Fill(mElectron->mcPt, mElectron->mcEta, mElectron->adc0,weight1);
								if(mElectron->nFitPts>=nHitsFitCut&&mElectron->dca<1)
								{
									hHt2Adc0vsPt_1->Fill(mElectron->pt,mElectron->adc0,weight1);
								}
							}
						}
					}
					float deltaeta = mElectron->mcEta - mElectron2->mcEta;
					float deltaphi = mElectron->mcPhi - mElectron2->mcPhi;
					while(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
					while(deltaphi<0) deltaphi += 2.*TMath::Pi();
					while(deltaphi>TMath::Pi()) deltaphi = deltaphi - 2*TMath::Pi();
					double deltaR =0;
					if(jpsiMc.Rapidity()<1&&jpsiMc.Rapidity()>-1.)
					{
						hdeltaEtavsMCPt->Fill(deltaeta, jpsiMc.Pt());
						hdeltaPhivsMCPt->Fill(deltaphi, jpsiMc.Pt());
						deltaR = TMath::Sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
						hdeltaRvsMCPt->Fill(deltaR, jpsiMc.Pt());
					}
				} else {
					cout<<" !!! something wrong with MC tracks !!!!!!"<<endl;
					continue;
				}//all of the MC J/psi below are good!!!!!!!!
				

				//////////////////////////////////////////
				Double_t mcPt_sigma = function_sigma->Eval(jpsiMc.Pt());
				//////////////////////////////////////////

			/*	if(mElectron->id>=0 && mElectron2->id>=0 &&mElectron->pt>0 && mElectron2->pt>0 ){ 
					jpsiRc = ePosRc + eNegRc;
					TLorentzVector ePosRest = ePosRc;
  				 	ePosRest.Boost(-jpsiRc.Px()/jpsiRc.E(),-jpsiRc.Py()/jpsiRc.E(),-jpsiRc.Pz()/jpsiRc.E());
					Double_t dtheta = jpsiRc.Angle(ePosRest.Vect());
					Double_t costheta = TMath::Cos(dtheta);			
					hHtJpsiPol->Fill(dtheta);
					hHtJpsiPolThetaPt->Fill(costheta,jpsiRc.Pt());
				}
			*/
				if(mElectron->id>=0&&mElectron2->id>=0)
				{
					cout<<"id========="<<mElectron->geantId<<"and"<<mElectron2->geantId<<endl;
					hRcJpsiPtY_B->Fill(jpsiRc.Rapidity(), jpsiRc.Pt(),weight1);
					hRcJpsiPtY_BOr->Fill(jpsiRc.Rapidity(), jpsiRc.Pt());
					bool Quailtyflag[2] = {kFALSE, kFALSE};
					double eta1 = mElectron->eta;
					double dca1 = mElectron->dca;
					double nHitsFit1 = mElectron->nFitPts;
					double eta2 = mElectron2->eta;
					double dca2 = mElectron2->dca;
					double nHitsFit2 = mElectron2->nFitPts;
					if((eta1<1&&eta1>-1) && (dca1<3) && (nHitsFit1<25))
					{
						Quailtyflag[0] = kTRUE;
					}
					if((eta2<1&&eta2>-1) && (dca2<3) && (nHitsFit1<25))
					{
						Quailtyflag[1] = kTRUE;
					}
					if(Quailtyflag[0]||Quailtyflag[1])hmcPtvsrcPt->Fill(jpsiRc.Pt(),jpsiMc.Pt());
					if(jpsiRc.Pt()>(jpsiMc.Pt()*(1.+3*mcPt_sigma)))continue;
					hmcPtvsrcPt_Cut->Fill(jpsiRc.Pt(),jpsiMc.Pt());
				}

				if(mElectron->id>=0&&mElectron2->id>=0){
					float deltaeta = mElectron->eta - mElectron2->eta;
					float deltaphi = mElectron->phi - mElectron2->phi;
					while(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
					while(deltaphi<0) deltaphi += 2.*TMath::Pi();
					while(deltaphi>TMath::Pi()) deltaphi = deltaphi - 2*TMath::Pi();
					double deltaR;
					if(jpsiRc.Rapidity()<1&&jpsiRc.Rapidity()>-1)
					{
						hdeltaEtavsRCPt->Fill(deltaeta, jpsiRc.Pt());
						hdeltaPhivsRCPt->Fill(deltaphi, jpsiRc.Pt());
						deltaR = TMath::Sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
						hdeltaRvsRCPt->Fill(deltaR, jpsiRc.Pt());
					}
					//cout<<mElectron->id<<endl;
					Double_t rcPt = jpsiRc.Pt();
					Double_t rcY = jpsiRc.Rapidity();
					Double_t rcPhi = jpsiRc.Phi();
					hRcJpsiPt_Or->Fill(rcPt);
					// hweightPt_1->Fill(rcPt,weight1);
					hRcJpsiPt->Fill(rcPt,weight1);

					hRcJpsiY->Fill(rcY);
					hRcJpsiPtY->Fill(rcY, rcPt,weight1);
					hRcJpsiPtY_Or->Fill(rcY, rcPt);
					hRcJpsiPhi->Fill(rcPhi,weight1);
					hRcJpsiPtM->Fill(rcPt, jpsiRc.M(),weight1);

					hRcJpsiPtDiff_rc->Fill(rcPt,rcPt/jpsiMc.Pt());

					if(fabs(rcY)<1){
						Double_t mcPt = jpsiMc.Pt();
						Double_t mcY = jpsiMc.Rapidity();
						Double_t mcPhi = jpsiMc.Phi();
						hRcJpsiPtDiff->Fill(mcPt, rcPt/mcPt);

						hRcJpsiYDiff->Fill(mcPt, mcY, rcY-mcY);
						hRcJpsiPhiDiff->Fill(mcPt, mcY, rcPhi-mcPhi);

					}//end of |eta|<1 cut
				}//end of RC electron loop
				hJpsiMc->Fill(jpsiMc.Rapidity(), jpsiMc.Pt());
				hJpsi3DMc->Fill(jpsiMc.Rapidity(), vz, jpsiMc.Pt());

				if(mElectron->id>=0&&mElectron2->id>=0) {//good reconstructed electrons
					Double_t eta1 = mElectron->eta;
					Double_t dca1 = mElectron->dca;
					Double_t nsigma1 = myGaus_1->GetRandom();
					Double_t nHitsFit1 = mElectron->nFitPts;
					Double_t e1 = mElectron->energy;
					Double_t adc01 = mElectron->adc0;
					Double_t dsmAdc01 = mElectron->dsmAdc0;
					Double_t p1 = mElectron->p;
					Double_t pe1 = (e1>0)? p1/e1:9999;
					Double_t pt1 = mElectron->pt;
					Double_t nEta1 = mElectron->nEta;
					Double_t nPhi1 = mElectron->nPhi;
					Double_t zDist1 = mElectron->zDist;
					Double_t phiDist1 = mElectron->phiDist;
					bool isTPC1[4], isTrg1[4], isAdc01[4], isPE1[4], isNSMD1[4], isDist1[4];
					for(int iht=0;iht<4;iht++) {
						isTPC1[iht] = kFALSE;
						isTrg1[iht] = kFALSE;
						isAdc01[iht] = kFALSE;
						isPE1[iht] = kFALSE;
						isNSMD1[iht] = kFALSE;
						isDist1[iht] = kFALSE;
					}
					//TPC track
					if(eta1>=mTpceEtaCut[0]&&eta1<=mTpceEtaCut[1] &&
							nHitsFit1>=mTpceHitsFitCut &&
							dca1<=mTpceDcaCut &&
							nsigma1>mTpcenSigmaElectronCut[0]&&nsigma1<mTpcenSigmaElectronCut[1] &&
							pt1<30.&&p1<30.) {
						for(int iht=0;iht<4;iht++) {
							if(pt1>=mTpcePtCut[iht]&&p1>=mTpcePCut[iht]) isTPC1[iht] = kTRUE;
						}
					}
					//EMC Track
					if(nHitsFit1>=mTpceHitsFitCut &&
							dca1<=mEmceDcaCut &&
							eta1>=mEmceEtaCut[0]&&eta1<=mEmceEtaCut[1] &&
							nsigma1>mEmcenSigmaElectronCut[0]&&nsigma1<mEmcenSigmaElectronCut[1] &&
							pt1<30) {
						for(int iht=0;iht<4;iht++) {
							if(pt1>mEmcePtCut[iht]&&e1>0&&dsmAdc01>mDsmAdcCut[iht]) {
								isTrg1[iht] = kTRUE;
								if(adc01>mEmceAdcCut[iht]) {
									isAdc01[iht] = kTRUE;
									if(pe1>mEmcePECut[0]&&pe1<mEmcePECut[1]) {
										isPE1[iht] = kTRUE;
										if(nEta1>=mEmcenEtaCut&&nPhi1>=mEmcenPhiCut) {
											isNSMD1[iht] = kTRUE;
											if(zDist1>mEmceZDistCut[0]&&zDist1<mEmceZDistCut[1]&&
													phiDist1>mEmcePhiDistCut[0]&&phiDist1<mEmcePhiDistCut[1]) {
												isDist1[iht] = kTRUE;
											}
										}
									}
								}
							}
						}//for loop
					}


					Double_t eta2 = mElectron2->eta;
					Double_t dca2 = mElectron2->dca;
					Double_t nsigma2 = myGaus->GetRandom();
					Double_t nHitsFit2 = mElectron2->nFitPts;
					Double_t e2 = mElectron2->energy;
					Double_t adc02 = mElectron2->adc0;
					Double_t dsmAdc02 = mElectron2->dsmAdc0;
					Double_t p2 = mElectron2->p;
					Double_t pe2 = (e2>0)? p2/e2:9999;
					Double_t pt2 = mElectron2->pt;
					Double_t nEta2 = mElectron2->nEta;
					Double_t nPhi2 = mElectron2->nPhi;
					Double_t zDist2 = mElectron2->zDist;
					Double_t phiDist2 = mElectron2->phiDist;
					bool isTPC2[4], isTrg2[4], isAdc02[4], isPE2[4], isNSMD2[4], isDist2[4];
					for(int iht=0;iht<4;iht++) {
						isTPC2[iht] = kFALSE;
						isTrg2[iht] = kFALSE;
						isAdc02[iht] = kFALSE;
						isPE2[iht] = kFALSE;
						isNSMD2[iht] = kFALSE;
						isDist2[iht] = kFALSE;
					}
					//TPC track
					if(eta2>=mTpceEtaCut[0]&&eta2<=mTpceEtaCut[1] &&
							nHitsFit2>=mTpceHitsFitCut &&
							dca2<=mTpceDcaCut &&
							nsigma2>mTpcenSigmaElectronCut[0]&&nsigma2<mTpcenSigmaElectronCut[1] &&
							pt2<30.&&p2<30.) {
						for(int iht=0;iht<4;iht++) {
							if(pt2>=mTpcePtCut[iht]&&p2>=mTpcePCut[iht]) isTPC2[iht] = kTRUE;
						}
					}
					//EMC Track
					if(nHitsFit2>=mTpceHitsFitCut &&
							dca2<=mEmceDcaCut &&
							eta2>=mEmceEtaCut[0]&&eta2<=mEmceEtaCut[1] &&
							nsigma2>mEmcenSigmaElectronCut[0]&&nsigma2<mEmcenSigmaElectronCut[1] &&
							pt2<30) {
						for(int iht=0;iht<4;iht++) {
							if(pt2>mEmcePtCut[iht]&&e2>0&&dsmAdc02>mDsmAdcCut[iht]) {
								isTrg2[iht] = kTRUE;
								if(adc02>mEmceAdcCut[iht]) {
									isAdc02[iht] = kTRUE;
									if(pe2>mEmcePECut[0]&&pe2<mEmcePECut[1]) {
										isPE2[iht] = kTRUE;
										if(nEta2>=mEmcenEtaCut&&nPhi2>=mEmcenPhiCut) {
											isNSMD2[iht] = kTRUE;
											if(zDist2>mEmceZDistCut[0]&&zDist2<mEmceZDistCut[1]&&
													phiDist2>mEmcePhiDistCut[0]&&phiDist2<mEmcePhiDistCut[1]) {
												isDist2[iht] = kTRUE;
											}
										}
									}
								}
							}
						}//for loop
					}

					// Electron/Positron pairs
					TLorentzVector ePos(0.,0.,0.,0.),eNeg(0.,0.,0.,0.);
					for(int iht=0;iht<4;iht++){
						if(isTPC1[iht]==kFALSE|| isTrg1[iht]==kFALSE || isAdc01[iht]==kFALSE || isPE1[iht]==kFALSE || isNSMD1[iht]==kFALSE || isDist1[iht]==kFALSE) continue; 	
						if(mElectron->geantId==2) {cout<< "===========>";
							ePos.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
						}
						else if(mElectron->geantId==3){cout<<"<========";
							   	eNeg.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
						}
					}
					for(int iht=0;iht<4;iht++){
						if(isTPC2[iht]==kFALSE || isTrg2[iht]==kFALSE || isAdc02[iht]==kFALSE && isPE2[iht]==kFALSE && isNSMD2[iht]==kFALSE && isDist2[iht]==kFALSE) continue;	
						if(mElectron2->geantId==2)
							ePos.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
						else if(mElectron2->geantId==3)
							eNeg.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
					}
					if(ePos.Pt() > 0 && eNeg.Pt() > 0) {jpsiRcRaw= ePos + eNeg;
					if (jpsiRcRaw.M() >= mJpsiMassCut && 
						jpsiRcRaw.Pt() >=mHadronPtCut[0] && jpsiRcRaw.Pt() <= mHadronPtCut[1] && 
						jpsiRcRaw.Eta() >= mHadronEtaCut[0] && jpsiRcRaw.Eta() <= mHadronEtaCut[1] ) jpsiRc = jpsiRcRaw;  
					
					TLorentzVector ePosRest = ePos;
					ePosRest.Boost(-jpsiRc.Px()/jpsiRc.E(),-jpsiRc.Py()/jpsiRc.E(),-jpsiRc.Pz()/jpsiRc.E());
					Double_t dtheta = jpsiRc.Angle(ePosRest.Vect());
					Double_t costheta = TMath::Cos(dtheta);         
					hHtJpsiPol->Fill(dtheta);
					hHtJpsiPolThetaPt->Fill(costheta,jpsiRc.Pt());
					}

					Double_t rcPt = jpsiRc.Pt();
					Double_t rcY = jpsiRc.Rapidity();
					Double_t rcM = jpsiRc.M();
					//ht0
					if((isTrg1[0]&&isTPC2[0])||(isTrg2[0]&&isTPC1[0])||(isTrg1[0]&&isTrg2[0])) {
						hHt0JpsiTrg->Fill(rcY,rcPt,weight1);
						//hHt0Jpsi3DTrg->Fill(rcY,vz,rcPt);
					}
					if((isAdc01[0]&&isTPC2[0])||(isAdc02[0]&&isTPC1[0])||(isAdc01[0]&&isAdc02[0])) {
						hHt0JpsiAdc0->Fill(rcY,rcPt,weight1);
						hHt0Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
					}
					if((isPE1[0]&&isTPC2[0])||(isPE2[0]&&isTPC1[0])||(isPE1[0]&&isPE2[0])) {
						hHt0JpsiPE->Fill(rcY,rcPt,weight1);
						hHt0Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
						hHt0JpsiMassPE->Fill(rcY,rcPt,rcM,weight1);
					}
					if((isNSMD1[0]&&isTPC2[0])||(isNSMD2[0]&&isTPC1[0])||(isNSMD1[0]&&isNSMD2[0])) {
						hHt0JpsiNSMD->Fill(rcY,rcPt,weight1);
						hHt0Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
					}
					if((isDist1[0]&&isTPC2[0])||(isDist2[0]&&isTPC1[0])||(isDist1[0]&&isDist2[0])) {
						hHt0JpsiDist->Fill(rcY,rcPt,weight1);
						hHt0JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
						hHt0Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
					}


					//ht2
					if((isTrg1[2]&&isTPC2[2])||(isTrg2[2]&&isTPC1[2])||(isTrg1[2]&&isTrg2[2])) {
						hHt2JpsiTrg->Fill(rcY,rcPt,weight1);
						//hHt2Jpsi3DTrg->Fill(rcY,vz,rcPt);
					}
					if((isAdc01[2]&&isTPC2[2])||(isAdc02[2]&&isTPC1[2])||(isAdc01[2]&&isAdc02[2])) {
						hHt2JpsiAdc0->Fill(rcY,rcPt,weight1);
						hHt2Jpsi3DAdc0->Fill(rcY,vz,rcPt,weight1);
					}
					if((isPE1[2]&&isTPC2[2])||(isPE2[2]&&isTPC1[2])||(isPE1[2]&&isPE2[2])) {

						hHt2JpsiPE_Or->Fill(rcY,rcPt);
						// hweightPt_2->Fill(rcPt,weight1);
						// hweightPt_0->Fill(mcPt,weight1);
						hHt2JpsiPE->Fill(rcY,rcPt,weight1);
						hHt2JpsiMassPE_1->Fill(rcY,rcPt,rcM,weight1);
						hHt2JpsiMassPE->Fill(rcY,rcPt,rcM);
						hHt2Jpsi3DPE->Fill(rcY,vz,rcPt,weight1);
					}
					if((isNSMD1[2]&&isTPC2[2])||(isNSMD2[2]&&isTPC1[2])||(isNSMD1[2]&&isNSMD2[2])) {
						hHt2JpsiNSMD->Fill(rcY,rcPt,weight1);
						hHt2Jpsi3DNSMD->Fill(rcY,vz,rcPt,weight1);
					}
					if((isDist1[2]&&isTPC2[2])||(isDist2[2]&&isTPC1[2])||(isDist1[2]&&isDist2[2])) {
						hHt2JpsiDist->Fill(rcY,rcPt,weight1);
						hHt2JpsiMassDist->Fill(rcY,rcPt,rcM,weight1);
						hHt2Jpsi3DDist->Fill(rcY,vz,rcPt,weight1);
					}

					//ht3
				}

				/*
				   if(mElectron->energy>0) {
				   cout<<"(pt, eta, phi)- = "<<mMcJpsi->eminusPt<<"  "<<mMcJpsi->eminusEta<<"  "<<mMcJpsi->eminusPhi<<endl;
				   cout<<"(EmcE, neta, nphi, zdist, phidist, mindist) = "<<mMcJpsi->eminusEmcE<<"  "<<mMcJpsi->eminusnEta<<"  "<<mMcJpsi->eminusnPhi<<"  "<<mMcJpsi->eminusZDist<<"  "<<mMcJpsi->eminusPhiDist<<"  "<<mMcJpsi->eminusMinDist<<endl;
				   }
				   */
			}//partner loop
		}//electron/positron loop
		hNJpsi->Fill(nJpsi);
	}//event loop

	char buf[1024];
	sprintf(buf,"%s_cent_%d_%d.root", OutFile, mCentCut1, mCentCut2);
	TFile *f = new TFile(buf,"recreate");
	f->cd();
	cout<<"writing to "<<buf<<endl;

	hMcJpsiPol->GetXaxis()->SetTitle("#theta");

	hMcJpsiPolThetaPt->GetXaxis()->SetTitle("cos(#theta)");
	hMcJpsiPolThetaPt->GetYaxis()->SetTitle("J/#psi Pt");
	hMcJpsiPolThetaPt->SetOption("COLZ");

	hHtJpsiPol->GetXaxis()->SetTitle("#theta");

	hHtJpsiPolThetaPt->GetXaxis()->SetTitle("cos(#theta)");
	hHtJpsiPolThetaPt->GetYaxis()->SetTitle("j/#psi Pt");
	hHtJpsiPolThetaPt->SetOption("COLZ");

	hMcJpsiPol->Write();
	hMcJpsiPolThetaPt->Write();
	hHtJpsiPol->Write();
	hHtJpsiPolThetaPt->Write();
	hNreal->Write();

	f->Close();
	cout<<"done"<<endl;
}

//=======================================
int getCentrality(int refmult) {
	int mCentrality = -1;
	int   cent[] = { 7,9,11,12,14,15,17,19,23};//run10 AuAu200
	if(     refmult <= cent[0]) mCentrality = 0;
	else if(refmult <= cent[1]) mCentrality = 1;
	else if(refmult <= cent[2]) mCentrality = 2;
	else if(refmult <= cent[3]) mCentrality = 3;
	else if(refmult <= cent[4]) mCentrality = 4;
	else if(refmult <= cent[5]) mCentrality = 5;
	else if(refmult <= cent[6]) mCentrality = 6;
	else if(refmult <= cent[7]) mCentrality = 7;
	else if(refmult <= cent[8]) mCentrality = 8;
	else                        mCentrality = 9;

	return mCentrality;
}

