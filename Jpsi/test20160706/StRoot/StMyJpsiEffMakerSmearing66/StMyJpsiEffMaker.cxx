/***************************************************************************
 *
 **************************************************************************/
#include "StMyJpsiEffMaker.h"

#include "StEventTypes.h"
#include "StEvent/StEvent.h"
#include "TChain.h"
#include "StMyElectronMaker/StMyElectron.h"
#include "StMyElectronMaker/StMyElectronEvent.h"
#include "TLorentzVector.h"
#include "StRtsTable.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TH3F.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TDatime.h"

#include "cuts.h"

#define EMASS 0.000511
#define nHitsFitCut 20

#define PMASS 0.938272
ClassImp(StMyJpsiEffMaker);

//_____________________________________________________________
StMyJpsiEffMaker::StMyJpsiEffMaker(const char *name, TChain *chain, Int_t uncertainty_init):StMaker("myJpsiEff",name)
{
	myChain = new TChain("mcT");
	myChain->Add(chain);
	evCnt = 1;
	myEvent = new StMyElectronEvent();
	myChain->SetBranchAddress("mcE",&myEvent);
	mEtaMin = -1.;
	mEtaMax = 1.;
	mdEta = 0.1;


	for(int i=0;i<20;i++)
		for(int j=0;j<6;j++){
			mTofEffParsPos[i][j] = 0.;
			mTofEffParsNeg[i][j] = 0.;
		}
	mRan = new TRandom3();
	uncertainty=uncertainty_init;

	LOG_DEBUG << "StMyJpsiEffMaker::ctor"  << endm;
}

//_____________________________________________________________
StMyJpsiEffMaker::~StMyJpsiEffMaker() 
{ }

//______________________________________________________________
void StMyJpsiEffMaker::Clear(Option_t* option)
{
}
//_____________________________________________________________
Int_t StMyJpsiEffMaker::Init()
{

	if(uncertainty>=1 && uncertainty<=11){
		mDoSmearing = true;
		mSmearingFac = 0.0065+0.0001*uncertainty;
	}

	if(uncertainty==12) Aplus = 0.3825;
	if(uncertainty==13) Aminus = 0.3825;
	if(uncertainty==14) Bplus = 0.1233;
	if(uncertainty==15) Bminus = 0.1233;
	if(uncertainty==16) Aplus = 0.5408;
	if(uncertainty==17) Aminus = 0.5408;
	if(uncertainty==18) Bplus = 0.1743;
	if(uncertainty==19) Bminus = 0.1743;

	if(uncertainty==20) mTpceHitsFitCut = 25;
	if(uncertainty==21) mTpceHitsFitCut = 19;
	if(uncertainty==22) mTpceHitsFitCut = 22;
	if(uncertainty==23) mTpceHitsFitCut = 18;

	if(uncertainty==24) meanplus = 1.;
	if(uncertainty==25) meanminus = 1.;
	if(uncertainty==26) sigmaplus = 1.;
	if(uncertainty==27) sigmaminus = 1.;
	if(uncertainty==28) POL = 1;

	if(uncertainty==29) dsmadcfactor = 1.05;
	if(uncertainty==30) dsmadcfactor = 0.95;

	if(uncertainty==31) meanbeta = 1.;
	if(uncertainty==32) meanbeta = -1.;
	if(uncertainty==33) sigmabeta = 1.;
	if(uncertainty==34) sigmabeta = -1.;

	if(uncertainty==35) mEmcePECut[0] = 0.2, mEmcePECut[1] = 2.2;

	if(uncertainty==36) mTpceHitsDedxCut = 15.;
	if(uncertainty==37) mTpceDcaCut = 3.; 

	TDatime* time = new TDatime();

	betarootfile = new TFile("/star/data01/pwg/siwei/Jpsi/TOF_1_beta_mean_sigma.root");
	betamean = (TH1F*)betarootfile->Get("Tof_mean");
	betasigma = (TH1F*)betarootfile->Get("Tof_sigma");
	nsigmarootfile = new TFile("Nsigma.root");
	mean = (TH1F*)nsigmarootfile->Get("mh1mean");
	sigma = (TH1F*)nsigmarootfile->Get("mh1sigma");

	if(POL==0)meanfit = new TF1("meanfit","[0]",0.,16.5);
	else meanfit = new TF1("meanfit","[0]+[1]*x",0.,16.5);
	if(POL==0)sigmafit = new TF1("sigmafit","[0]",0.,16.5);
	else sigmafit = new TF1("sigmafit","[0]+[1]*x",0.,16.5);
	mean->Fit("meanfit","0Q");
	sigma->Fit("sigmafit","Q0");

	double para1[2];
	para1[0] = meanfit->GetParameter(0)+meanplus*meanfit->GetParError(0)-meanminus*meanfit->GetParError(0);
	para1[1] = sigmafit->GetParameter(0)+sigmaplus*sigmafit->GetParError(0)-sigmaminus*sigmafit->GetParError(0);

	double para2[2];
	para2[0] = meanfit->GetParameter(0)+meanplus*meanfit->GetParError(0)-meanminus*meanfit->GetParError(0);
	para2[1] = sigmafit->GetParameter(0)+sigmaplus*sigmafit->GetParError(0)-sigmaminus*sigmafit->GetParError(0);

	betafit = new TF1("betafit","[0]",0,4);

	char buf[1024];
	sprintf(buf,"rootfile%d/%s_sys%d.root",time->GetDate(),"OutFile",uncertainty);
	f = new TFile(buf,"recreate");
	f->cd();

	myGaus = new TF1("myGaus","gaus",-6,6);
	myGaus->SetParameters(1,para2[0],para2[1]);

	myGaus_1 = new TF1("myGaus_1","gaus",-6,6);
	myGaus_1->SetParameters(1,para1[0],para1[1]);	

	betaGaus1 = new TF1("betaGaus1","gaus",0.9,1.1);
	betaGaus2 = new TF1("betaGaus2","gaus",0.9,1.1);

	function_sigma = new TF1("function_sigma","[0]+[1]*x+[2]*x*x",0,30);	   
	function_sigma->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);

	function_tofeff = new TF1("function_tofeff","[0]*exp(-pow([1]/x,[2]))",0,30);
	function_tofeff->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);

	ifstream inf("tofeff/Eminus_TofEff_all_7.root_tofeffEMCMat_err.txt");
	cout<<"e-"<<endl;
	for(int i=0;i<20;i++){
		for(int j=0;j<6;j++) inf>>mTofEffParsNeg[i][j];
		for(int j=0;j<6;j++) cout<<mTofEffParsNeg[i][j]<<",";
		cout<<endl;
	}
	cout<<endl;
	inf.close();

	cout<<"e+"<<endl;
	inf.open("tofeff/Eplus_TofEff_all_7.root_tofeffEMCMat_err.txt");
	for(int i=0;i<20;i++){
		for(int j=0;j<6;j++) inf>>mTofEffParsPos[i][j];
		for(int j=0;j<6;j++) cout<<mTofEffParsPos[i][j]<<",";
		cout<<endl;
	}
	cout<<endl;
	inf.close();

	//hMCElectronPt = new TH1D("mcElectronPt","input electron pt",300,0,30);
	testhist = new TH1F("test","test",30,0,30);

	hCommonhitsvsRCPt = new TH2D("hCommonhitsvsRCPt","commonhits vs RC pT;tpc commonHits;RC p_{T} (GeV/c)",50,0,50,300,0,300);
	hCommonhitsvsMCPt = new TH2D("hCommonhitsvsMCPt","commonhits vs MC pT;tpc commonHits;MC p_{T} (GeV/c)",50,0,50,300,0,300);

	hJpsiPtCosThetaInvM = new TH3F("hJpsiPtCosThetaInvM","J/#psi Pt; Cos(#theta); Invariant mass",120,0,30,10,-1,1,20,2,4);
	hJpsiPtCosThetaInvM->Sumw2();
	/*
	   hJpsiCosThetaInvMPt = new TH3F("hJpsiCosThetaInvMPt","hJpsiCosThetaInvMPt",40,-1,1,40,2,4,120,0,30);
	   hJpsiCosThetaInvMPtCS = new TH3F("hJpsiCosThetaInvMPtCS","hJpsiCosThetaInvMPt",40,-1,1,40,2,4,120,0,30);
	   hJpsiCosThetaInvMPt1 = new TH3F("hJpsiCosThetaInvMPt1","hJpsiCosThetaInvMPt1",40,-1,1,40,2,4,120,0,30);
	   hJpsiCosThetaInvMPtCS1 = new TH3F("hJpsiCosThetaInvMPtCS1","hJpsiCosThetaInvMPt1",40,-1,1,40,2,4,120,0,30);

	   hJpsiCosThetaInvMPt->Sumw2();
	   hJpsiCosThetaInvMPtCS->Sumw2();
	   hJpsiCosThetaInvMPt1->Sumw2();
	   hJpsiCosThetaInvMPtCS1->Sumw2();
	   */
	hJpsiPhiInvMPt = new TH3F("hJpsiPhiInvMPt","hJpsiPhiInvMPt",40,-1,1,40,2,4,120,0,30);
	hJpsiPhiInvMPtCS = new TH3F("hJpsiPhiInvMPtCS","hJpsiPhiInvMPt",40,-1,1,40,2,4,120,0,30);
	hJpsiPhiInvMPt1 = new TH3F("hJpsiPhiInvMPt1","hJpsiPhiInvMPt1",40,-1,1,40,2,4,120,0,30);
	hJpsiPhiInvMPtCS1 = new TH3F("hJpsiPhiInvMPtCS1","hJpsiPhiInvMPt1",40,-1,1,40,2,4,120,0,30);

	hJpsiPhiInvMPt->Sumw2();
	hJpsiPhiInvMPtCS->Sumw2();
	hJpsiPhiInvMPt1->Sumw2();
	hJpsiPhiInvMPtCS1->Sumw2();

	//	hMBdsmAdcInvMPt = new TH3F("hMBdsmAdcInvMPt","hMBdsmAdcInvMPt",65,0,65,40,2,4,120,0,30);
	//	hMBdsmAdcInvMPtBG = new TH3F("hMBdsmAdcInvMPtBG","hMBdsmAdcInvMPtBG",65,0,65,40,2,4,120,0,30);
	//	hMBAdcInvMPt = new TH3F("hMBAdcInvMPt","hMBAdcInvMPt",800,0,800,40,2,4,120,0,30);
	//	hMBAdcInvMPtBG = new TH3F("hMBAdcInvMPtBG","hMBAdcInvMPtBG",800,0,800,40,2,4,120,0,30);
	//	hMBdsmAdcInvMPt->Sumw2();
	//	hMBdsmAdcInvMPtBG->Sumw2();
	//	hMBAdcInvMPt->Sumw2();
	//	hMBAdcInvMPtBG->Sumw2();
	/*
	   hHT0dsmAdcInvMPt = new TH3F("hHT0dsmAdcInvMPt","hHT0dsmAdcInvMPt",65,0,65,40,2,4,120,0,30);
	   hHT0dsmAdcInvMPtBG = new TH3F("hHT0dsmAdcInvMPtBG","hHT0dsmAdcInvMPtBG",65,0,65,40,2,4,120,0,30);
	   hHT0AdcInvMPt = new TH3F("hHT0AdcInvMPt","hHT0AdcInvMPt",800,0,800,40,2,4,120,0,30);
	   hHT0AdcInvMPtBG = new TH3F("hHT0AdcInvMPtBG","hHT0AdcInvMPtBG",800,0,800,40,2,4,120,0,30);
	   hHT0dsmAdcInvMPt->Sumw2();
	   hHT0dsmAdcInvMPtBG->Sumw2();
	   hHT0AdcInvMPt->Sumw2();
	   hHT0AdcInvMPtBG->Sumw2();

	   hHT1dsmAdcInvMPt = new TH3F("hHT1dsmAdcInvMPt","hHT1dsmAdcInvMPt",65,0,65,40,2,4,120,0,30);
	   hHT1dsmAdcInvMPtBG = new TH3F("hHT1dsmAdcInvMPtBG","hHT1dsmAdcInvMPtBG",65,0,65,40,2,4,120,0,30);
	   hHT1AdcInvMPt = new TH3F("hHT1AdcInvMPt","hHT1AdcInvMPt",800,0,800,40,2,4,120,0,30);
	   hHT1AdcInvMPtBG = new TH3F("hHT1AdcInvMPtBG","hHT1AdcInvMPtBG",800,0,800,40,2,4,120,0,30);
	   hHT1dsmAdcInvMPt->Sumw2();
	   hHT1dsmAdcInvMPtBG->Sumw2();
	   hHT1AdcInvMPt->Sumw2();
	   hHT1AdcInvMPtBG->Sumw2();

	   hHT2dsmAdcInvMPt = new TH3F("hHT2dsmAdcInvMPt","hHT2dsmAdcInvMPt",65,0,65,40,2,4,120,0,30);
	   hHT2dsmAdcInvMPtBG = new TH3F("hHT2dsmAdcInvMPtBG","hHT2dsmAdcInvMPtBG",65,0,65,40,2,4,120,0,30);
	   hHT2AdcInvMPt = new TH3F("hHT2AdcInvMPt","hHT2AdcInvMPt",800,0,800,40,2,4,120,0,30);
	   hHT2AdcInvMPtBG = new TH3F("hHT2AdcInvMPtBG","hHT2AdcInvMPtBG",800,0,800,40,2,4,120,0,30);
	   hHT2dsmAdcInvMPt->Sumw2();
	   hHT2dsmAdcInvMPtBG->Sumw2();
	   hHT2AdcInvMPt->Sumw2();
	   hHT2AdcInvMPtBG->Sumw2();
	   */
	hJpsiCosThetaPhiPt1 = new TH3F("hJpsiCosThetaPhiPt1","hJpsiCosThetaPhiPt1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hMBJpsiCosThetaPhiPt1 = new TH3F("hMBJpsiCosThetaPhiPt1","hMBJpsiCosThetaPhiPt1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT0JpsiCosThetaPhiPt1 = new TH3F("hHT0JpsiCosThetaPhiPt1","hHT0JpsiCosThetaPhiPt1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT1JpsiCosThetaPhiPt1 = new TH3F("hHT1JpsiCosThetaPhiPt1","hHT1JpsiCosThetaPhiPt1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT2JpsiCosThetaPhiPt1 = new TH3F("hHT2JpsiCosThetaPhiPt1","hHT2JpsiCosThetaPhiPt1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);

	hJpsiCosThetaPhiPt1->Sumw2();
	hMBJpsiCosThetaPhiPt1->Sumw2();
	hHT0JpsiCosThetaPhiPt1->Sumw2();
	hHT1JpsiCosThetaPhiPt1->Sumw2();
	hHT2JpsiCosThetaPhiPt1->Sumw2();

	hJpsiCosThetaPhiPtCS1 = new TH3F("hJpsiCosThetaPhiPtCS1","hJpsiCosThetaPhiPtCS1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hMBJpsiCosThetaPhiPtCS1 = new TH3F("hMBJpsiCosThetaPhiPtCS1","hMBJpsiCosThetaPhiPtCS1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT0JpsiCosThetaPhiPtCS1 = new TH3F("hHT0JpsiCosThetaPhiPtCS1","hHT0JpsiCosThetaPhiPtCS1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT1JpsiCosThetaPhiPtCS1 = new TH3F("hHT1JpsiCosThetaPhiPtCS1","hHT1JpsiCosThetaPhiPtCS1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT2JpsiCosThetaPhiPtCS1 = new TH3F("hHT2JpsiCosThetaPhiPtCS1","hHT2JpsiCosThetaPhiPtCS1",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);

	hJpsiCosThetaPhiPtCS1->Sumw2();
	hMBJpsiCosThetaPhiPtCS1->Sumw2();
	hHT0JpsiCosThetaPhiPtCS1->Sumw2();
	hHT1JpsiCosThetaPhiPtCS1->Sumw2();
	hHT2JpsiCosThetaPhiPtCS1->Sumw2();
	Clear("");
	return kStOK;
}

//_____________________________________________________________
Int_t StMyJpsiEffMaker::InitRun(Int_t runnumber)
{
	return kStOK;
}

//_____________________________________________________________
Int_t StMyJpsiEffMaker::FinishRun(Int_t runnumber)
{

	return kStOK;
}


//-------------------------------------------------------------
Int_t StMyJpsiEffMaker::Finish()
{
	f->Write();
	f->Close();
	Clear("");
	return kStOK;
}
//_____________________________________________________________
/*!
 * This method is to obtain the btofCollection from StEvent.
 * If StEvent is in the chain, retrieve it; if no StEvent in the chain,
 * a new StEvent is created.
 */
//_____________________________________________________________
Int_t StMyJpsiEffMaker::Make()
{
	myChain->GetEntry(evCnt++);
	mRan->SetSeed(evCnt);
	if(!myEvent) return kStOk;
	if(myEvent->eventID()<=0) return kStOk;
	Double_t vz = myEvent->vertexZ();
	if(vz<mVzCut[0]||vz>mVzCut[1]) return kStOk;

	TLorentzVector JpsiMc(0.,0.,0.,0.), ePosMc(0.,0.,0.,0.), eNegMc(0.,0.,0.,0.);
	TLorentzVector JpsiRc(0.,0.,0.,0.), ePosRc(0.,0.,0.,0.), eNegRc(0.,0.,0.,0.);
	Int_t nJpsi = 0;
	for(int j=0;j<myEvent->nReal();j++){
		mElectron = (StMyElectron*) myEvent->real()->UncheckedAt(j);
		if(mElectron->pGeantId!=160) continue;
		if(mElectron->mcId<0) continue;
		hCommonhitsvsMCPt->Fill(mElectron->tpcCommonHits,mElectron->mcPt);
		hCommonhitsvsRCPt->Fill(mElectron->tpcCommonHits,mElectron->pt);
		bool tag = kFALSE;
		for(int k=0;k<myEvent->nReal();k++){
			mElectron2 =(StMyElectron*) myEvent->real()->UncheckedAt(k);
			if(mElectron2->pGeantId!=160) continue;
			if(mElectron2->mcId<0) continue;
			if(mElectron2->mcId==mElectron->mcId) continue;
			//			if((mElectron->geantId!=2 || mElectron2->geantId!=3) && (mElectron->geantId!=3 || mElectron2->geantId!=2)) continue;	
			if(mElectron->geantId==2 && mElectron2->geantId==3){;}
			else if(mElectron->geantId==3 && mElectron2->geantId==2){;}
			else {continue;}

			Double_t deta = mElectron->mcY - mElectron2->mcY;
			Double_t dphi = mElectron->mcPhi - mElectron2->mcPhi;
			while(dphi>2*TMath::Pi()) dphi -= 2.*TMath::Pi();
			while(dphi<0) dphi += 2.*TMath::Pi();
			while(dphi>TMath::Pi()) dphi = dphi -2*TMath::Pi();
			Double_t dReta = mElectron->eta - mElectron2->eta;
			Double_t dRphi = mElectron->phi - mElectron2->phi;
			while(dRphi>2*TMath::Pi()) dRphi-=2.*TMath::Pi();
			while(dRphi<0) dRphi += 2.*TMath::Pi();
			while(dRphi>TMath::Pi()) dRphi = dRphi - 2.*TMath::Pi();
			if(TMath::Abs(deta)<0.1 && TMath::Abs(dphi)<0.5 && mElectron2->pId!=mElectron->pId) tag = kTRUE;
		}
		if(tag) continue;
		testhist->Fill(0);

		for(int k=j+1; k<myEvent->nReal();k++){
			mElectron2 =(StMyElectron*) myEvent->real()->UncheckedAt(k);
			if(mElectron2->pGeantId!=160) continue;
			if(mElectron2->mcId<0) continue;
			if(mElectron2->mcId==mElectron->mcId) continue;
			if(mElectron2->pId!=mElectron->pId) continue;
			if(mElectron->geantId==2 && mElectron2->geantId==3){
				ePosMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
				eNegMc.SetPtEtaPhiM(mElectron2->mcPt,mElectron2->mcEta, mElectron2->mcPhi, EMASS);
				ePosRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
				eNegRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);	
			}
			else if(mElectron->geantId==3 && mElectron2->geantId==2){
				eNegMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
				ePosMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
				eNegRc.SetPtEtaPhiM(mElectron->pt, mElectron->eta, mElectron->phi, EMASS);
				ePosRc.SetPtEtaPhiM(mElectron2->pt, mElectron2->eta, mElectron2->phi, EMASS);
			}
			else {continue;}
			if(mDoSmearing){
				TRandom *rcRand1 = new TRandom();
				TRandom *rcRand2 = new TRandom();
				double rcPt1 = (mElectron->pt)*(1+rcRand1->Gaus(0,mSmearingFac*(mElectron->pt)));
				double rcPt2 = (mElectron2->pt)*(1+rcRand2->Gaus(0,mSmearingFac*(mElectron2->pt)));
				double mcPt1 = (mElectron->mcPt)*(1+rcRand1->Gaus(0,mSmearingFac*(mElectron->mcPt)));
				double mcPt2 = (mElectron2->mcPt)*(1+rcRand2->Gaus(0,mSmearingFac*(mElectron2->mcPt)));
				if(mElectron->geantId==2 && mElectron2->geantId==3){
//					ePosMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					ePosMc.SetPtEtaPhiM(mcPt1, mElectron->mcEta, mElectron->mcPhi, EMASS);
//					eNegMc.SetPtEtaPhiM(mElectron2->mcPt,mElectron2->mcEta, mElectron2->mcPhi, EMASS);
					eNegMc.SetPtEtaPhiM(mcPt2,mElectron2->mcEta, mElectron2->mcPhi, EMASS);
					ePosRc.SetPtEtaPhiM(rcPt1, mElectron->eta, mElectron->phi, EMASS);
					eNegRc.SetPtEtaPhiM(rcPt2, mElectron2->eta, mElectron2->phi, EMASS);	
				}
				else if(mElectron->geantId==3 && mElectron2->geantId==2){
//					eNegMc.SetPtEtaPhiM(mElectron->mcPt, mElectron->mcEta, mElectron->mcPhi, EMASS);
					eNegMc.SetPtEtaPhiM(mcPt1, mElectron->mcEta, mElectron->mcPhi, EMASS);
//					ePosMc.SetPtEtaPhiM(mElectron2->mcPt, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
					ePosMc.SetPtEtaPhiM(mcPt2, mElectron2->mcEta, mElectron2->mcPhi, EMASS);
					eNegRc.SetPtEtaPhiM(rcPt1, mElectron->eta, mElectron->phi, EMASS);
					ePosRc.SetPtEtaPhiM(rcPt2, mElectron2->eta, mElectron2->phi, EMASS);
				}
			}
			JpsiMc = ePosMc + eNegMc;
			if(ePosRc.Pt()>0 && eNegRc.Pt()>0) JpsiRc = ePosRc + eNegRc;
			nJpsi++;

			//	Double_t weight1 = 4.32*TMath::Power(1+(JpsiMc.Pt()/4.10)*(JpsiMc.Pt()/4.10), -6)*(JpsiMc.Pt())*TMath::Exp(-0.5*(JpsiMc.Rapidity()*JpsiMc.Rapidity())/(1.416*1.416));
			Double_t weight1 = (A+Aplus-Aminus)*TMath::Power(1+(JpsiMc.Pt()/(B+Bplus-Bminus))*(JpsiMc.Pt()/(B+Bplus-Bminus)), -6)*(JpsiMc.Pt());
			if(rapidity)weight1 = weight1*TMath::Exp(-0.5*(JpsiMc.Rapidity()*JpsiMc.Rapidity())/(1.416*1.416));


			float deltaeta = mElectron->mcEta -mElectron2->mcEta;
			float deltaphi = mElectron->mcPhi - mElectron2->mcPhi;
			while(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
			while(deltaphi<0) deltaphi += 2.*TMath::Pi();
			while(deltaphi>TMath::Pi()) deltaphi = deltaphi -2*TMath::Pi();
			double deltaR =0;
			if(JpsiMc.Rapidity()<mPairYCut[0] || JpsiMc.Rapidity()>mPairYCut[1]) continue;

			TLorentzVector Proton1(0.,0.,100.,100),Proton2(0.,0.,-100.,100);						
			TLorentzVector Zaxis(0.,0.,0.,0.),Yaxis(0.,0.,0.,0.),Xaxis(0.,0.,0.,0.);
			TVector3 XX(0.,0.,0.),YY(0.,0.,0.),ZZ(0.,0.,0.);
			TVector3 XXHX(0.,0.,0.),YYHX(0.,0.,0.),ZZHX(0.,0.,0.);

			Proton1.Boost(-JpsiMc.Px()/JpsiMc.E(),-JpsiMc.Py()/JpsiMc.E(),-JpsiMc.Pz()/JpsiMc.E());
			Proton2.Boost(-JpsiMc.Px()/JpsiMc.E(),-JpsiMc.Py()/JpsiMc.E(),-JpsiMc.Pz()/JpsiMc.E());

			YYHX = JpsiMc.Vect().Cross(Proton1.Vect());
			XXHX = YYHX.Cross(JpsiMc.Vect());

			Yaxis.SetPx(Proton1.Py()*Proton2.Pz()-Proton1.Pz()*Proton2.Py());
			Yaxis.SetPy(Proton1.Pz()*Proton2.Px()-Proton1.Px()*Proton2.Pz());
			Yaxis.SetPz(Proton1.Px()*Proton2.Py()-Proton1.Py()*Proton2.Px());

			ZZ = Proton1.Vect()*(1/(Proton1.Vect()).Mag())-Proton2.Vect()*(1/(Proton2.Vect()).Mag());

			YY = Proton1.Vect().Cross(Proton2.Vect());
			Xaxis = Proton1;
			Xaxis = Zaxis;
			XX = Proton1.Vect()*(1/(Proton1.Vect()).Mag())+Proton2.Vect()*(1/(Proton2.Vect()).Mag());

			TLorentzVector ePosMcRest = ePosMc;
			TLorentzVector ePosRcRest = ePosRc;
			ePosMcRest.Boost(-JpsiMc.Px()/JpsiMc.E(),-JpsiMc.Py()/JpsiMc.E(),-JpsiMc.Pz()/JpsiMc.E());
			ePosRcRest.Boost(-JpsiRc.Px()/JpsiRc.E(),-JpsiRc.Py()/JpsiRc.E(),-JpsiRc.Pz()/JpsiRc.E());
			Double_t dtheta = JpsiMc.Angle(ePosMcRest.Vect());
			Double_t costheta = TMath::Cos(dtheta);
			Double_t sintheta = TMath::Sin(dtheta);

			Double_t dtheta_CS = ZZ.Angle(ePosMcRest.Vect());
			Double_t dphi_CS = TMath::ATan2(ePosMcRest.Vect().Dot(YY.Unit()),ePosMcRest.Vect().Dot(XX.Unit()));
			Double_t dphi_HX = TMath::ATan2(ePosMcRest.Vect().Dot(YYHX.Unit()),ePosMcRest.Vect().Dot(XXHX.Unit()));

			Double_t PtEdge[7] = {0,2,3,4,6,8,14};
			Int_t npt;
			for(npt=1;npt<6;npt++){
				if(JpsiMc.Pt()>=PtEdge[npt] && JpsiMc.Pt()<PtEdge[npt+1]) {
					polarization = lambda[npt]+dopol*lambda_err[npt]; 	
					polarizationphi = lambdaphi[npt]+dopolphi*lambdaphi_err[npt];
					if(polarization<-1) polarization = -1.;
					if(polarization>1) polarization = 1.;
					if(polarizationphi<-1) polarizationphi = -1.;
					if(polarizationphi>1) polarizationphi = 1.;
				}
			}
			weight1 = weight1*(1+polarization*costheta*costheta+polarizationphi*sintheta*sintheta*TMath::Cos(2*dphi_HX));	

			hJpsiPtCosThetaInvM->Fill(JpsiRc.Pt(),TMath::Cos(dtheta),JpsiRc.M());
			cout<<"Jpsi M="<<JpsiMc.M()<<"   "<<"Jpsi pt = "<<JpsiMc.Pt()<<"   "<<"weight1="<<weight1<<endl;
			hJpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
			hJpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);

			if(mElectron->id>=0 && mElectron2->id>=0){
				bool Qualityflag[2] ={kFALSE, kFALSE};
				double eta1 = mElectron->eta;
				double phi1 = mElectron->phi;
				double dca1 = mElectron->dca;
				double nHitsFit1 = mElectron->nFitPts;
				double eta2 = mElectron2->eta;
				double phi2 = mElectron2->phi;
				double dca2 = mElectron2->dca;
				double nHitsFit2 = mElectron2->nFitPts;
			}
			if(mElectron->id>=0 && mElectron2->id>=0){
				Double_t eta1 = mElectron->eta;
				Double_t phi1 = mElectron->phi;
				Double_t dca1 = mElectron->dca;
				Double_t nHitsFit1 = mElectron->nFitPts;
				Double_t nMaxPts1 = mElectron->nMaxPts;
				Double_t e1 = mElectron->energy0;	
				Double_t adc01 = mElectron->adc0;
				Double_t dsmAdc01 = mElectron->dsmAdc0;
				Double_t p1 = mElectron->p;
				Double_t pe1 = (e1>0.1)? p1/e1:9999;
				Double_t pt1 = mElectron->pt;
				if(mDoSmearing) pt1=pt1*(1.+mRan->Gaus(0,mSmearingFac*pt1));	
				Double_t nEta1 = mElectron->nEta;
				Double_t nPhi1 = mElectron->nPhi;
				Double_t zDist1 = mElectron->zDist;
				Double_t phiDist1 = mElectron->phiDist;
				Double_t nHitsdedx1 = mElectron->nDedxPts;
				Double_t nsigma1 = myGaus_1->GetRandom();

				double polpara[2][2];
				if(POL==1){
					polpara[0][0]=meanfit->GetParameter(0);
					polpara[0][1]=meanfit->GetParameter(1);
					polpara[1][0]=sigmafit->GetParameter(0);
					polpara[1][1]=sigmafit->GetParameter(1);
					myGaus_1->SetParameters(1,polpara[0][0]+polpara[0][1]*pt1,polpara[1][0]+polpara[1][1]*pt1);
					cout<<"         "<<myGaus_1->GetParameter(1)<<"         "<<myGaus_1->GetParameter(2)<<endl;
					nsigma1 = myGaus_1->GetRandom();
				}
				bool isEmc1 = kFALSE,isTpc1[4],isTOF1 = kFALSE,isTrg1[4];
				for(int iht=0;iht<4;iht++) {
					isTrg1[iht] = kFALSE;
					isTpc1[iht] = kFALSE;
				}
				int charge1 = 0;
				if(mElectron->geantId==2) charge1 = 1;
				if(mElectron->geantId==3) charge1 = -1;
				double tofEff1 = getTOFeff(charge1, pt1, eta1);
				double beta1para[2][2];
				double beta1=0;
				double pEff1 =p1;
				if(p1>1.5) pEff1 =1.5;
				beta1para[0][0]=betamean->GetBinContent(betamean->FindBin(pEff1));
				beta1para[0][1]=betamean->GetBinError(betamean->FindBin(pEff1));
				beta1para[1][0]=betasigma->GetBinContent(betasigma->FindBin(pEff1));
				beta1para[1][1]=betasigma->GetBinError(betasigma->FindBin(pEff1));
				betaGaus1->SetParameters(1,beta1para[0][0]+meanbeta*beta1para[0][1],beta1para[1][0]+sigmabeta*beta1para[1][1]);
				beta1=betaGaus1->GetRandom();
				cout<<"beta1============"<<beta1<<endl;

				Double_t eta2 = mElectron2->eta;
				Double_t phi2 = mElectron2->phi;
				Double_t dca2 = mElectron2->dca;
				Double_t nHitsFit2 = mElectron2->nFitPts;
				Double_t nMaxPts2 = mElectron2->nMaxPts;
				Double_t e2 = mElectron2->energy0;
				Double_t adc02 = mElectron2->adc0;
				Double_t dsmAdc02 = mElectron2->dsmAdc0;
				Double_t p2 = mElectron2->p;
				Double_t pe2 = (e2>0.1)? p2/e2:9999;
				Double_t pt2 = mElectron2->pt;
				if(mDoSmearing) pt2 = pt2*(1+mRan->Gaus(0,mSmearingFac*pt2));
				Double_t nEta2 = mElectron2->nEta;
				Double_t nPhi2 = mElectron2->nPhi;
				Double_t zDist2 = mElectron2->zDist;
				Double_t phiDist2 = mElectron2->phiDist;
				Double_t nHitsdedx2 = mElectron2->nDedxPts;
				Double_t nsigma2 = myGaus->GetRandom();

				if(POL==1){
					myGaus->SetParameters(1,polpara[0][0]+polpara[0][1]*pt2,polpara[1][0]+polpara[1][1]*pt2);
					nsigma2=myGaus->GetRandom();
				}
				bool isTpc2[4], isEmc2 = kFALSE,isTOF2 = kFALSE,isTrg2[4];
				for(int iht=0;iht<4;iht++){
				   	isTrg2[iht] = kFALSE;
					isTpc2[iht] = kFALSE;
				}
				int charge2 = 0;
				if(mElectron2->geantId==2) charge2 = 1;
				if(mElectron2->geantId==3) charge2 = -1;
				double tofEff2 = getTOFeff(charge2, pt2, eta2);
				double beta2para[2][2];
				double beta2=0;
				double pEff2 = p2;
				if(p2>1.5) pEff2 =1.5;
				beta2para[0][0]=betamean->GetBinContent(betamean->FindBin(pEff2));
				beta2para[0][1]=betamean->GetBinError(betamean->FindBin(pEff2));
				beta2para[1][0]=betasigma->GetBinContent(betasigma->FindBin(pEff2));
				beta2para[1][1]=betasigma->GetBinError(betasigma->FindBin(pEff2));
				betaGaus2->SetParameters(1,beta2para[0][0]+meanbeta*beta2para[0][1],beta2para[1][0]+sigmabeta*beta2para[1][1]);
				beta2=betaGaus2->GetRandom();
				cout<<"beta2========"<<beta2<<endl;

				if(nHitsFit1>=mTpceHitsFitCut &&
						nHitsFit1/nMaxPts1>=mTpceHitsRatio &&
						dca1<mTpceDcaCut &&
						eta1>=mTpceEtaCut[0] && eta1<=mTpceEtaCut[1] &&
						nsigma1>mTpceLoosenSigmaElectronCut[0] && nsigma1<mTpceLoosenSigmaElectronCut[1] &&
						nHitsdedx1>=mTpceHitsDedxCut && 
						mElectron->tpcCommonHits>=10 && 
						pt1<30. &&
						pt1>0.2){
					for(int iht=0;iht<4;iht++){					
						if(pt1>=mTpcePtCut[iht] && p1>=mTpcePCut[iht])isTpc1[iht] = kTRUE;
					}
						testhist->Fill(1);
					if(pe1>0.3 && pe1<1.5 && pt1>mEmcePtMin && nsigma1>=mTpcenSigmaElectronCut[0] && nsigma1<=mTpcenSigmaElectronCut[1]) {
						isEmc1 = kTRUE;
						testhist->Fill(6);	
					}
					if(mRan->Uniform(0,1)<tofEff1 && beta1>=mTpceBetaCut[0] && beta1<=mTpceBetaCut[1] && nsigma1>mTpcenSigmaElectronCut[0] && nsigma1<mTpcenSigmaElectronCut[1]){
						isTOF1 = kTRUE;
						testhist->Fill(7);
					}
					if(pt1>2.5 && dsmAdc01>11 && e1>0 && adc01>mEmceAdcCut[0]*dsmadcfactor && pe1>0.3 && pe1<1.5){
						isTrg1[0] = kTRUE;
						testhist->Fill(9);
					}
					if(pt1>3.6 && dsmAdc01>15 && e1>0 && adc01>mEmceAdcCut[1]*dsmadcfactor && pe1>0.3 && pe1<1.5){
						isTrg1[1] = kTRUE;
						testhist->Fill(10);
					}
					if(pt1>4.3 && dsmAdc01>18 && e1>0 && adc01>mEmceAdcCut[2]*dsmadcfactor && pe1>0.3 && pe1<1.5){
						isTrg1[2] = kTRUE;
						testhist->Fill(11);
					}
				}

				if(nHitsFit2>=mTpceHitsFitCut &&
						nHitsFit2/nMaxPts2>=mTpceHitsRatio &&
						dca2<mTpceDcaCut &&
						eta2>=mTpceEtaCut[0] && eta2<=mTpceEtaCut[1] &&
						nsigma2>=mTpceLoosenSigmaElectronCut[0] && nsigma2<=mTpceLoosenSigmaElectronCut[1] &&
						nHitsdedx2>=mTpceHitsDedxCut &&
						mElectron2->tpcCommonHits>=10 &&
						pt2<30. &&
						pt2>0.2){
					for(int iht=0;iht<4;iht++){
						if(pt2>=mTpcePtCut[iht] && p2>=mTpcePCut[iht])isTpc2[iht] = kTRUE;
					}
						testhist->Fill(2);
					if(pe2>0.3 && pe2<1.5 && pt2>mEmcePtMin && nsigma2>=mTpcenSigmaElectronCut[0] && nsigma2<=mTpcenSigmaElectronCut[1]) {
						isEmc2 = kTRUE;
							testhist->Fill(22);
					}
						if(mRan->Uniform(0,1)<tofEff2 && beta2>=mTpceBetaCut[0] && beta2<=mTpceBetaCut[1] && nsigma2>mTpcenSigmaElectronCut[0] && nsigma2<mTpcenSigmaElectronCut[1]){
						isTOF2 = kTRUE;
						testhist->Fill(21);
						}
					if(pt2>2.5 && dsmAdc02>11 && e2>0 && adc02>mEmceAdcCut[0]*dsmadcfactor && pe2>0.3 && pe2<1.5){ 
						isTrg2[0] = kTRUE;
						testhist->Fill(12);
					}
		
					if(pt2>3.6 && dsmAdc02>15 && e2>0 && adc02>mEmceAdcCut[1]*dsmadcfactor && pe2>0.3 && pe2<1.5){
						   	isTrg2[1] = kTRUE;
							testhist->Fill(13);
						}
					if(pt2>4.3 && dsmAdc02>18 && e2>0 && adc02>mEmceAdcCut[2]*dsmadcfactor && pe2>0.3 && pe2<1.5) {
						isTrg2[2] = kTRUE;
						testhist->Fill(14);
					}
					}		

				testhist->Fill(26);	
					if(JpsiRc.M()>3.0 && JpsiRc.M()<3.2){
						if((isTpc1[0] && isTpc2[0]) || (isTpc2[0] && isEmc1) || (isTpc1[0] && isEmc2) || (isEmc1 && isEmc2)) {
							hMBJpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
							hMBJpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);
						}

						if((isTrg1[0] && isEmc1 && isTpc2[0]) || (isTrg2[0] && isEmc2 && isTpc1[0])) {
							cout<<"mcPt = "<<JpsiMc.Pt()<<"weight1="<<weight1<<" mass = "<<JpsiMc.M()<<endl;
							hHT0JpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
							hHT0JpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);
							testhist->Fill(27);
						}
						if((isTrg1[1] && isEmc1 && isTpc2[1])||(isTrg2[1] && isEmc2 && isTpc1[1])) {
							hHT1JpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
							hHT1JpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);
						}
						if((isEmc1 && isTpc2[2] && isTrg1[2])||(isEmc2 && isTpc1[2] && isTrg2[2])) {
							hHT2JpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
							hHT2JpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);
						}
					}
				}
			}
			}
			return kStOk;
			}

			Double_t StMyJpsiEffMaker::getTOFeff(int charge, double pt, double eta){
				int ieta = (eta-mEtaMin)/mdEta;
				if(eta>mEtaMin && eta<mEtaMax){
					if(ieta<0||ieta>20) cout<<"WARN: eta bin is not within [-1,1]"<<endl;
					if(charge==1){
						function_tofeff->SetParameters(mTofEffParsPos[ieta][0]+tofmatching*mTofEffParsPos[ieta][1],mTofEffParsPos[ieta][2]+tofmatching*mTofEffParsPos[ieta][3],mTofEffParsPos[ieta][4]+tofmatching*mTofEffParsPos[ieta][5]);
					}else if(charge==-1){
						function_tofeff->SetParameters(mTofEffParsNeg[ieta][0]+tofmatching*mTofEffParsNeg[ieta][1],mTofEffParsNeg[ieta][2]+tofmatching*mTofEffParsNeg[ieta][3],mTofEffParsNeg[ieta][4]+tofmatching*mTofEffParsNeg[ieta][5]);
					}

					if(charge==1 || charge==-1){
						return function_tofeff->Eval(pt);
					}else{
						return 0.;
					}
				}else{
					return 0.;
				}
			}
