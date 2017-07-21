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

Double_t CrystalBall2(Double_t *x, Double_t *par);

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

	TDatime* time = new TDatime();

	mDoSmearing = true;

	if(uncertainty==3) DeltaB = 0.0133904;
	if(uncertainty==4) DeltaB = -0.0133904;

	if(uncertainty==10) mTpceHitsDedxCut = 15.;
	if(uncertainty==11) mTpceDcaCut = 3.;

	if(uncertainty==12) mTpceHitsFitCut = 25;
	if(uncertainty==13) mTpceHitsFitCut = 19;
	if(uncertainty==14) mTpceHitsFitCut = 22;
	if(uncertainty==15) mTpceHitsFitCut = 18;

	if(uncertainty==16) dsmadcfactor = 0.05;
	if(uncertainty==17) dsmadcfactor = -0.05;

	if(uncertainty==18) deltameanbeta = 1.;
	if(uncertainty==19) deltameanbeta = -1.;
	if(uncertainty==20) deltasigmabeta = 1.;
	if(uncertainty==21) deltasigmabeta = -1.;

	if(uncertainty==22) mEmcePECut[0] =0.2, mEmcePECut[1] = 2.2;


	betarootfile = new TFile("/star/data01/pwg/siwei/Jpsi/TOF_1_beta_mean_sigma.root");
	betamean = (TH1F*)betarootfile->Get("Tof_mean");
	betasigma = (TH1F*)betarootfile->Get("Tof_sigma");
	betafit = new TF1("betafit","[0]",0,4);

	betaGaus1 = new TF1("betaGaus1","gaus",0.9,1.1);
	betaGaus2 = new TF1("betaGaus2","gaus",0.9,1.1);

	function_tofeff = new TF1("function_tofeff","[0]*exp(-pow([1]/x,[2]))",0,30);
	function_tofeff->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);

	double para[2];

	nsigmarootfile = new TFile("nsigmae_fit.root","read");
	meanfit = (TF1*)nsigmarootfile->Get("mean1");
	sigmafit = (TF1*)nsigmarootfile->Get("sigma1");
	para[0] = meanfit->GetParameter(0);
	para[1] = sigmafit->GetParameter(0);

	cout<<"nsigma parameters ==="<<para[0]<<"    "<<para[1]<<endl;

	if(uncertainty>=5 && uncertainty<=6){	
		meanfit = (TF1*)nsigmarootfile->Get(Form("mean%d",uncertainty-3));
		para[0] = meanfit->GetParameter(0);
	}
	if(uncertainty>=7 && uncertainty<=8){
		sigmafit = (TF1*)nsigmarootfile->Get(Form("sigma%d",uncertainty-5));
		para[1] = sigmafit->GetParameter(0);	
	}

	if(uncertainty==9){
		meanfit = (TF1*)nsigmarootfile->Get("mean4");
		sigmafit = (TF1*)nsigmarootfile->Get("sigma4");
	}

	char buf[1024];
	sprintf(buf,"rootfile%d/%s_sys%d.root",time->GetDate(),"OutFile",uncertainty);
	f = new TFile(buf,"recreate");
	f->cd();
	myGaus = new TF1("myGaus","gaus",-6,6);
	myGaus->SetParameters(1,para[0],para[1]);

	myGaus_1 = new TF1("myGaus_1","gaus",-6,6);
	myGaus_1->SetParameters(1,para[0],para[1]);


	function_sigma = new TF1("function_sigma","[0]+[1]*x+[2]*x*x",0,30);	   
	function_sigma->SetParameters(1.77250e-2,3.18836e-3,1.68829e-3);

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

	inf.open("ptsmearing/reso_fit.txt");
	double reso[2],resoErr[2];
	for(int i=0;i<2;i++){
		inf>>reso[i]>>resoErr[i];
	}
	inf.close();

	fReso = new TF1("fReso","sqrt([0]*[0]*x*x+[1]*[1])",0,20);
	fReso->SetParameters(reso);
	if(uncertainty==1){
		fReso->SetParameter(1,reso[0]+resoErr[0]);
	}
	else if(uncertainty==2){
		fReso->SetParameter(1,reso[0]-resoErr[0]);
	}
	fReso->SetNpx(1000);

	inf.open("ptsmearing/dpt_CBfit.txt");
	double pars[7],parsErr[7];
	for(int i=0;i<7;i++){
		inf>>pars[i]>>parsErr[i];
	}
	inf.close();

	fmomShape = new TF1("fmomShape", CrystalBall2, -1., 1., 7);
	fmomShape->SetParameters(pars);
	fmomShape->SetNpx(1000);

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

	hMBJpsiPtInvM = new TH2F("hMBJpsiPtInvM","hMBJpsiPtInvM;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT0JpsiPtInvM = new TH2F("hHT0JpsiPtInvM","hHT0JpsiPtInvM;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT1JpsiPtInvM = new TH2F("hHT1JpsiPtInvM","hHT1JpsiPtInvM;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT2JpsiPtInvM = new TH2F("hHT2JpsiPtInvM","hHT2JpsiPtInvM;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);

	hMBJpsiPtInvM->Sumw2();
	hHT0JpsiPtInvM->Sumw2();
	hHT1JpsiPtInvM->Sumw2();
	hHT2JpsiPtInvM->Sumw2();

	hMBJpsiPtInvM2 = new TH2F("hMBJpsiPtInvM2","hMBJpsiPtInvM2;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT0JpsiPtInvM2 = new TH2F("hHT0JpsiPtInvM2","hHT0JpsiPtInvM2;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT1JpsiPtInvM2 = new TH2F("hHT1JpsiPtInvM2","hHT1JpsiPtInvM2;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);
	hHT2JpsiPtInvM2 = new TH2F("hHT2JpsiPtInvM2","hHT2JpsiPtInvM2;p_{T} GeV/c;m_{ee} GeV/c^2",120,0,30,40,2,4);

	hMBJpsiPtInvM2->Sumw2();
	hHT0JpsiPtInvM2->Sumw2();
	hHT1JpsiPtInvM2->Sumw2();
	hHT2JpsiPtInvM2->Sumw2();


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

	// 2eID observation histograms    
	hJpsiCosThetaPhiPt2 = new TH3F("hJpsiCosThetaPhiPt2","hJpsiCosThetaPhiPt2",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hMBJpsiCosThetaPhiPt2 = new TH3F("hMBJpsiCosThetaPhiPt2","hMBJpsiCosThetaPhiPt2",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT0JpsiCosThetaPhiPt2 = new TH3F("hHT0JpsiCosThetaPhiPt2","hHT0JpsiCosThetaPhiPt2",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT1JpsiCosThetaPhiPt2 = new TH3F("hHT1JpsiCosThetaPhiPt2","hHT1JpsiCosThetaPhiPt2",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hHT2JpsiCosThetaPhiPt2 = new TH3F("hHT2JpsiCosThetaPhiPt2","hHT2JpsiCosThetaPhiPt2",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);

	hJpsiCosThetaPhiPt2->Sumw2();
	hMBJpsiCosThetaPhiPt2->Sumw2();
	hHT0JpsiCosThetaPhiPt2->Sumw2();
	hHT1JpsiCosThetaPhiPt2->Sumw2();
	hHT2JpsiCosThetaPhiPt2->Sumw2();
	// 2eID observation histograms

	ht0trigpt = new TH1F("ht0trigpt","ht0trigpt",120,0,30);	
	ht1trigpt = new TH1F("ht1trigpt","ht1trigpt",120,0,30);
	ht2trigpt = new TH1F("ht2trigpt","ht2trigpt",120,0,30);

	ht0trigpt->Sumw2();
	ht1trigpt->Sumw2();
	ht2trigpt->Sumw2();

	ht0partnerpt = new TH1F("ht0partnerpt","ht0partnerpt",120,0,30);	
	ht1partnerpt = new TH1F("ht1partnerpt","ht1partnerpt",120,0,30);
	ht2partnerpt = new TH1F("ht2partnerpt","ht2partnerpt",120,0,30);

	ht0partnerpt->Sumw2();
	ht1partnerpt->Sumw2();
	ht2partnerpt->Sumw2();

	ht0trigpoe = new TH1F("ht0trigpoe","ht0trigpoe",100,0,4);
	ht1trigpoe = new TH1F("ht1trigpoe","ht1trigpoe",100,0,4);
	ht2trigpoe = new TH1F("ht2trigpoe","ht2trigpoe",100,0,4);

	ht0trigpoe->Sumw2();
	ht1trigpoe->Sumw2();
	ht2trigpoe->Sumw2();

	ht0partnerpoe = new TH1F("ht0partnerpoe","ht0partnerpoe",100,0,4);
	ht1partnerpoe = new TH1F("ht1partnerpoe","ht1partnerpoe",100,0,4);
	ht2partnerpoe = new TH1F("ht2partnerpoe","ht2partnerpoe",100,0,4);

	ht0partnerpoe->Sumw2();
	ht1partnerpoe->Sumw2();
	ht2partnerpoe->Sumw2();

	ht0trigadc0dsmadc = new TH2F("ht0trigadc0dsmadc","ht0trigadc0dsmadc",800,0,800,100,0,40);
	ht1trigadc0dsmadc = new TH2F("ht1trigadc0dsmadc","ht1trigadc0dsmadc",800,0,800,100,0,40);
	ht2trigadc0dsmadc = new TH2F("ht2trigadc0dsmadc","ht2trigadc0dsmadc",800,0,800,100,0,40);
	ht0trigadc0dsmadc->Sumw2();
	ht1trigadc0dsmadc->Sumw2();
	ht2trigadc0dsmadc->Sumw2();

	ht0partneradc0dsmadc = new TH2F("ht0partneradc0dsmadc","ht0partneradc0dsmadc",800,0,800,100,0,40);
	ht1partneradc0dsmadc = new TH2F("ht1partneradc0dsmadc","ht1partneradc0dsmadc",800,0,800,100,0,40);
	ht2partneradc0dsmadc = new TH2F("ht2partneradc0dsmadc","ht2partneradc0dsmadc",800,0,800,100,0,40);
	ht0partneradc0dsmadc->Sumw2();
	ht1partneradc0dsmadc->Sumw2();
	ht2partneradc0dsmadc->Sumw2();




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
			//			Double_t dReta = mElectron->eta - mElectron2->eta;
			//			Double_t dRphi = mElectron->phi - mElectron2->phi;
			/*			while(dRphi>2*TMath::Pi()) dRphi-=2.*TMath::Pi();
						while(dRphi<0) dRphi += 2.*TMath::Pi();
						while(dRphi>TMath::Pi()) dRphi = dRphi - 2.*TMath::Pi();*/
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

			//			TRandom *rcRand1 = new TRandom();
			//			TRandom *rcRand2 = new TRandom();

			double rcPt1 = smearElecPt(mElectron->pt,fReso,fmomShape);
			double rcPt2 = smearElecPt(mElectron2->pt,fReso,fmomShape);

			double mcPt1 = mElectron->mcPt;
			double mcPt2 = mElectron2->mcPt;

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
			JpsiMc = ePosMc + eNegMc;
			if(ePosRc.Pt()>0 && eNegRc.Pt()>0) JpsiRc = ePosRc + eNegRc;
			nJpsi++;

			//	Double_t weight1 = 4.32*TMath::Power(1+(JpsiMc.Pt()/4.10)*(JpsiMc.Pt()/4.10), -6)*(JpsiMc.Pt())*TMath::Exp(-0.5*(JpsiMc.Rapidity()*JpsiMc.Rapidity())/(1.416*1.416));
			//			Double_t weight1 = (A+Aplus-Aminus)*TMath::Power(1+(JpsiMc.Pt()/(B+DeltaB))*(JpsiMc.Pt()/(B+DeltaB)), -6)*(JpsiMc.Pt());
			Double_t weight1 = A*TMath::Power(1+(JpsiMc.Pt()/(B+DeltaB))*(JpsiMc.Pt()/(B+DeltaB)), -6)*(JpsiMc.Pt());
			weight1 *= TMath::Exp(-0.5*(JpsiMc.Rapidity()*JpsiMc.Rapidity())/(1.416*1.416));

			/*
			   float deltaeta = mElectron->mcEta -mElectron2->mcEta;
			   float deltaphi = mElectron->mcPhi - mElectron2->mcPhi;
			   while(deltaphi>2*TMath::Pi()) deltaphi -= 2.*TMath::Pi();
			   while(deltaphi<0) deltaphi += 2.*TMath::Pi();
			   while(deltaphi>TMath::Pi()) deltaphi = deltaphi -2*TMath::Pi();
			   */
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
			//			Double_t sintheta = TMath::Sin(dtheta);

			Double_t dtheta_CS = ZZ.Angle(ePosMcRest.Vect());
			Double_t dphi_CS = TMath::ATan2(ePosMcRest.Vect().Dot(YY.Unit()),ePosMcRest.Vect().Dot(XX.Unit()));
			Double_t dphi_HX = TMath::ATan2(ePosMcRest.Vect().Dot(YYHX.Unit()),ePosMcRest.Vect().Dot(XXHX.Unit()));

			hJpsiPtCosThetaInvM->Fill(JpsiRc.Pt(),TMath::Cos(dtheta),JpsiRc.M());
			cout<<"Jpsi M="<<JpsiMc.M()<<"   "<<"Jpsi pt = "<<JpsiMc.Pt()<<"   "<<"weight1="<<weight1<<endl;
			hJpsiCosThetaPhiPt1->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
			hJpsiCosThetaPhiPtCS1->Fill(TMath::Cos(dtheta_CS),dphi_CS,JpsiMc.Pt(),weight1);
			/*
			   if(mElectron->id>=0 && mElectron2->id>=0){
			//				bool Qualityflag[2] ={kFALSE, kFALSE};
			double eta1 = mElectron->eta;
			//				double phi1 = mElectron->phi;
			double dca1 = mElectron->dca;
			double nHitsFit1 = mElectron->nFitPts;
			double eta2 = mElectron2->eta;
			double phi2 = mElectron2->phi;
			double dca2 = mElectron2->dca;
			double nHitsFit2 = mElectron2->nFitPts;
			}
			*/
			if(mElectron->id>=0 && mElectron2->id>=0){
				Double_t eta1 = mElectron->eta;
				//				Double_t phi1 = mElectron->phi;
				Double_t dca1 = mElectron->dca;
				Double_t nHitsFit1 = mElectron->nFitPts;
				Double_t nMaxPts1 = mElectron->nMaxPts;
				Double_t e1 = mElectron->energy0;	
				Double_t adc01 = mElectron->adc0;
				Double_t dsmAdc01 = mElectron->dsmAdc0;
				Double_t p1 = mElectron->p;
				Double_t pe1 = (e1>0.1)? p1/e1:9999;
				Double_t pt1 = rcPt1;
//				Double_t pt1 = mElectron->pt;		
				//				Double_t nEta1 = mElectron->nEta;
				//				Double_t nPhi1 = mElectron->nPhi;
				//				Double_t zDist1 = mElectron->zDist;
				//				Double_t phiDist1 = mElectron->phiDist;
				Double_t nHitsdedx1 = mElectron->nDedxPts;
				Double_t nsigma1 = myGaus_1->GetRandom();

				if(uncertainty==9) myGaus_1->SetParameters(1,meanfit->Eval(pt1),sigmafit->Eval(pt1));

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
				betaGaus1->SetParameters(1,beta1para[0][0]+deltameanbeta*beta1para[0][1],beta1para[1][0]+deltasigmabeta*beta1para[1][1]);
				beta1=betaGaus1->GetRandom();
				cout<<"beta1============"<<beta1<<endl;

				Double_t eta2 = mElectron2->eta;
				//				Double_t phi2 = mElectron2->phi;
				Double_t dca2 = mElectron2->dca;
				Double_t nHitsFit2 = mElectron2->nFitPts;
				Double_t nMaxPts2 = mElectron2->nMaxPts;
				Double_t e2 = mElectron2->energy0;
				Double_t adc02 = mElectron2->adc0;
				Double_t dsmAdc02 = mElectron2->dsmAdc0;
				Double_t p2 = mElectron2->p;
				Double_t pe2 = (e2>0.1)? p2/e2:9999;
				Double_t pt2 = rcPt2;
//				Double_t pt2 = mElectron2->pt;
				//				Double_t nEta2 = mElectron2->nEta;
				//				Double_t nPhi2 = mElectron2->nPhi;
				//				Double_t zDist2 = mElectron2->zDist;
				//				Double_t phiDist2 = mElectron2->phiDist;
				Double_t nHitsdedx2 = mElectron2->nDedxPts;
				Double_t nsigma2 = myGaus->GetRandom();

				if(uncertainty==9) myGaus->SetParameters(1,meanfit->Eval(pt2),sigmafit->Eval(pt2));

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
				betaGaus2->SetParameters(1,beta2para[0][0]+deltameanbeta*beta2para[0][1],beta2para[1][0]+deltasigmabeta*beta2para[1][1]);
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
						testhist->Fill(1);
						if(pt1>=mTpcePtCut[iht] && p1>=mTpcePCut[iht]){
							isTpc1[iht] = kTRUE;
							testhist->Fill(6);	
						}
					}
					if(pe1>mEmcePECut[0] && pe1<mEmcePECut[1] && pt1>mEmcePtMin && nsigma1>=mTpcenSigmaElectronCut[0] && nsigma1<=mTpcenSigmaElectronCut[1]) {
						isEmc1 = kTRUE;
					}
					if(mRan->Uniform(0,1)<tofEff1 && beta1>=mTpceBetaCut[0] && beta1<=mTpceBetaCut[1] && nsigma1>mTpcenSigmaElectronCut[0] && nsigma1<mTpcenSigmaElectronCut[1]){
						isTOF1 = kTRUE;
						testhist->Fill(7);
					}
					if(pt1>2.5 && dsmAdc01>11 && e1>0 && adc01>mEmceAdcCut[0]*dsmadcfactor && pe1>mEmcePECut[0] && pe1<mEmcePECut[1]){
						isTrg1[0] = kTRUE;
						testhist->Fill(9);
					}
					if(pt1>3.6 && dsmAdc01>15 && e1>0 && adc01>mEmceAdcCut[1]*dsmadcfactor && pe1>mEmcePECut[0] && pe1<mEmcePECut[1]){
						isTrg1[1] = kTRUE;
						testhist->Fill(10);
					}
					if(pt1>4.3 && dsmAdc01>18 && e1>0 && adc01>mEmceAdcCut[2]*dsmadcfactor && pe1>mEmcePECut[0] && pe1<mEmcePECut[1]){
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
						testhist->Fill(2);
						if(pt2>=mTpcePtCut[iht] && p2>=mTpcePCut[iht]){
							isTpc2[iht] = kTRUE;
							testhist->Fill(22);
						}
					}
					if(pe2>mEmcePECut[0] && pe2<mEmcePECut[1] && pt2>mEmcePtMin && nsigma2>=mTpcenSigmaElectronCut[0] && nsigma2<=mTpcenSigmaElectronCut[1]) {
						isEmc2 = kTRUE;
					}
					if(mRan->Uniform(0,1)<tofEff2 && beta2>=mTpceBetaCut[0] && beta2<=mTpceBetaCut[1] && nsigma2>mTpcenSigmaElectronCut[0] && nsigma2<mTpcenSigmaElectronCut[1]){
						isTOF2 = kTRUE;
						testhist->Fill(21);
					}
					if(pt2>2.5 && dsmAdc02>11 && e2>0 && adc02>mEmceAdcCut[0]*dsmadcfactor && pe2>mEmcePECut[0] && pe2<mEmcePECut[1]){ 
						isTrg2[0] = kTRUE;
						testhist->Fill(12);
					}

					if(pt2>3.6 && dsmAdc02>15 && e2>0 && adc02>mEmceAdcCut[1]*dsmadcfactor && pe2>mEmcePECut[0] && pe2<mEmcePECut[1]){
						isTrg2[1] = kTRUE;
						testhist->Fill(13);
					}
					if(pt2>4.3 && dsmAdc02>18 && e2>0 && adc02>mEmceAdcCut[2]*dsmadcfactor && pe2>mEmcePECut[0] && pe2<mEmcePECut[1]) {
						isTrg2[2] = kTRUE;
						testhist->Fill(14);
					}
				}		

				if((isTpc1[0] && isTpc2[0] && isTOF1) || (isTpc1[0] && isTpc2[0] && isTOF2) || (isTpc2[0] && isEmc1) || (isTpc1[0] && isEmc2) || (isEmc1 && isEmc2)) {
					hMBJpsiPtInvM->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}

				if((isTrg1[0] && isEmc1 && isTpc2[0]) || (isTrg2[0] && isEmc2 && isTpc1[0])) {
					hHT0JpsiPtInvM->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}

				if((isTrg1[1] && isEmc1 && isTpc2[1])||(isTrg2[1] && isEmc2 && isTpc1[1])) {
					hHT1JpsiPtInvM->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}

				if((isEmc1 && isTpc2[2] && isTrg1[2])||(isEmc2 && isTpc1[2] && isTrg2[2])) {
					hHT2JpsiPtInvM->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}

				if(isTOF1 && isTOF2){
					hMBJpsiPtInvM2->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}
				if((isTrg1[0] && isEmc1 && isTOF2) || (isTrg2[0] && isEmc2 && isTOF1) || (isTrg1[0] && isEmc1 && isEmc2) || (isTrg2[0] && isEmc2 && isEmc1)){
					hHT0JpsiPtInvM2->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}
				if((isTrg1[1] && isEmc1 && isTOF2) || (isTrg2[1] && isEmc2 && isTOF1) || (isTrg1[1] && isEmc1 && isEmc2) || (isTrg2[1] && isEmc2 && isEmc1)){
					hHT1JpsiPtInvM2->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}
				if((isTrg1[2] && isEmc1 && isTOF2) || (isTrg2[2] && isEmc2 && isTOF1) || (isTrg1[2] && isEmc1 && isEmc2) || (isTrg2[2] && isEmc2 && isEmc1)){
					hHT2JpsiPtInvM2->Fill(JpsiRc.Pt(),JpsiRc.M(),weight1);
				}

				testhist->Fill(26);	
				if(JpsiRc.M()>3.0 && JpsiRc.M()<3.2){
					if((isTpc1[0] && isTpc2[0] && isTOF1) || (isTpc1[0] && isTpc2[0] && isTOF2) || (isTpc2[0] && isEmc1) || (isTpc1[0] && isEmc2) || (isEmc1 && isEmc2)) {
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

					// 2eID observation histogram 
					if((isTOF1 && isTOF2) || (isEmc1 && isEmc2) || (isEmc1 && isTOF2) || (isEmc2 && isTOF1)){
	//					hMBJpsiCosThetaPhiPt2->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
						hMBJpsiCosThetaPhiPt2->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
					}
					if((isTrg1[0] && isEmc1 && isTOF2) || (isTrg2[0] && isEmc2 && isTOF1) || (isTrg1[0] && isEmc1 && isEmc2) || (isTrg2[0] && isEmc2 && isEmc1)){
						hHT0JpsiCosThetaPhiPt2->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
						if(isTrg1[0]){
							ht0trigpt->Fill(pt1,weight1);
							ht0partnerpt->Fill(pt2,weight1);
							ht0trigadc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht0partneradc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht0trigpoe->Fill(pe1,weight1);
							ht0partnerpoe->Fill(pe2,weight1);
						}
						else if(isTrg2[0]){
							ht0trigpt->Fill(pt2,weight1);
							ht0partnerpt->Fill(pt1,weight1);
							ht0trigadc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht0partneradc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht0trigpoe->Fill(pe2,weight1);
							ht0partnerpoe->Fill(pe1,weight1);
						}
					}
					if((isTrg1[1] && isEmc1 && isTOF2) || (isTrg2[1] && isEmc2 && isTOF1) || (isTrg1[1] && isEmc1 && isEmc2) || (isTrg2[1] && isEmc2 && isEmc1)){
						hHT1JpsiCosThetaPhiPt2->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
						if(isTrg1[1]){
							ht1trigpt->Fill(pt1,weight1);
							ht1partnerpt->Fill(pt2,weight1);
							ht1trigadc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht1partneradc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht1trigpoe->Fill(pe1,weight1);
							ht1partnerpoe->Fill(pe2,weight1);
						}
						else if(isTrg2[1]){
							ht1trigpt->Fill(pt2,weight1);
							ht1partnerpt->Fill(pt1,weight1);
							ht1trigadc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht1partneradc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht1trigpoe->Fill(pe2,weight1);
							ht1partnerpoe->Fill(pe1,weight1);

						}
					}
					if((isTrg1[2] && isEmc1 && isTOF2) || (isTrg2[2] && isEmc2 && isTOF1) || (isTrg1[2] && isEmc1 && isEmc2) || (isTrg2[2] && isEmc2 && isEmc1)){
						hHT2JpsiCosThetaPhiPt2->Fill(costheta,dphi_HX,JpsiMc.Pt(),weight1);
						if(isTrg1[2]){
							ht2trigpt->Fill(pt1,weight1);
							ht2partnerpt->Fill(pt2,weight1);
							ht2trigadc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht2partneradc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht2trigpoe->Fill(pe1,weight1);
							ht2partnerpoe->Fill(pe2,weight1);
						}
						else if(isTrg2[2]){
							ht2trigpt->Fill(pt2,weight1);
							ht2partnerpt->Fill(pt1,weight1);
							ht2trigadc0dsmadc->Fill(adc02,dsmAdc02,weight1);
							ht2partneradc0dsmadc->Fill(adc01,dsmAdc01,weight1);
							ht2trigpoe->Fill(pe2,weight1);
							ht2partnerpoe->Fill(pe1,weight1);
						}
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

Double_t StMyJpsiEffMaker::smearElecPt(Double_t ept, TF1 *fdPtSig, TF1 *momShape){
	Double_t epsig = fdPtSig->Eval(ept);
	Double_t smpt = ept*(1. + momShape->GetRandom()*epsig/0.01);
	return smpt;
}

Double_t CrystalBall2(Double_t *x, Double_t *par)
{
	Double_t N = par[0];
	Double_t mu = par[1];
	Double_t s = par[2];
	Double_t n1 = par[3];
	Double_t alpha1 = par[4];
	Double_t n2 = par[5];
	Double_t alpha2 = par[6];

	Double_t A = TMath::Power(n1/fabs(alpha1), n1) * TMath::Exp(-alpha1*alpha1/2.);
	Double_t B = n1/fabs(alpha1) - fabs(alpha1);

	Double_t C = TMath::Power(n2/fabs(alpha2), n2) * TMath::Exp(-alpha2*alpha2/2.);
	Double_t D = n2/fabs(alpha2) - fabs(alpha2);

	Double_t norm = (x[0]-mu)/s;

	if(norm < -alpha1) {
		return N * A * TMath::Power(B-norm, -n1);
	} else if(norm < alpha2) {
		return N * TMath::Exp(-0.5*norm*norm);
	} else {
		return N * C * TMath::Power(D+norm, -n2);
	}
}
