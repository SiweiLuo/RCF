#ifndef StMyJpsiEffMaker_H
#define StMyJpsiEffMaker_H

/***************************************************************************
 *
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/
#include "StMaker.h"

#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TChain;
class TLorentzVector;
class StMyElectronEvent;
class StMyElectron;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
class TF1;
class TProfile;
class TRandom3;

class StMyJpsiEffMaker:public StMaker
{
	public:

		/// Default constructor
		StMyJpsiEffMaker(const char *name="myJpsiEff", TChain *chain=0, Int_t uncertainty=0);

		~StMyJpsiEffMaker() ;

		void   Clear(Option_t* option="");
		Int_t  Init();
		Int_t  InitRun(Int_t);
		Int_t  FinishRun(Int_t);
		Int_t  Finish();
		Int_t  Make();
		Int_t  getCentrality(Int_t);
		Double_t getTOFeff(Int_t, Double_t, Double_t);
	private:
		TChain *myChain;
		Int_t  evCnt;
		StMyElectronEvent *myEvent;
		StMyElectron *mElectron;
		StMyElectron *mElectron2;
		StMyElectron *allElectron;
		Int_t uncertainty;

		TRandom3 *mRan;
		TRandom3 *mRan1;

		Double_t mEtaMin;
		Double_t mEtaMax;
		Double_t mdEta;
		Double_t mTofEffParsPos[20][6];
		Double_t mTofEffParsNeg[20][6];
		
		TFile *f;
		TFile *nsigmarootfile;
		TFile *betarootfile;

		TH1F *hnFitPts;
		TH1F *hnFitPtsOnFitPtsMax;
		TH1F *hDCA;
		TH1F *hETA;
		TH1F *hPt;
		TH1F *hP;


		TH1F *hMcJpsiPol;
		TH2F *hMcJpsiPolThetaPt;
		TH1F *hMcJpsiCosTheta;
		TH2F *hMcJpsiCosThetaPt;
		TH2F *hMcJpsiMassPt;

		TH2F *hMcJpsiThetaPt;
		TH2F *hRcJpsiThetaPt;
		TH2F *hEidJpsiThetaPt;

		TH2F *hMcJpsiPhiPt;
		TH2F *hRcJpsiPhiPt;
		TH2F *hEidJpsiPhiPt;
		TH2F *hMBJpsiPhiPt;
		TH2F *hHt0JpsiPhiPt;
		TH2F *hHt1JpsiPhiPt;
		TH2F *hHt2JpsiPhiPt;

		TH2F *hMcJpsiPhiPt1;
		TH2F *hRcJpsiPhiPt1;
		TH2F *hEidJpsiPhiPt1;
		TH2F *hMBJpsiPhiPt1;
		TH2F *hHt0JpsiPhiPt1;
		TH2F *hHt1JpsiPhiPt1;
		TH2F *hHt2JpsiPhiPt1;

		TH2F *hMcJpsiCosThetaPt1;
		TH2F *hRcJpsiCosThetaPt1;
		TH2F *hEidJpsiCosThetaPt1;
		TH2F *hMBJpsiCosThetaPt1;
		TH2F *hHt0JpsiCosThetaPt1;
		TH2F *hHt1JpsiCosThetaPt1;
		TH2F *hHt2JpsiCosThetaPt1;

		TH2F *hRcJpsiMassPt;

		TH1F *hHt0JpsiPol;
//		TH2F *hHt0JpsiPolThetaPt;
//		TH2F *hHt0JpsiPolThetaPt1;
//		TH2F *hHt0JpsiPolThetaPt2;
//		TH2F *hHt0JpsiPolThetaPt3;
//		TH2F *hHt0JpsiPolThetaPt4;
		TH2F *hMBJpsiCosThetaPt;
		TH2F *hHt0JpsiCosThetaPt;
		TH2F *hHt1JpsiCosThetaPt;
		TH2F *hMBJpsiThetaPt;
		TH2F *hHt0JpsiThetaPt;
		TH2F *hHt1JpsiThetaPt;

		TH2F *hMcJpsiThetaPt1;
		TH2F *hRcJpsiThetaPt1;
		TH2F *hEidJpsiThetaPt1;
		TH2F *hMBJpsiThetaPt1;
		TH2F *hHt0JpsiThetaPt1;
		TH2F *hHt1JpsiThetaPt1;
		TH2F *hHt2JpsiThetaPt1;

		TH2F *hHt0JpsiMassPt;
		TH2F *hHt2JpsiCosThetaPt;
		TH2F *hHt2JpsiMassPt;
		TH2F *hHt2JpsiThetaPt;

		TH2F *hEidJpsiCosThetaPt;
		
		TH1F *hHt0JpsiPolEff;
		TH2F *hHt0JpsiPolThetaPtEff;
		TH1F *hHt0JpsiCosTheta;
		TH2F *hHt0JpsiCosThetaPtEff;

		TH1F *hHt2JpsiPol;
		TH2F *hHt2JpsiPolThetaPt;
		TH1F *hHt2JpsiPolEff;

		TH1F *hMcVertexZ;
		TH2F *hMcVertexXY;
		TH1F *hRcVertexZ;
		TH1F *hVertexZdiff;
		TH1F *hRefMult;
		TH1F *hRefMultCut;
		TH2F *hMCdeltaEtavsdeltaPhi;
		TH2F *hRCdeltaEtavsdeltaPhi;
		TH2F *hRcJpsiCosThetaPt;

		TH1F *hNJpsi;
		TH1D *hMcJpsiPt;
		TH1D *hMcJpsiPt_Or;

		TH1F *hMcJpsiY;
		TH2D *hMcJpsiPtY;
		TH2D *hMcJpsiPtY_Or;
		TH1F *hMcJpsiPhi;
		TH2D *hMcJpsiMPt;
		TH2D *hMcJpsiMPt_1;

		TH1D *hMCElectronPt;

		TH2D *hmcPtvsrcPt;
		TH2D *hmcPtvsrcPt_Cut;

		TH1D *hRcJpsiPt;

		TH1D *hRcJpsiPt_Or;
		TH1F *hRcJpsiY;
		TH2D *hRcJpsiPtY;
		TH2D *hRcJpsiPtY_Or;
		TH2D *hRcJpsiPtY_B;
		TH2D *hRcJpsiPtY_BOr;
		TH1F *hRcJpsiPhi;
		TH2D *hRcJpsiMPt;
		TH2D *hRcJpsiPtDiff;
		TH2D *hRcJpsiPtDiff_rc;

		TH3F *hRcJpsiYDiff;
		TH3F *hRcJpsiPhiDiff;

		TH2D *hJpsiMc;

		TH2D *hHt0JpsiTrg;
		TH2D *hHt0JpsiAdc0;
		TH2D *hHt0JpsiPE;
		TH2D *hHt0JpsiNSMD;
		TH2D *hHt0JpsiDist;

		TH2D *hHt2JpsiTrg;
		TH2D *hHt2JpsiAdc0;
		TH2D *hHt2JpsiPE;
		TH2D *hHt2JpsiPE_Or;

		TH2D *hHt2JpsiNSMD;
		TH2D *hHt2JpsiDist;

		TH3F *hHt0JpsiMassPE;
		TH3F *hHt0JpsiMassDist;
		TH3F *hHt2JpsiMassPE;
		TH3F *hHt2JpsiMassPE_1;
		TH3F *hHt2JpsiMassPe;
		TH3F *hHt2JpsiMassDist;

		TH3F *hJpsi3DMc;

		TH3F *hHt0Jpsi3DAdc0;
		TH3F *hHt0Jpsi3DPE;
		TH3F *hHt0Jpsi3DNSMD;
		TH3F *hHt0Jpsi3DDist;

		TH3F *hHt2Jpsi3DAdc0;
		TH3F *hHt2Jpsi3DPE;
		TH3F *hHt2Jpsi3DNSMD;
		TH3F *hHt2Jpsi3DDist;

		TH2D *hCommonhitsvsRCPt;
		TH2D *hCommonhitsvsMCPt;

		TH2D *hdeltaEtavsMCPt;
		TH2D *hdeltaEtavsRCPt;
		TH2D *hdeltaPhivsMCPt;
		TH2D *hdeltaPhivsRCPt;
		TH2D *hdeltaRvsRCPt;
		TH2D *hdeltaRvsMCPt;

		TH2F *hHt2Adc0vsPt;
		TH2F *hHt2Adc0vsrcPt;
		TH2F *hHt2Adc0vsPt_1;
		TH3F *hHt2Adc0vsPtY;
		TH2F *hTest;
		TH1F *hPt1;
		TH1F *hPt2;
		TH1F *hEta1;
		TH1F *hEta2;
		TH1F *hPhi1;
		TH1F *hPhi2;

		TH2F *hMcPositronPhiPt;
		TH2F *hMBPositronPhiPt;
		TH2F *hHt0PositronPhiPt;
		TH2F *hHt1PositronPhiPt;
		TH2F *hHt2PositronPhiPt;
	
		TH2F *hMcPositronThetaPtCS;
		TH2F *hMBPositronThetaPtCS;
		TH2F *hHt0PositronThetaPtCS;
		TH2F *hHt1PositronThetaPtCS;
		TH2F *hHt2PositronThetaPtCS;

		TH2F *hMcPositronPhiPtCS;
		TH2F *hMBPositronPhiPtCS;
		TH2F *hHt0PositronPhiPtCS;
		TH2F *hHt1PositronPhiPtCS;
		TH2F *hHt2PositronPhiPtCS;

		TH2F *hMcJpsiThetaPtCS;
		TH2F *hRcJpsiThetaPtCS;
		TH2F *hEidJpsiThetaPtCS;
		TH2F *hMBJpsiThetaPtCS;
		TH2F *hHt0JpsiThetaPtCS;
		TH2F *hHt1JpsiThetaPtCS;
		TH2F *hHt2JpsiThetaPtCS;

		TH2F *hMcJpsiThetaPtCS1;
		TH2F *hRcJpsiThetaPtCS1;
		TH2F *hEidJpsiThetaPtCS1;
		TH2F *hMBJpsiThetaPtCS1;
		TH2F *hHt0JpsiThetaPtCS1;
		TH2F *hHt1JpsiThetaPtCS1;
		TH2F *hHt2JpsiThetaPtCS1;

		TH2F *hMcJpsiCosThetaPtCS;
		TH2F *hRcJpsiCosThetaPtCS;
		TH2F *hEidJpsiCosThetaPtCS;
		TH2F *hMBJpsiCosThetaPtCS;
		TH2F *hHt0JpsiCosThetaPtCS;
		TH2F *hHt1JpsiCosThetaPtCS;
		TH2F *hHt2JpsiCosThetaPtCS;

		TH2F *hMcJpsiCosThetaPtCS1;
		TH2F *hRcJpsiCosThetaPtCS1;
		TH2F *hEidJpsiCosThetaPtCS1;
		TH2F *hMBJpsiCosThetaPtCS1;
		TH2F *hHt0JpsiCosThetaPtCS1;
		TH2F *hHt1JpsiCosThetaPtCS1;
		TH2F *hHt2JpsiCosThetaPtCS1;


		TH2F *hRcJpsiPhiPtCS;
		TH2F *hEidJpsiPhiPtCS;
		TH2F *hMcJpsiPhiPtCS;
		TH2F *hMBJpsiPhiPtCS;
		TH2F *hHt0JpsiPhiPtCS;
		TH2F *hHt1JpsiPhiPtCS;
		TH2F *hHt2JpsiPhiPtCS;

		TH2F *hRcJpsiPhiPtCS1;
		TH2F *hEidJpsiPhiPtCS1;
		TH2F *hMcJpsiPhiPtCS1;
		TH2F *hMBJpsiPhiPtCS1;
		TH2F *hHt0JpsiPhiPtCS1;
		TH2F *hHt1JpsiPhiPtCS1;
		TH2F *hHt2JpsiPhiPtCS1;


		TH2F *hMcJpsiEtaPt;
		TH2F *hRcJpsiEtaPt;
		TH2F *hEidJpsiEtaPt;
		TH2F *hMBJpsiEtaPt;
		TH2F *hHt0JpsiEtaPt;
		TH2F *hHt1JpsiEtaPt;
		TH2F *hHt2JpsiEtaPt;

		TH2F *hMcJpsiEtaPtCS;
		TH2F *hRcJpsiEtaPtCS;
		TH2F *hEidJpsiEtaPtCS;
		TH2F *hMBJpsiEtaPtCS;
		TH2F *hHt0JpsiEtaPtCS;
		TH2F *hHt1JpsiEtaPtCS;
		TH2F *hHt2JpsiEtaPtCS;


		TH2F *hMcJpsiEtaPt1;
		TH2F *hRcJpsiEtaPt1;
		TH2F *hEidJpsiEtaPt1;
		TH2F *hMBJpsiEtaPt1;
		TH2F *hHt0JpsiEtaPt1;
		TH2F *hHt1JpsiEtaPt1;
		TH2F *hHt2JpsiEtaPt1;

		TH2F *hMcJpsiEtaPtCS1;
		TH2F *hRcJpsiEtaPtCS1;
		TH2F *hEidJpsiEtaPtCS1;
		TH2F *hMBJpsiEtaPtCS1;
		TH2F *hHt0JpsiEtaPtCS1;
		TH2F *hHt1JpsiEtaPtCS1;
		TH2F *hHt2JpsiEtaPtCS1;

		TH2F *hMcJpsiRapidityPt;
		TH2F *hRcJpsiRapidityPt;
		TH2F *hEidJpsiRapidityPt;
		TH2F *hMBJpsiRapidityPt;
		TH2F *hHt0JpsiRapidityPt;
		TH2F *hHt1JpsiRapidityPt;
		TH2F *hHt2JpsiRapidityPt;

		TH2F *hMcJpsiRapidityPt1;
		TH2F *hRcJpsiRapidityPt1;
		TH2F *hEidJpsiRapidityPt1;
		TH2F *hMBJpsiRapidityPt1;
		TH2F *hHt0JpsiRapidityPt1;
		TH2F *hHt1JpsiRapidityPt1;
		TH2F *hHt2JpsiRapidityPt1;

		TH1F *hMcJpsiRapidity1;
		TH1F *hRcJpsiRapidity1;
		TH1F *hEidJpsiRapidity1;
		TH1F *hMBJpsiRapidity1;
		TH1F *hHt0JpsiRapidity1;
		TH1F *hHt1JpsiRapidity1;
		TH1F *hHt2JpsiRapidity1;

		TH3F *hJpsiPtCosThetaInvM;
		TH2F *hMcJpsiPtCosThetaUN;
		TH2F *hMBJpsiPtCosThetaUN;
		TH2F *hHt0JpsiPtCosThetaUN;
		TH2F *hHt1JpsiPtCosThetaUN;
		TH2F *hHt2JpsiPtCosThetaUN;

		TH2F *hMBJpsiPtCosThetaUN1;
		TH2F *hHt0JpsiPtCosThetaUN1;
		TH2F *hHt1JpsiPtCosThetaUN1;
		TH2F *hHt2JpsiPtCosThetaUN1;

		TH2F *McJpsiCosThetaPt;
		TH2F *RcJpsiCosThetaPt;
		TH2F *MBJpsiCosThetaPt;
		TH2F *HT0JpsiCosThetaPt;
		TH2F *HT1JpsiCosThetaPt;
		TH2F *HT2JpsiCosThetaPt;	

		TH1F *hEnergyP;
		TH1F *TOFtest;

		TH1F *hLongitudinal;
		TH1F *hTransverse;

		TH2F *hRcAdcPt;
		TH2F *hMBAdcPt;
		TH2F *hHt0AdcPt;
		TH2F *hHt1AdcPt;
		TH2F *hHt2AdcPt;

		TH1F *nsigmacuteff;
		TH1F *betamean;
		TH1F *betasigma;
		TH1F *mean;
		TH1F *sigma;

		TF1 *myGaus;
		TF1 *myGaus_1;
		TF1 *function_sigma;
		TF1 *function_tofeff;
		TF1 *meanfit;
		TF1 *sigmafit;
		TF1 *betafit;
		TF1 *betaGaus1;
		TF1 *betaGaus2;
	
		TH2F *AdcDsmAdc1;
		TH1F *AdcDsmAdc2;	

		TH3F *hMcJpsiCosThetaInMPt;
//		TH3F *hJpsiCosThetaInvPtBG;

		TH3F *hJpsiCosThetaInvMPt;
		TH3F *hJpsiCosThetaInvMPtCS;
		TH3F *hJpsiPhiInvMPt;
		TH3F *hJpsiPhiInvMPtCS;
		TH3F *hJpsiCosThetaInvMPt1;
		TH3F *hJpsiCosThetaInvMPtCS1;
		TH3F *hJpsiPhiInvMPt1;
		TH3F *hJpsiPhiInvMPtCS1;

		TH3F *hMBdsmAdcInvMPt;
		TH3F *hMBdsmAdcInvMPtBG;
		TH3F *hMBAdcInvMPt;
		TH3F *hMBAdcInvMPtBG;
	
		TH3F *hHT0dsmAdcInvMPt;
		TH3F *hHT0dsmAdcInvMPtBG;
		TH3F *hHT0AdcInvMPt;
		TH3F *hHT0AdcInvMPtBG;

		TH3F *hHT1dsmAdcInvMPt;
		TH3F *hHT1dsmAdcInvMPtBG;
		TH3F *hHT1AdcInvMPt;
		TH3F *hHT1AdcInvMPtBG;

		TH3F *hHT2dsmAdcInvMPt;
		TH3F *hHT2dsmAdcInvMPtBG;
		TH3F *hHT2AdcInvMPt;
		TH3F *hHT2AdcInvMPtBG;

		TH3F *hJpsiCosThetaPhiPt1;
		TH3F *hMBJpsiCosThetaPhiPt1;
		TH3F *hHT0JpsiCosThetaPhiPt1;
		TH3F *hHT1JpsiCosThetaPhiPt1;
		TH3F *hHT2JpsiCosThetaPhiPt1;

		TH3F *hJpsiCosThetaPhiPtCS1;
		TH3F *hMBJpsiCosThetaPhiPtCS1;
		TH3F *hHT0JpsiCosThetaPhiPtCS1;
		TH3F *hHT1JpsiCosThetaPhiPtCS1;
		TH3F *hHT2JpsiCosThetaPhiPtCS1;


		TH1F *testhist;

		/// cvs
		virtual const char *GetCVS() const
		{
			static const char cvs[]="Tag $Name:  $Id: built " __DATE__ " " __TIME__ ; return cvs;
		}

		ClassDef(StMyJpsiEffMaker, 1)    ///StMyJpsiEffMaker - class to fille the StEvent from DAQ reader
};

#endif
