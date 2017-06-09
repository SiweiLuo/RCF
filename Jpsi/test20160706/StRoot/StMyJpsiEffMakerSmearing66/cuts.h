const Double_t mVzCut[2] = {-50., 50.};

const Double_t mEmcePtCut[4] = {2.5,3.6,4.3,5.7};
const Double_t mEmceAdcCut[4] = {180,250,300, 400};
const Int_t mDsmAdcCut[4] = {11,15,18,25};

Double_t mEmcePECut[2] = {0.3,1.5};
const Double_t mEmcePtMin = 1.;

const Int_t mCommonHitCut = 10;
const Double_t mTpcePtCut[4] = {0.2,0.2,0.2,0.2};
const Double_t mTpcePCut[4] = {0.2,0.2,0.2,0.2};
const Double_t mTpceEtaCut[2] = {-1,1};
Double_t mTpceDcaCut = 1.;
Double_t mTpceHitsFitCut = 20;
Double_t mTpceHitsDedxCut = 11;
const Double_t mTpceHitsRatio = 0.52;

//const Double_t mTpceLoosenSigmaElectronCut[2] = {-3,3.};
const Double_t mTpceLoosenSigmaElectronCut[2] = {-2.5,3.};

const Double_t mTpcenSigmaElectronCut[2] = {-1.9,3};

const Double_t mTpceBetaCut[2] = {0.97,1.03};

const Double_t mPairYCut[2] = {-1,1};
const Double_t mJpsiMassCut = 1.0;

Double_t mSmearingFac = 0.;
bool mDoSmearing = false;

Double_t polarization = 0.;
Double_t polarizationphi = 0.;

Double_t A = 4.32, Aplus =0., Aminus =0.;
Double_t B = 4.10, Bplus =0., Bminus =0.;
const Int_t rapidity = 1;

Double_t meanplus = 0.,meanminus = 0.;
Double_t sigmaplus = 0.,sigmaminus = 0.;

const Double_t Htplus = 0.,Htminus = 0.;
const Double_t betaplus = 0., betaminus = 0.;
const Double_t tofmatching=0.;
Double_t meanbeta=0.,sigmabeta=0.;
Int_t POL=0;
Double_t dsmadcfactor=1.;

Double_t dopol = 0., dopolphi = 0.;
//const Double_t lambda[6] = {0.,-0.0931998 ,-0.126578,-0.199871,0.559828,0.349133};
//const Double_t lambda_err[6] = {0.,0.301887,0.243976,0.186846,0.381409,0.682504};
//const Double_t lambdaphi[6] = {0.,0.333039,0.441182,-0.472056,0.642627,10.3113};
//const Double_t lambdaphi_err[6] ={0.,0.547362,1.51316,0.288551,0.974175,1.52668};

const Double_t lambda[6] = {0.,0.,0.,0.,0.,0.};
const Double_t lambda_err[6] = {0.,0.,0.,0.,0.,0.};
const Double_t lambdaphi[6] = {0.,0.,0.,0.,0.,0.};
const Double_t lambdaphi_err[6] ={0.,0.,0.,0.,0.,0.};










