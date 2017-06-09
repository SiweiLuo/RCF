#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLatex.h"

#include "../ana/ptBins.h"
#include "TTimer.h"
#include "jpsi.h"

#define eMass 0.000511
#define jpsiMass 3.096916 

TRandom3 *mRan = 0;

TH1D *hMcMass;

//------------------------------------------------------
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


void getReso(){
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(1);

   TH2F *dPtvsmcPt = (TH2F*)fsim->Get("dPtvsmcPt");
   dPtvsmcPt->RebinY(10);
   dPtvsmcPt->RebinX(2);

   TFile *fdata = new TFile("../ana/out/jpsi_mass_HT2_2TrkPid.root");
   TH1D *hMAll = (TH1D*)fdata->Get("hJpsiMAll");
   TH1D *hm[nJpsiPt];
   for(int i=0;i<nJpsiPt;i++){
      hm[i] = (TH1D*)fdata->Get("hJpsiMAll");
   }

   TFile *fth = new TFile("../ana/out/jpsi_xsec_final.root");
   TGraphErrors *gCemZY = (TGraphErrors*)fth->Get("gCemZY");

   double dptMin = -0.02, dptMax = 0.02;
   int b1 = dPtvsmcPt->GetYaxis()->FindBin(dptMin);
   int b2 = dPtvsmcPt->GetYaxis()->FindBin(dptMax)-1;
   dPtvsmcPt->GetYaxis()->SetRange(b1,b2);
   dPtvsmcPt->FitSlicesY();
   TCanvas *c1 = new TCanvas();
   TH1D *hReso = (TH1D*)gDirectory->Get(Form("%s_2",dPtvsmcPt->GetName()));
   hReso->SetAxisRange(0,13,"X");
   hReso->SetAxisRange(0,0.05,"Y");
   hReso->Draw();
   TF1 *fReso = new TF1("fReso","sqrt([0]*[0]*x*x+[1]*[1])",hReso->GetXaxis()->GetXmin(),hReso->GetXaxis()->GetXmax());
   fReso->SetParameters(0.003,0.006);
   fReso->SetRange(0.2,15);
   hReso->Fit(fReso,"REI");
   hReso->Draw();
   fReso->SetLineColor(2);
   fReso->Draw("same");

   c1->Print("fig/reso_fit.pdf");

   ofstream outf("txt/reso_fit.txt");
   outf<<fReso->GetParameter(0)<<" \t "<<fReso->GetParError(0)<<" \t "<<fReso->GetParameter(1)<<" \t "<<fReso->GetParError(1)<<endl;
   outf<<"function sqrt([0]*[0]*x*x+[1]*[1])."<<endl;
   outf.close();
}

//------------------------------------------------------
void fitNormReso(int flag=1){

   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1);
   const char *fname = "OutFile_wTree_cent_0_9_resoNorm.root";
   TFile *fsim = new TFile(fname);   
   TH2F *dPtvsmcPtNorm = (TH2F*)fsim->Get("dPtvsmcPtNorm");
   TH1D *hResoNorm = (TH1D*)dPtvsmcPtNorm->ProjectionY("hResoNorm");
   TH1D *hResoNormByTree = (TH1D*)hResoNorm->Clone("hResoNormByTree");
   hResoNormByTree->Reset();
   TF1 *fCB2 = new TF1("fCB2", CrystalBall2, -1., 1., 7);
   fCB2->SetParNames("NS", "#mu", "#sigma_{p_{T}}/p_{T}", "n1", "#alpha1", "n2", "#alpha2");
   fCB2->SetParameters(1.e5, -0.01, 0.01, 1.5, 1.5, 1.5, 1.5 );
   fCB2->SetParLimits(1, -0.1, 0.1);
   fCB2->SetParLimits(2, 0, 1.);
   fCB2->SetParLimits(3, 0., 20.);
   fCB2->SetParLimits(4, 0., 20.);
   fCB2->SetParLimits(5, 0., 20.);
   fCB2->SetParLimits(6, 0., 20.);
   if(flag==2){
      ifstream inf("txt/reso_fit.txt");
      double reso[2],resoErr[2];
      for(int i=0;i<2;i++){
         inf>>reso[i]>>resoErr[i];
      }
      inf.close();

      TF1 *fReso = new TF1("fReso","sqrt([0]*[0]*x*x+[1]*[1])",0,20);
      fReso->SetParameters(reso);
      fReso->SetNpx(1000);

      TChain *chain = new TChain("jpsi");
      chain->Add(fname);
      jpsi *t = new jpsi(chain);
      for(int i=0;i<chain->GetEntries();++i){
         t->GetEntry(i);
         if(i%100000==0) cout<<i<<" th entry ..."<<endl;
         float mcPt = t->mcPt; 
         float mcEta = t->mcEta; 
         float mcPhi = t->mcPhi; 
         float mcM = t->mcM; 
         float mce1Pt = t->mce1Pt;
         float mce2Pt = t->mce2Pt;

         float rce1Pt = t->rce1Pt;
         float rce1Eta = t->rce1Eta;
         float rce1Phi = t->rce1Phi;
         float rce2Pt = t->rce2Pt;
         float rce2Eta = t->rce2Eta;
         float rce2Phi = t->rce2Phi;

         if(mce1Pt>0&&rce2Pt>0){
            float dpt = (rce2Pt-mce1Pt)/mce1Pt;
            float norm = fReso->Eval(mce1Pt);
            hResoNormByTree->Fill(dpt*0.01/norm);
         }
      }
   }
   TCanvas *c1 = new TCanvas();
   //for(int i=0;i<hResoNorm->GetNbinsX();i++){
   //   hResoNorm->SetBinError(i+1,0);
   //}
   if(flag==1){
      hResoNorm->UseCurrentStyle();
      hResoNorm->Rebin(5);
      hResoNorm->SetAxisRange(-0.2,0.2,"X");
      hResoNorm->SetAxisRange(1,1e7,"Y");
      hResoNorm->SetMarkerStyle(1);
      hResoNorm->Draw("e1");
      hResoNorm->Fit(fCB2,"R");
   }
   if(flag==2){
      hResoNormByTree->UseCurrentStyle();
      hResoNormByTree->SetAxisRange(-0.2,0.2,"X");
      hResoNormByTree->SetAxisRange(1,1e7,"Y");
      hResoNormByTree->SetMarkerStyle(1);
      hResoNormByTree->Fit(fCB2,"R");
   }
   gPad->SetLogy();
   fCB2->SetRange(-0.2,0.2);
   //fCB2->SetRange(-0.2,0.2);
   fCB2->Draw("same");
   c1->Print(Form("fig/dpt_CBfit_flag_%d.pdf",flag));

   ofstream outf("txt/dpt_CBfit.txt");
   outf<<fCB2->GetParameter(0)<<" \t "<<fCB2->GetParError(0)<<
      " \t "<<fCB2->GetParameter(1)<<" \t "<<fCB2->GetParError(1)<<
      " \t "<<fCB2->GetParameter(2)<<" \t "<<fCB2->GetParError(2)<<
      " \t "<<fCB2->GetParameter(3)<<" \t "<<fCB2->GetParError(3)<<
      " \t "<<fCB2->GetParameter(4)<<" \t "<<fCB2->GetParError(4)<<
      " \t "<<fCB2->GetParameter(5)<<" \t "<<fCB2->GetParError(5)<<
      " \t "<<fCB2->GetParameter(6)<<" \t "<<fCB2->GetParError(6)<<endl;
   outf<<"function double crystal ball par parErr"<<endl;
   outf.close();

}

//-------------------------------------------------------
TLorentzVector myBoost(TLorentzVector parent, TLorentzVector daughter){
   float betax = parent.Px()/parent.E();
   float betay = parent.Py()/parent.E();
   float betaz = parent.Pz()/parent.E();
   daughter.Boost(betax,betay,betaz);
   return daughter;
}


//-------------------------------------------------------
TLorentzVector twoBodyDecay(TLorentzVector parent, Double_t dmass) {

   Double_t e = parent.M()/2.;
   Double_t p = sqrt(e*e-dmass*dmass);
   Double_t costheta = mRan->Uniform(-1.,1.);
   Double_t phi = mRan->Uniform(0,TMath::Pi()*2);
   Double_t pz = p*costheta;
   Double_t px = p*sqrt(1.-costheta*costheta)*cos(phi);
   Double_t py = p*sqrt(1.-costheta*costheta)*sin(phi);
   TLorentzVector daughter(px,py,pz,e);
   return myBoost(parent,daughter);
}

void smearPt(int nRan, TH1D *rcPt, TH2F *simRcMPt, TF1 *fdPtSig, TF1 *momShape){
   for(int i=0;i<nRan;i++){
      Double_t rap = mRan->Uniform(-1.,1.); 
      Double_t phi = mRan->Uniform(0.,2.*TMath::Pi());
      Double_t pt  = rcPt->GetRandom();
      Double_t m   = jpsiMass;
      Double_t mt = sqrt(pt*pt+m*m);
      Double_t pz = mt*TMath::SinH(rap);
      Double_t ptot = sqrt(pt*pt+pz*pz);
      Double_t eta = 0.5*log((ptot+pz)/(ptot-pz+1e-20));


      if(ptot==0){
         cout<<"p==0 continue"<<endl;
         continue;
      }

      TLorentzVector parent;
      TLorentzVector parentSmear;
      parent.SetPtEtaPhiM(pt,eta,phi,m);

      TLorentzVector daughterP(0,0,0,0);
      TLorentzVector daughterN(0,0,0,0);
      // two body for omega, phi, jpsi
      daughterP = twoBodyDecay(parent,eMass);
      daughterN = parent-daughterP;

      Double_t eppt  = daughterP.Pt();
      Double_t epeta = daughterP.Eta();
      Double_t epphi = daughterP.Phi();

      Double_t epsig = 0.01;
      epsig = fdPtSig->Eval(eppt);
      Double_t smeppt  = eppt*(1. + momShape->GetRandom()*epsig/0.01);

      Double_t empt  = daughterN.Pt();
      Double_t emeta = daughterN.Eta();
      Double_t emphi = daughterN.Phi();

      Double_t emsig = 0.01;
      emsig = fdPtSig->Eval(empt);
      Double_t smempt  = empt*(1. + momShape->GetRandom()*emsig/0.01);

      TLorentzVector smdaughterP(0,0,0,0);
      TLorentzVector smdaughterN(0,0,0,0);
      smdaughterP.SetPtEtaPhiM(smeppt,epeta,epphi,eMass);
      smdaughterN.SetPtEtaPhiM(smempt,emeta,emphi,eMass);

      parentSmear = smdaughterP + smdaughterN;
      simRcMPt->Fill(parentSmear.M(),pt);
   }
}
void smearPt(int nRan, TH1D *rcPt, TH1D *simRcM, TH1D *simRcPt, TF1 *fdPtSig, TF1 *momShape){

   TH2F *h2D = new TH2F("h2D","h2D",simRcM->GetNbinsX(),simRcM->GetXaxis()->GetXmin(),simRcM->GetXaxis()->GetXmax(),simRcPt->GetNbinsX(),simRcPt->GetXaxis()->GetXmin(),simRcPt->GetXaxis()->GetXmin());
   smearPt(nRan,rcPt,h2D,fdPtSig,momShape);
   TH1D *hM = (TH1D*)h2D->ProjectionX("hM");
   simRcM->Reset();
   simRcM->Add(hM);
   hM->Reset();
   TH1D *hPt = (TH1D*)h2D->ProjectionY("hPt");
   simRcPt->Reset();
   simRcPt->Add(hPt);
   hPt->Reset();
}

void smearPtFromTree(TH2F *hRcMPt, TF1 *fdPtSig, TF1 *momShape, const char *fname, const char *tname = "jpsi"){

   TFile *fin = new TFile(fname);
   TChain *chain = new TChain(tname);
   chain->Add(fname);
   jpsi *t = new jpsi(chain);
   int n = chain->GetEntries();
   cout<<"nEntries = "<<n<<endl;
   //for(int i=0;i<10000;++i){
   for(int i=0;i<n;++i){
      t->GetEntry(i);
      if(i%100000==0) cout<<i<<" th entry ..."<<endl;
      float mcPt = t->mcPt; 
      float mcEta = t->mcEta; 
      float mcPhi = t->mcPhi; 
      float mcM = t->mcM; 
      float mce1Pt = t->mce1Pt;
      float mce1Eta = t->mce1Eta;
      float mce1Phi = t->mce1Phi;
      float mce2Pt = t->mce2Pt;
      float mce2Eta = t->mce2Eta;
      float mce2Phi = t->mce2Phi;

      float rce1Pt = t->rce1Pt;
      float rce1Eta = t->rce1Eta;
      float rce1Phi = t->rce1Phi;
      float rce2Pt = t->rce2Pt;
      float rce2Eta = t->rce2Eta;
      float rce2Phi = t->rce2Phi;

      Double_t e1sig = 0.01;
      e1sig = fdPtSig->Eval(mce1Pt);
      Double_t sme1pt  = mce1Pt*(1. + momShape->GetRandom()*e1sig/0.01);

      Double_t e2sig = 0.01;
      e2sig = fdPtSig->Eval(mce2Pt);
      Double_t sme2pt  = mce2Pt*(1. + momShape->GetRandom()*e2sig/0.01);

      TLorentzVector smdaughter1(0,0,0,0);
      TLorentzVector smdaughter2(0,0,0,0);
      smdaughter1.SetPtEtaPhiM(sme1pt,mce1Eta,mce1Phi,eMass);
      smdaughter2.SetPtEtaPhiM(sme2pt,mce2Eta,mce2Phi,eMass);
      TLorentzVector fourMomSmear = smdaughter1 + smdaughter2;
      if(rce1Pt>0&&rce2Pt>0){
         hRcMPt->Fill(fourMomSmear.M(),mcPt); 
      }
   }
}


void smearPtFromTree(TH1D *hRcM, TH1D *hRcPt, TF1 *fdPtSig, TF1 *momShape, const char *fname, const char *tname = "jpsi"){

   TH2F *h2D = new TH2F("h2D","h2D",hRcM->GetNbinsX(),hRcM->GetXaxis()->GetXmin(),hRcM->GetXaxis()->GetXmax(),hRcPt->GetNbinsX(),hRcPt->GetXaxis()->GetXmin(),hRcPt->GetXaxis()->GetXmax());
   smearPtFromTree(h2D,fdPtSig,momShape,fname,tname);
   TH1D *hM = (TH1D*)h2D->ProjectionX("hM");
   hRcM->Reset();
   hRcM->Add(hM);
   //hRcM->Reset();
   TH1D *hPt = (TH1D*)h2D->ProjectionY()->Clone("hPt");
   hRcPt->Reset();
   h2D->Draw("col");
   hRcPt->Add(hPt);
   //hPt->Reset();
}

//------------------------------------------------------
void compareInput(int flag = 1){ //1 is toy MC, 2 is embedding

   ifstream inf("txt/reso_fit.txt");
   double reso[2],resoErr[2];
   for(int i=0;i<2;i++){
      inf>>reso[i]>>resoErr[i];
   }
   inf.close();

   TF1 *fReso = new TF1("fReso","sqrt([0]*[0]*x*x+[1]*[1])",0,15);
   fReso->SetParameters(reso);
   fReso->SetNpx(1000);

   ifstream inf1("txt/dpt_CBfit.txt");
   double pars[7],parsErr[7];
   for(int i=0;i<7;i++){
      inf1>>pars[i]>>parsErr[i];
   }
   inf1.close();

   TF1 *fmomShape = new TF1("momShape", CrystalBall2, -1., 1., 7);
   fmomShape->SetParameters(pars);
   fmomShape->SetNpx(1000);

   TFile *fdata = new TFile("../ana/out/jpsi_mass_HT2_2TrkPid.root");
   TH1D *hMAll = (TH1D*)fdata->Get("hJpsiMAll");

   TFile *fsim = new TFile("OutFile_wTree_cent_0_9_resoNorm.root");   
   //TH2F *rcJpsiMPt = (TH2F*)fsim->Get("rcJpsiMPt_1"); // pt weighted
   TH2F *rcJpsiMPt = (TH2F*)fsim->Get("rcJpsiMPt"); // not weighted
   rcJpsiMPt->Sumw2();
   TH1D *rcM = (TH1D*)rcJpsiMPt->ProjectionX("rcM"); // 
   TH1D *rcPt = (TH1D*)rcJpsiMPt->ProjectionY("rcPt"); // 

   //TF1 *funTPCEffPt = new TF1("funTPCEffPt","[0]*exp(-pow([1]/x,[2]))+[3]*exp(-pow((x-[4])/[5],2))",0.,20);
   //funTPCEffPt->SetParameters(0.7503,0.1513,1.122,0.04132,0.3967,0.377);
   //TH1D *simRcM = (TH1D*)rcM->Clone("simRcM");
   TH1D *simRcM = new TH1D("simRcM","",400,2.,4.);
   simRcM->Reset();
   TH1D *simRcPt = (TH1D*)rcPt->Clone("simRcPt");
   simRcPt->Reset();

   if(flag==1){
      smearPt(1000000, rcPt, simRcM, simRcPt, fReso, fmomShape);
   } else if(flag==2){
      cout<<"smearing pt for tree"<<endl;
      smearPtFromTree(simRcM, simRcPt, fReso, fmomShape, "OutFile_wTree_cent_0_9_resoNorm.root", "jpsi");
   }

   TCanvas *c1 = new TCanvas("c1","c1",1000,400);
   c1->Divide(2,1);
   c1->cd(1);
   gPad->SetLogy();

   //rcPt->SetMarkerStyle(20);
   //rcPt->SetMarkerColor(1);
   //rcPt->SetLineColor(1);
   //rcPt->Draw("pe1");
   //simRcPt->SetLineColor(2);
   //simRcPt->Draw("histsame");

   float xmin = 2.8, xmax = 3.4;

   rcM->SetAxisRange(xmin,xmax,"X");
   simRcM->SetAxisRange(xmin,xmax,"X");
   hMAll->SetAxisRange(xmin,xmax,"X");

   rcM->Scale(1./rcM->Integral(rcM->FindBin(xmin),rcM->FindBin(xmax)-1));
   simRcM->Scale(1./simRcM->Integral(simRcM->FindBin(xmin),simRcM->FindBin(xmax)-1)*rcM->GetBinWidth(1)/simRcM->GetBinWidth(1));

   hMAll->Scale(1./hMAll->Integral(hMAll->FindBin(xmin),hMAll->FindBin(xmax)-1)*rcM->GetBinWidth(1)/hMAll->GetBinWidth(1));
   rcM->SetMarkerStyle(20);
   rcM->SetMarkerColor(1);
   rcM->SetLineColor(1);
   rcM->Draw("pe1");
   hMAll->SetMarkerStyle(21);
   hMAll->SetMarkerColor(4);
   hMAll->SetLineColor(4);
   hMAll->Draw("pe1same");
   simRcM->SetLineColor(2);
   simRcM->SetLineWidth(3);
   simRcM->Draw("histsame");
   TLegend *leg = new TLegend(0.25,0.2,0.5,0.4);
   leg->AddEntry(rcM,"Embedding","lp");
   leg->AddEntry(simRcM,"Simulation","l");
   leg->AddEntry(hMAll,"Data","p");
   leg->Draw();

   c1->cd(2);
   gPad->SetLogy(0);
   rcM->Draw("pe1");
   hMAll->Draw("pe1same");
   simRcM->Draw("histsame");
   c1->Print(Form("fig/compareJpsiM_mc_%d.pdf",flag));
}

Double_t funMass(double *x, double *par){

   if(!hMcMass) return 0;
   int   ptbin   = hMcMass->FindBin(x[0]);
   if(x[0]<hMcMass->GetBinCenter(ptbin)&&ptbin>1) ptbin = ptbin-1;
   double ptlw   = hMcMass->GetBinCenter(ptbin);
   double ptup   = hMcMass->GetBinCenter(ptbin+1);

   double a1 = hMcMass->GetBinContent(ptbin);
   double a2 = hMcMass->GetBinContent(ptbin+1);
   double p0 = (a1-a2)/(ptlw-ptup);
   double p1 = (a2*ptlw-a1*ptup)/(ptlw-ptup);

   return par[0]*(p0*x[0]+p1);
}

//------------------------------------------------------
//void compareData(){
//   TFile *fdata = new TFile("../ana/out/jpsi_mass_HT2_2TrkPid.root");
//   TH1D *hMAll = (TH1D*)fdata->Get("hJpsiMAll");
//   TH1D *hm[nJpsiPt];
//   for(int i=0;i<nJpsiPt;i++){
//      hm[i] = (TH1D*)fdata->Get("hJpsiMAll");
//   }
//
//   TFile *fth = new TFile("../ana/out/jpsi_xsec_final.root");
//   TGraphErrors *gCemZY = (TGraphErrors*)fth->Get("gCemZY");
//}

//------------------------------------------------------
void fitDataM(TH1D *hSim, TH1D *hData, Double_t &norm, Double_t &chi2, Int_t &ndf){
   hMcMass = (TH1D*)hSim->Clone("hMcMass");
   TF1 *fM = new TF1("fM",funMass,2.5,3.3,1);
   fM->SetParameter(0,1.);

   hData->Fit(fM,"R");
   norm = fM->GetParameter(0);
   chi2 = fM->GetChisquare();
   ndf = fM->GetNDF();
   hMcMass->Reset();
}

//------------------------------------------------------
void getBestResoPars(float amin = 0., float amax = 0.004, const int nIters = 20, int flag=1){
   double step = (amax-amin)/nIters;
   if(step<=0){
      cout<<"wrong inputs!"<<endl;
      return;
   }
   TH2F *hSimRcMPt[nIters];
   TH1D *hSimRcM[nIters];
   Double_t a[nIters];
   ifstream inf("txt/reso_fit.txt");
   double reso[2],resoErr[2];
   for(int i=0;i<2;i++){
      inf>>reso[i]>>resoErr[i];
   }
   inf.close();
   Double_t b = reso[1];

   ifstream inf1("txt/dpt_CBfit.txt");
   double pars[7],parsErr[7];
   for(int i=0;i<7;i++){
      inf1>>pars[i]>>parsErr[i];
   }
   inf1.close();

   TF1 *fmomShape = new TF1("momShape", CrystalBall2, -1., 1., 7);
   fmomShape->SetParameters(pars);
   fmomShape->SetNpx(1000);

   TFile *fdata = new TFile("../ana/out/jpsi_mass_HT2_2TrkPid.root");
   TH1D *hMAll = (TH1D*)fdata->Get("hJpsiMAll");

   TFile *fsim = new TFile("OutFile_wTree_cent_0_9_resoNorm.root");   
   TH2F *rcJpsiMPt = (TH2F*)fsim->Get("rcJpsiMPt");
   rcJpsiMPt->Sumw2();
   TH1D *rcM = (TH1D*)rcJpsiMPt->ProjectionX("rcM"); // flat pT
   TH1D *rcPt = (TH1D*)rcJpsiMPt->ProjectionY("rcPt"); // 

   TF1 *fReso = new TF1("fReso","sqrt([0]*[0]*x*x+[1]*[1])",0,20);
   fReso->SetNpx(1000);
   Double_t chi2[nIters];
   Int_t ndf[nIters];
   Double_t fitPar[nIters];
   Double_t chi2Min = 1e9;
   Double_t ndfMin = 1e9;
   Int_t bestFitIdx = -1;
   for(int i=0;i<nIters;++i){
      a[i] = amin + (i+1)*step;
      cout<<"generating rc mass for a = "<<a[i]<<" b = "<<b<<endl;
      fReso->SetParameters(a[i],b);

      if(flag==1){ //toy MC
         hSimRcMPt[i] = (TH2F*)rcJpsiMPt->Clone(Form("hSimRcMPt_%d",i));
         hSimRcMPt[i]->Reset();
         smearPt(100000,rcPt,hSimRcMPt[i],fReso,fmomShape);
      }
      if(flag==2){
         hSimRcMPt[i] = (TH2F*)rcJpsiMPt->Clone(Form("hSimRcMPt_%d",i));
         hSimRcMPt[i]->Reset();
         smearPtFromTree(hSimRcMPt[i],fReso,fmomShape,"OutFile_wTree_cent_0_9_resoNorm.root","jpsi");
      }
      hSimRcM[i] = (TH1D*)hSimRcMPt[i]->ProjectionX(Form("hSimRcM_%d",i));
   
      chi2[i] = 1e10;
      fitDataM(hSimRcM[i],hMAll,fitPar[i],chi2[i],ndf[i]);
      if(chi2[i]<chi2Min){ 
         chi2Min = chi2[i]; 
         ndfMin = ndf[i];
         bestFitIdx = i;
      }
   }

   TGraphErrors *gFit = new TGraphErrors(nIters,a,chi2,0,0);
   gFit->SetName("gFit");
   gFit->SetMarkerStyle(20);
   gFit->SetMarkerColor(2);
   gFit->GetXaxis()->SetTitle("a");
   gFit->GetYaxis()->SetTitle("#chi^{2}");
   TCanvas *c1 = new TCanvas("c1","c1",1200,600);
   c1->Divide(2,1);
   c1->cd(1);
   gFit->Draw("ap");
   TF1 *fg = new TF1("fg","pol4",a[1],a[nIters-2]);
   fg->SetParameters(700,-500,1e8,1,1);
   gFit->Fit(fg);
   Double_t chi2Best = fg->GetMinimum();
   Double_t aBest = fg->GetMinimumX();

   Double_t dchi2 = 0;
   Double_t chi2Test = chi2Best;
   Double_t aTest = aBest;
   while(dchi2<1){
      aTest += 1e-4;
      chi2Test = fg->Eval(aTest);
      dchi2 = fabs(chi2Test-chi2Best);
   }
   cout<<"aTest = "<<aTest<<endl;
   Double_t aBestErr = fabs(aTest - aBest);
   cout<<"best a = "<<aBest<<" +/- "<<aTest-aBest<<endl;

   c1->cd(2);
   hMAll->Draw();
   TH1D *hBest;
   if(bestFitIdx>=0){
      hBest = (TH1D*)hSimRcM[bestFitIdx]->Clone("hBest");
      //hBest->Scale(fitPar[bestFitIdx]*hBest->GetBinWidth(1)/hMAll->GetBinWidth(1));
      hBest->Scale(fitPar[bestFitIdx]);
      hBest->Draw("histsame");
   }
   TLatex *tt1 = new TLatex();
   tt1->SetNDC(1);
   tt1->DrawLatex(0.22,0.88,Form("Best fit:"));
   cout<<"bestFitIdx = "<<bestFitIdx<<endl;
   if(bestFitIdx>=0){

      tt1->DrawLatex(0.22,0.88-0.05,Form("a = %f",a[bestFitIdx]));
      tt1->DrawLatex(0.22,0.88-0.10,Form("b = %f",b));
      tt1->DrawLatex(0.22,0.88-0.15,Form("#chi^{2}/ndf = %3.1f/%d",chi2[bestFitIdx],ndf[bestFitIdx]));
   }
   c1->Print(Form("fig/bestFit_flag_%d.pdf",flag));

   ofstream outf(Form("txt/bestFit_flag_%d.txt",flag));
   outf<<"aTest = "<<aTest<<endl;
   outf<<"best a = "<<aBest<<" +/- "<<aTest-aBest<<endl;
   outf<<"b = "<<b<<endl;
   outf<<"bestFitIdx = "<<bestFitIdx<<endl;
   if(bestFitIdx>=0){
      outf<<"chi2/ndf = "<<chi2[bestFitIdx]<<"/"<<ndf[bestFitIdx]<<endl;
   }
   outf.close();

   TFile *fout = new TFile(Form("out/bestFit_flag_%d.root",flag),"recreate");
   
   for(int i=0;i<nIters;i++){
      hSimRcMPt[i]->Write();
      hSimRcM[i]->Write();
   }
   gFit->Write();
   fout->Close();
   
}

//------------------------------------------------------
void doSmear(int mode = 5){
   TTimer *timer = new TTimer();
   mRan = new TRandom3();
   unsigned long long tmp = (unsigned long long)timer->GetAbsTime();
   unsigned int stime = tmp/gRandom->Rndm() ;
   cout << " Random seed is " << stime <<endl;
   mRan->SetSeed(stime);

   if(mode==1) getReso();
   if(mode==2) fitNormReso(2);
   if(mode==3) compareInput(1);
   if(mode==4) getBestResoPars(0.0022,0.0046,50,2);
   if(mode==5){
      getReso();
      fitNormReso(2);
      compareInput(2);
      return;
      getBestResoPars(0.0022,0.0046,50,2);
   }
}
