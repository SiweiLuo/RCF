// Example on how to use the new Minimizer class in ROOT
//  Show usage with all the possible minimizers. 
// Minimize the Rosenbrock function (a 2D -function)
// This example is described also in 
// http://root.cern.ch/drupal/content/numerical-minimization#multidim_minim
// input : minimizer name + algorithm name
// randomSeed: = <0 : fixed value: 0 random with seed 0; >0 random with given seed 
//
//Author: L. Moneta Dec 2010

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

#define xNBins 100
#define yNBins 100

double RosenBrock(const double *xx )
{
  const Double_t p1 = xx[0];//lambda_theta
  const Double_t p2 = xx[1];//lambda_phi

	TF2* sigma = new TF2("sigma","1+[0]*y*y+[1]*(1-y*y)*TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi(),-1,1);	
	sigma->SetParameters(p1,p2);

	TRandom3 rndeff;
	TH2F* effhist = (TH2F*)eff2D[ifile][itrig][ipt][iframe][2]->Clone();
	effhist->SetName(Form("eff_2dhist_%d_%d_%d_%d",ifile,itrig,ipt,iframe));
	for(int i=0;i<10000;i++){
		sigma->GetRandom2(phi,theta);
		double eff = effhist->GetBinContent(effhist->GetXaxis()->FindBin(phi),effhist->GetYaxis()->FindBin(theta));
		if(rndeff.Uniform()<eff) random2->Fill(phi,theta);	
	}
	random2->Scale(1./random2->Integral());

	Int_t THETANBIN=10,PHINBIN=10;
	double chisquare;
	for(int thetabin=1;thetabin<THETANBIN+1;thetabin++){
		for(int phibin=1;phibin<PHINBIN+1;phibin++){
			if(sample2->GetBinContent(thetabin,phibin)>=1e-3) chisquare += TMath::Power((random2->GetBinContent(thetabin,phibin)-sample2->GetBinContent(thetabin,phibin))/(sample2->GetBinError(thetabin,phibin)),2);	
		}
	}

  return chisquare;
}
 
int NumericalMinimization(const char * minName = "Minuit2",
                          const char *algoName = "" , 
                          int randomSeed = -1)
{
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are: 
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   ROOT::Math::Minimizer* min = 
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
   min->SetMaxIterations(10000);  // for GSL 
   min->SetTolerance(0.001);
   min->SetPrintLevel(1);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type 
   ROOT::Math::Functor f(&RosenBrock,2); 
   double step[2] = {0.01,0.01};
   // starting point
    
   double variable[2] = { -1.,1.2};
   if (randomSeed >= 0) { 
      TRandom2 r(randomSeed);
      variable[0] = r.Uniform(-1,1);
      variable[1] = r.Uniform(-1,1);
   }
 
   min->SetFunction(f);
 
   // Set the free variables to be minimized!
   min->SetVariable(0,"x",variable[0], step[0]);
   min->SetVariable(1,"y",variable[1], step[1]);
 
   // do the minimization
   min->Minimize(); 
 
   const double *xs = min->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
             << min->MinValue()  << std::endl;

   // expected minimum is 0
   if ( min->MinValue()  < 1.E-4  && f(xs) < 1.E-4) 
      std::cout << "Minimizer " << minName << " - " << algoName 
                << "   converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer " << minName << " - " << algoName 
                << "   failed to converge !!!" << std::endl;
      Error("NumericalMinimization","fail to converge");
   }
 
   return 0;
}

void correcteddata(int sys2 = 0){
	int selectedfile,file;
	if(sys2>=0 && sys2<=NFILE) file=sys2,selectedfile=sys2+1;
	else return;
	
	for(;file<selectedfile;file++){
		efficiencyfile[file] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/OutFile_cent_0_9_%d.root",file),"read"); 
		for(int trig=0;trig<NTRIG;trig++){//marker
			if((sys2>=20 && sys2<=23) || (sys2>=29 && sys2<=30) || sys2==35) rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,sys2),"read");//marker
			else rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160210_Barbara/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,0),"read");//marker
			for(int pt=0;pt<NPT;pt++){
				for(int frame=0;frame<NFRAME;frame++){
					canvas[trig][pt][frame] = new TCanvas(Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),Form("%d_%s_%d_%d_frame%d",file,trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1],frame),3200,1600);
					canvas[trig][pt][frame]->Divide(4,3);
					canvas[trig][pt][frame]->cd(1);
					if(frame==0){
						rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPt"));
						rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtBG"));
						rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_hx");
						rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_hx");
					}
					else{
						rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][0]->Get("hJpsiCosThetaPhiPtCS"));
						rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][0]->Get("hJpsiCosThetaPhiPtCSBG"));
						rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_cs");
						rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_cs");
					}
					rawdata3D[file][trig][pt][frame][0]->Draw();

					canvas[trig][pt][frame]->cd(2);
					rawdata3D[file][trig][pt][frame][1]->Draw();

					canvas[trig][pt][frame]->cd(3);
					rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike-like");
					if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_hx");
					else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike-like_cs");
					rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);
					rawdata3D[file][trig][pt][frame][2]->Draw();

					int max = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
					int min = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	

					if(frame==0){
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPt1");
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
					}
					else{	
						eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
						eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
						eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPtCS1");
						eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
					}

					rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);
					rawdata3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					rawdata3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
					eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));

					canvas[trig][pt][frame]->cd(4);
					rawdata2D[file][trig][pt][frame][0] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
					rawdata2D[file][trig][pt][frame][0]->SetName(Form("rawdata_%s_pT%d_%d",trigName[trig].Data(),PtEdge[pt],PtEdge[pt+1]));
					rawdata2D[file][trig][pt][frame][0]->RebinX(NREBIN)->RebinY(NREBIN);	
					rawdata2D[file][trig][pt][frame][0]->Draw("colz");

					canvas[trig][pt][frame]->cd(5);
					eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
					eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);

					eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
					eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);

					eff2D[file][trig][pt][frame][2] = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
					eff2D[file][trig][pt][frame][2]->SetName(Form("eff_ratio_%d_%d_%d_%d",file,trig,pt,frame));
					eff2D[file][trig][pt][frame][2]->Divide(eff2D[file][trig][pt][frame][1]);
					eff2D[file][trig][pt][frame][2]->Draw("colz");

					if(pt==0) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf(",trigName[trig].Data(),frame));			
					else if(pt==NPT-1) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf)",trigName[trig].Data(),frame));			
					else canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf",trigName[trig].Data(),frame));			
				}
			}
		}
	}
}


