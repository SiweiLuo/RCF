#include <iostream> 
#include "TCanvas.h"
#include "TGaxis.h"

TString variable[4] = {"#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}","#lambda_{inv}"};
TString component[8] = {"default","p_{T} smearing","Dedx","Dca","NHitsFit","Dsmadc","1/beta","nsigmae"};
TString FRM[2]={"HX","CS"};

void uncertainty_pt(){
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);	

	Double_t uncertainty[8][4][4][2];// component, pt, variable
	TLegend *uncertaintyleg[4][2];
	ifstream file("sPlot3D_uncertainty_com.txt");	

	int pt=0,var=0,frame=0;
	/*
	   while(file>>uncertainty[0][pt][var][frame]>>uncertainty[1][pt][var][frame]>>uncertainty[2][pt][var][frame]>>uncertainty[3][pt][var][frame]>>uncertainty[4][pt][var][frame]>>uncertainty[5][pt][var][frame]>>uncertainty[6][pt][var][frame]>>uncertainty[7][pt][var][frame]){
	   if(var%3==0 && var!=0) var=-1,frame++;
	   if(pt%3==0 && pt!=0) pt=-1,var++;
	   pt++;
	   cout<<" pt ==== "<<pt<<"  var ==== "<<var<<"   frame == "<<frame<<endl;
	   }	
	   */

	TFile *errfile[2][6][2];
	Double_t ydef[4][4][2];

	for(int pt=0,trig=0,ipt=1,itrig=1;pt<4;pt++){
		for(int frame=0;frame<2;frame++){
			for(int var=0;var<4;var++)
			{
				file>>uncertainty[0][pt][var][frame]>>uncertainty[1][pt][var][frame]>>uncertainty[2][pt][var][frame]>>uncertainty[3][pt][var][frame]>>uncertainty[4][pt][var][frame]>>uncertainty[5][pt][var][frame]>>uncertainty[6][pt][var][frame]>>uncertainty[7][pt][var][frame];
				ipt=pt+1;
				
				if(pt>1)trig=1,itrig=3;
				errfile[trig][pt][frame] = new TFile(Form("~/polresults/20160707/functional/splot_3D_20171005_2eid_inv30/splot_3D_%d_%d_%d_0_20171005.root",itrig,ipt,frame),"read");
				if(var!=3)ydef[pt][var][frame] = ((TH3F*)errfile[trig][pt][frame]->Get("hlambda"))->GetMean(var+1);
				else ydef[pt][var][frame] = ((TH1F*)errfile[trig][pt][frame]->Get("hlambda_inv"))->GetMean();
			}
		}
	}

	TGraphErrors *uncertaintyvspt[9][4][4][2];//component, variable, pt,frame
	TGraphErrors *sysvspt[4][4][2];
	//	Double_t y[4],yerr[4];

	//	Double_t x[4] = {2.5,3.5,5,7};	
	//	Double_t xerr[4] = {0.5,0.5,1.,1.};
	//	Double_t ydef[4] = {0.1,0.1,0.1,0.1};// read ydef from file;

	TCanvas *c = new TCanvas("uncertaintyc","uncertaintyc",4000,2000);
	//c->SetGrid(false);
	c->Divide(4,2);

	TCanvas *sysc = new TCanvas("sysc","sysc",4000,2000);
	sysc->Divide(4,2);

	Double_t x,xerr,y,yerr;
	Double_t RMS[4][4][2],rms;

	for(int ifrm=0;ifrm<2;ifrm++){
		for(int ivar=0;ivar<4;ivar++){
			c->cd(ivar+4*ifrm+1);
			uncertaintyleg[ivar][ifrm] = new TLegend(0.6,0.6,0.9,0.9);

			for(int ipt=0;ipt<4;ipt++){
				RMS[ipt][ivar][ifrm]=0.;
				for(int com=0;com<8;com++){
					cout<<" com === "<<com<<" ivar = "<<ivar<<" ipt === "<<ipt<<" ifrm == "<<ifrm<<endl;
					cout<<"uncertainty ==== "<<uncertainty[com][ipt][ivar][ifrm]<<endl;

					if(ipt==0)x = 2.5,xerr=0.5; 
					else if(ipt==1) x=3.5,xerr=0.5,;
					else if(ipt==2) x=5.,xerr=1.;
					else if(ipt==3) x=7.,xerr=1.;
					y = ydef[ipt][ivar][ifrm]+uncertainty[com][ipt][ivar][ifrm];
					RMS[ipt][ivar][ifrm] += uncertainty[com][ipt][ivar][ifrm]*uncertainty[com][ipt][ivar][ifrm];
					uncertaintyvspt[com][ivar][ipt][ifrm]=new TGraphErrors(1,&x,&y,&xerr,0);
					uncertaintyvspt[com][ivar][ipt][ifrm]->SetMarkerColor(com+2);
					uncertaintyvspt[com][ivar][ipt][ifrm]->SetLineColor(com+2);
					if(com==0 && ipt==0){
						uncertaintyvspt[com][ivar][ipt][ifrm]->Draw("ap");
						uncertaintyvspt[com][ivar][ipt][ifrm]->GetYaxis()->SetRangeUser(-2,2);
						uncertaintyvspt[com][ivar][ipt][ifrm]->GetXaxis()->SetLimits(1,9);
						uncertaintyvspt[com][ivar][ipt][ifrm]->SetTitle(Form("lambda parameters in %s ;p_{T} GeV/c;%s",FRM[ifrm].Data(),variable[ivar].Data()));
					}
					else uncertaintyvspt[com][ivar][ipt][ifrm]->Draw("psame");
					if(ipt==0)uncertaintyleg[ivar][ifrm]->AddEntry(uncertaintyvspt[com][ivar][ipt][ifrm],Form("%s",component[com].Data()),"lep");
				}
				uncertaintyleg[ivar][ifrm]->Draw("same");

				RMS[ipt][ivar][ifrm] = RMS[ipt][ivar][ifrm]/7.;
				RMS[ipt][ivar][ifrm] = TMath::Sqrt(RMS[ipt][ivar][ifrm]);
				rms = RMS[ipt][ivar][ifrm];
				def = ydef[ipt][ivar][ifrm];
				sysvspt[ipt][ivar][ifrm] = new TGraphErrors(1,&x,&def,&xerr,&rms);
			}
		}
	}

	for(int ifrm=0;ifrm<2;ifrm++){
		for(int ivar=0;ivar<4;ivar++){
			sysc->cd(ivar+4*ifrm+1);		
			for(int ipt=0;ipt<4;ipt++){
				sysvspt[ipt][ivar][ifrm]->SetFillStyle(1001);
				sysvspt[ipt][ivar][ifrm]->SetFillColorAlpha(kGray,0.35);
				if(ipt==0){
					sysvspt[ipt][ivar][ifrm]->GetYaxis()->SetRangeUser(-2,2);
					sysvspt[ipt][ivar][ifrm]->GetXaxis()->SetLimits(1,9);
					sysvspt[ipt][ivar][ifrm]->SetTitle(Form("lambda parameters in %s;p_{T} GeV/c;%s",FRM[ifrm].Data(),variable[ivar].Data()));
					sysvspt[ipt][ivar][ifrm]->Draw("ape2");
				}
				else sysvspt[ipt][ivar][ifrm]->Draw("e2 psame");
			}	
		}
	}	
	c->SaveAs("~/WWW/Systematic/uncertainty_pt.pdf");
	sysc->SaveAs("~/WWW/Systematic/uncertainty_pt_syspt.pdf");	
}
