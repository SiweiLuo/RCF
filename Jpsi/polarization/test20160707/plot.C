#include "TFile.h"
#include "TCanvas.h"
#include "TPDF.h"

#define NTRIG 3
#define NPT 6
#define NFRAME 2 
#define NPLOT 6
#define NFILE 40
#define NREBIN 4
//#define rebin 0 // 1

TFile* rootfile;
TFile* combinedfile;
TH2F* rawdata;
TH2F* eff;
TH2F* chi2;
TH2F* truth;
TH2F* bestfit;
TH2F* check1Dfit;

TH2F* template_theta;
TH2F* template_phi;
TH2F* template_theta_eff;
TH2F* template_phi_eff;

TH2F* truthcomparison;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TLegend* leg;

TFile *efficiencyfile[NFILE];
TH1F *eff1D[NFILE][NTRIG][NPT][NFRAME][3]; //pass  , total , ratio;
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][2]; // pass , total ;

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

void plot(int file=0,int rebin=1){
	gStyle->SetPadRightMargin(0.2);	

	TFile* outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	
	//	TFile* templates = new TFile("histograms20161107.root","read");
	TFile* templates = new TFile("histograms20161206.root","read");
	efficiencyfile[file] = new TFile(Form("rootfile/OutFile_sys%d.root",file),"read");					

	for(int trig=0;trig<3;trig++){
		for(int frame=0;frame<2;frame++){
			comparisonlambda[trig][3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
			comparisonlambda[trig][1+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));
			//				if(phase==2) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_inv_%d_1eID",trigName[trig].Data(),frame));
			for(int pt=0;pt<6;pt++){
				TCanvas* c1 = new TCanvas(Form("trig%d_frame%d_pt%d",trig,frame,pt),Form("trig%d_frame%d_pt%d",trig,frame,pt),2000,2000);
				c1->Divide(4,4);
				combinedfile = new TFile(Form("~/polresults/20160707/rootcombined/lambda_file%d_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");		
				rootfile = new TFile(Form("~/polresults/20160707/rootfiles/chi2histograms_%d/chi2histograms_%d_%d_%d_%d_1_1.root",file,file,trig,pt,frame),"read");
				chi2 = (TH2F*)combinedfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));
				int xx,yy,zz;
				chi2->GetMinimumBin(xx,yy,zz);
				bestfit = (TH2F*)combinedfile->Get(Form("template_file%d_trg%d_pt%d_frame%d_x%d_y%d",file,trig,pt,frame,xx,yy));
				//				truth = (TH2F*)combinedfile->Get(Form("theta_%d_phi_%d",xx-1,yy-1));
				truth = (TH2F*)combinedfile->Get(Form("template_%d_%d_%d_%d_theta_%d_phi_%d",file,trig,pt,frame,xx-1,yy-1));
				rawdata = (TH2F*)rootfile->Get(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));//marker add file
				rawdata->Sumw2();
				eff = (TH2F*)rootfile->Get(Form("eff_%d_%d_%d_%d",file,trig,pt,frame));	
				eff->Sumw2();

				Double_t ptx,lambdatheta,lambdaphi;

				cout<<"comparisonlambda "<<comparisonlambda[trig][3*frame]<<endl;

				comparisonlambda[trig][3*frame]->GetPoint(pt,ptx,lambdatheta);
				comparisonlambda[trig][1+3*frame]->GetPoint(pt,ptx,lambdaphi);

				cout<<"lambdatheta="<<lambdatheta<<"lambdaphi="<<lambdaphi<<endl;

				check1Dfit = (TH2F*)templates->Get(Form("theta_%d_phi_%d",chi2->GetXaxis()->FindBin(lambdatheta)-1,chi2->GetYaxis()->FindBin(lambdaphi)-1));
				if(check1Dfit!=0x0)	{
					if(rebin==1){
						check1Dfit->RebinX(4);
						check1Dfit->RebinY(4);
					}
					check1Dfit->Multiply(eff);
				}
				//					cout<<"lambdatheta="<<lambdatheta<<"lambdaphi"<<lambdaphi<<endl;

				if(chi2==0x0 || bestfit==0x0 || truth==0x0 || rawdata==0x0 || eff==0x0) {
					if(chi2==0x0) 		cout<<"chi2 error"<<endl;
					if(bestfit==0x0) 	cout<<"bestfit error"<<endl;
					if(truth==0x0) 		cout<<"templates error"<<endl;
					if(rawdata==0x0) 	cout<<"rawdata error"<<endl;
					if(eff==0x0) 		cout<<"efficiency error"<<endl;
					continue;
				}

				if(check1Dfit!=0x0) check1Dfit->Scale(rawdata->Integral()/check1Dfit->Integral());
				c1->cd(1);
				truth->Draw("colz");
				c1->cd(2);
				eff->Draw("colz");
				c1->cd(3);
				bestfit->Draw("colz");
				c1->cd(4);
				rawdata->Draw("colz");	

				c1->cd(5);
				template_theta = (TH2F*)templates->Get(Form("theta_%d_phi_%d",chi2->GetXaxis()->FindBin(lambdatheta)-1,50));
				if(template_theta!=0x0){
					template_theta->SetName(Form("theta_%d_%d_%d_%d",file,trig,pt,frame));
					template_theta->SetTitle(Form("template #lambda_{#theta}=%.2f #lambda_{#phi}=0",lambdatheta));
					if(rebin==1){
						template_theta->RebinX(4);
						template_theta->RebinY(4);
					}
					//					template_theta->Multiply(eff);
					template_theta->Draw("colz");
				}
				c1->cd(6);
				if(template_theta!=0x0){
					template_theta_eff = (TH2F*)template_theta->Clone("template_theta_eff");
					template_theta_eff->SetTitle(Form("template*eff #lambda_{#theta}=%.2f #lambda_{#phi}=0",lambdatheta));
					template_theta_eff->Multiply(eff);
					template_theta_eff->Draw("colz");
				}

				c1->cd(7);
				template_phi = (TH2F*)templates->Get(Form("theta_%d_phi_%d",50,chi2->GetYaxis()->FindBin(lambdaphi)-1));
				if(template_phi!=0x0){
					template_phi->SetName(Form("phi_%d_%d_%d_%d",file,trig,pt,frame));
					template_phi->SetTitle(Form("template #lambda_{#theta}=0 #lambda_{#phi}=%.2f",lambdaphi));
					if(rebin==1){
						template_phi->RebinX(4);
						template_phi->RebinY(4);
					}
					//					template_phi->Multiply(eff);
					template_phi->Draw("colz");
				}
				c1->cd(8);
				if(template_phi!=0x0){
					template_phi_eff = (TH2F*)template_phi->Clone("template_phi_eff");
					template_phi_eff->SetTitle(Form("template*eff #lambda_{#theta}=0 #lambda_{#phi}=%.2f",lambdaphi));
					template_phi_eff->Multiply(eff);
					template_phi_eff->Draw("colz");
				}

				TH1F* template_theta_ratio;
				TH1F* template_phi_ratio;

				c1->cd(9);
				leg = new TLegend(0.1,0.7,0.4,0.9);
				rawdata->ProjectionY("rawdata_py");
				rawdata_py->SetMarkerStyle(21);	
				rawdata_py->SetMarkerColor(kRed);
				rawdata_py->SetLineColor(kRed);
				rawdata_py->SetMarkerSize(1.2);
				double maxval1 = rawdata_py->GetMaximum();
				rawdata_py->SetMaximum(1.7*maxval1);
				rawdata_py->Draw("");
				bestfit->ProjectionY("bestfit_py");
				bestfit_py->SetMarkerStyle(22);
				bestfit_py->SetMarkerColor(kBlue);
				bestfit_py->SetLineColor(kBlue);
				bestfit_py->Draw("same");
				if(check1Dfit!=0x0){
					check1Dfit->ProjectionY("check1Dfit_py")->Draw("same");
					check1Dfit_py->SetMarkerStyle(24);

					if(template_theta_eff!=0x0 && template_theta_eff->ProjectionY("template_theta_eff_py")){
						cout<<"rawdata py ====="<<rawdata_py<<endl;
						//					if(template_theta_eff_py!=0x0){
						template_theta_eff_py->Scale(rawdata->Integral()/template_theta_eff->Integral());
						template_theta_eff_py->SetMarkerStyle(23);
						template_theta_eff_py->SetMarkerColor(kGreen);
						template_theta_eff_py->SetLineColor(kGreen);
						template_theta_eff_py->Draw("same");
						template_theta_ratio = (TH1F*)template_theta_eff_py->Clone("template_theta_ratio");
					}
					}
					leg->AddEntry(rawdata_py,"raw data","lpe");
					leg->AddEntry(bestfit_py,"2D fit","lpe");
					if(check1Dfit!=0x0)	{
						leg->AddEntry(check1Dfit_py,"1D fit","lpe");
						leg->AddEntry(template_theta_eff_py,"#lambda_{#phi} = 0","lpe");
					}
					leg->Draw("same");

					c1->cd(10);
					leg = new TLegend(0.1,0.7,0.4,0.9);
					rawdata->ProjectionX("rawdata_px");
					rawdata_px->SetMarkerStyle(21);
					rawdata_px->SetMarkerColor(kRed);
					rawdata_px->SetLineColor(kRed);
					rawdata_px->SetMarkerSize(1.2);
					rawdata_px->SetFillStyle(3003);
					double maxval2 = rawdata_px->GetMaximum();
					rawdata_px->SetMaximum(1.7*maxval2);
					rawdata_px->Draw("");
					bestfit->ProjectionX("bestfit_px");
					bestfit_px->SetMarkerStyle(22);
					bestfit_px->SetMarkerColor(kBlue);
					bestfit_px->SetLineColor(kBlue);
					bestfit_px->Draw("same");
					if(check1Dfit!=0x0){
						check1Dfit->ProjectionX("check1Dfit_px")->Draw("same");
						check1Dfit_px->SetMarkerStyle(24);
						template_phi_eff->ProjectionX("template_phi_eff_px");
						template_phi_eff_px->Scale(rawdata->Integral()/template_phi_eff->Integral());
						template_phi_eff_px->SetMarkerStyle(23);
						template_phi_eff_px->SetMarkerColor(kGreen);
						template_phi_eff_px->SetLineColor(kGreen);
						template_phi_eff_px->Draw("same");
						template_phi_ratio = (TH1F*)template_phi_eff_px->Clone("template_phi_ratio");
					}
					rawdata_px->Draw("E2 samep");
					leg->AddEntry(rawdata_px,"raw data","lpe");
					leg->AddEntry(bestfit_px,"2D fit","lpe");
					if(check1Dfit!=0x0){
						leg->AddEntry(check1Dfit_px,"1D fit","lpe");
						leg->AddEntry(template_phi_eff_px,"#lambda_{#theta} = 0","lpe");	
					}
					leg->Draw("same");

					c1->cd(11);
					//				TH1F* rawdata_py_ratio = (TH1F*)rawdata_py->Clone("raw ratio py");
					TH1F* rawdata_py_ratio = (TH1F*)bestfit_py->Clone("raw ratio py");
					rawdata_py_ratio->Divide(rawdata_py);	
					rawdata_py_ratio->SetMaximum(4);
					rawdata_py_ratio->SetMinimum(-0.5);
					rawdata_py_ratio->Draw();					
					if(template_theta_ratio!=0x0){
						template_theta_ratio->Divide(rawdata_py);
						template_theta_ratio->Draw("same");
					}
					if(check1Dfit!=0x0){
						TH1F* check1Dfit_py_ratio = (TH1F*)check1Dfit_py->Clone("check1dfit_py_ratio");
						check1Dfit_py_ratio->Divide(rawdata_py);
						check1Dfit_py_ratio->Draw("same");
					}			
					//				comparisonlambda[trig][3*frame]->Draw("ap");
					//				if(check1Dfit!=0x0)check1Dfit->Draw("colz");
					c1->cd(12);
					//				comparisonlambda[trig][1+3*frame]->Draw("ap");	
					//				TH1F* rawdata_px_ratio = (TH1F*)rawdata_px->Clone("raw ratio px");
					TH1F* rawdata_px_ratio = (TH1F*)bestfit_px->Clone("raw ratio px");
					rawdata_px_ratio->Divide(rawdata_px);	
					rawdata_px_ratio->SetMaximum(4);
					rawdata_px_ratio->SetMinimum(-0.5);
					rawdata_px_ratio->Draw();					
					if(template_phi_ratio!=0x0){
						template_phi_ratio->Divide(rawdata_px);
						template_phi_ratio->Draw("same");
					}
					if(check1Dfit!=0x0){
						TH1F* check1Dfit_px_ratio = (TH1F*)check1Dfit_px->Clone("check1dfit_px_ratio");
						check1Dfit_px_ratio->Divide(rawdata_px);
						check1Dfit_px_ratio->Draw("same");
					}			

					c1->cd(13);//raw data
					rawdata->ProjectionY("raw_py");
					raw_py->SetTitle("raw cos#theta");
					raw_py->Draw();


					c1->cd(14);
//					TH1F* eff_py = (TH1F*)eff->ProjectionY("eff_py");
//					eff_py->Draw();
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

					eff1D[file][trig][pt][frame][0] = (TH1F*)eff2D[file][trig][pt][frame][0]->ProjectionY("eff_pass_py");
					eff1D[file][trig][pt][frame][0]->SetName(Form("eff_1d_ratio_%d_%d_%d_%d",file,trig,pt,frame));
					eff1D[file][trig][pt][frame][0]->Divide((TH1F*)eff2D[file][trig][pt][frame][1]->ProjectionY("eff_total_py"));
					eff1D[file][trig][pt][frame][0]->SetTitle("cos#theta efficiency; cos#theta");	
					eff1D[file][trig][pt][frame][0]->Draw();	
//					eff2D[file][trig][pt][frame][0]->Draw();
//					eff3D[file][trig][pt][frame][0]->Draw();

					c1->cd(15);
					TF1* costheta = new TF1("costheta","[0]*(1+[1]*x*x)",-1,1);
					rawdata->ProjectionY("raw_py_divided");
					raw_py_divided->Divide(eff1D[file][trig][pt][frame][0]);
					raw_py_divided->SetTitle("corrected cos#theta");
					raw_py_divided->Draw();
					raw_py_divided->Fit("costheta");

					c1->cd(16);
					chi2->Draw("surf2");		
					/*
					   c1->cd(13);
					   TH2F* rawdata_ratio = (TH2F*)rawdata->Clone("raw ratio");
					   rawdata_ratio->Divide(bestfit);
					   rawdata_ratio->Draw("colz");
					 */

					c1->SaveAs(Form("~/polresults/20160707/pdf/sys_%d/trig%d_frame%d_pt%d.pdf",file,trig,frame,pt));

					if((TF2*)combinedfile->Get("mlefcn")==0x0 || (TF2*)combinedfile->Get("contfcn")==0x0 || (TGraph*)combinedfile->Get("small_cont_min_half")==0x0 || (TGraphErrors*)combinedfile->Get("lambda")==0x0 || (TH2F*)combinedfile->Get("difference")==0x0) continue;
					TCanvas* c2 = new TCanvas(Form("fitlambda_%d_%d_%d_%d",file,trig,pt,frame),Form("fitlambda_%d_%d_%d_%d",file,trig,pt,frame),3200,600);
					c2->Divide(4,1);
					c2->cd(1);
					chi2->Draw("surf2");
					c2->cd(2);	
					(TF2*)combinedfile->Get("mlefcn")->Draw("colz");
					(TGraph*)combinedfile->Get("small_cont_min_half")->Draw("same");
					(TGraphErrors*)combinedfile->Get("lambda")->Draw("same");
					c2->cd(3);
					(TF2*)combinedfile->Get("contfcn")->Draw("colz");
					(TGraph*)combinedfile->Get("small_cont_min_half")->Draw("same");
					(TGraphErrors*)combinedfile->Get("lambda")->Draw("same");
					c2->cd(4);
					(TH2F*)combinedfile->Get("difference")->Draw("colz");
					(TGraph*)combinedfile->Get("small_cont_min_half")->Draw("same");
					(TGraphErrors*)combinedfile->Get("lambda")->Draw("same");
					
					c2->SaveAs(Form("~/polresults/20160707/pdf/sys_%d/fitlambda_%d_%d_%d_%d.pdf",file,file,trig,pt,frame));
				}
			}
		}
	}
