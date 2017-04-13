#include <algorithm>

#define NFILE 40
#define NTRIG 3
#define NPT 6
#define NPHASE 2 // theta, phi
#define NFRAME 2
#define NPLOT 6
#define NREBIN 4

Float_t z1[40/NREBIN],z2[40/NREBIN],errorz1[40/NREBIN],errorz2[40/NREBIN],x[40/NREBIN],y[40/NREBIN];
Double_t chisquare,chisquare1,chisquare2;
Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht0","ht1","ht2"};

TFile* jpsifile;
TFile* crosscheckfile;
TFile* lambdas;

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME][2];// raw data, corrected data
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH1F *costheta_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to costheta pass , total 
TH1F *phi_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to phi pass , total 
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;

TH1F *costheta[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TH1F *phi[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TGraphAsymmErrors* efficiency1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];
TH1F * eff1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];

TGraphAsymmErrors *efficiency;
TFile *efficiencyfile[NFILE];
TFile *rawdatafile[NFILE][NTRIG];

TF1 *fit1[NPT][NFRAME][2]; //npt frame default/revised
TF1 *fit2[NPT][NFRAME][2];

double lambda_theta[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_theta_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_invariant[NFILE][NTRIG][NPT][NFRAME][3];// central value , statistic error , systematic error

double fit2dtheta[NFILE][NTRIG][NPT][NFRAME];
double fit2dtheta_err[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi_err[NFILE][NTRIG][NPT][NFRAME];
double fitparameters[NFILE][NTRIG][NPT][NFRAME][4];//theta , theta_err , phi , phi_err
double ptsmearing[NTRIG][NPT][NFRAME][NPHASE],weight[NTRIG][NPT][NFRAME][NPHASE],nhitsfit[NTRIG][NPT][NFRAME][NPHASE],nsigma[NTRIG][NPT][NFRAME][NPHASE],dsmadc[NTRIG][NPT][NFRAME][NPHASE],beta[NTRIG][NPT][NFRAME][NPHASE],poe[NTRIG][NPT][NFRAME][NPHASE],pol[NTRIG][NPT][NFRAME][NPHASE];
double systematic_error[NTRIG][NPT][NFRAME][NPHASE];

Double_t matrix[4][4];
Double_t covariant[NTRIG][NPT][NFRAME];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker
TCanvas *systematic[NTRIG];

double fParamVal[NFILE][NTRIG][NPT][NFRAME][4];
double fParamErr[NFILE][NTRIG][NPT][NFRAME][4];
double arglist[10];
Int_t ierflg=0;
Int_t fitflag[4][NPT][2][4];
TH1F* fitdata[NTRIG][NPT][NFRAME][NPHASE];// theta , phi

TGraphErrors* lambda_theta_hx;// marker write them as an array
TGraphErrors* lambda_phi_hx; 
TGraphErrors* lambda_theta_cs;
TGraphErrors* lambda_phi_cs; 
TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda_parameters_sys[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant

TLegend* legend;

TFile *outputfile;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

void do2Dcorrection(int uncertainty = 0){
	//	gStyle->SetPadRightMargin(0.2);
	outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	
	gStyle->SetOptStat(false);

	printscheme();
	systematic(uncertainty);
	lambdaparameters(uncertainty);
	//	cout<<"The calculation of J/psi polarization is done"<<endl;
}

Double_t angular(Double_t *x, Double_t *par){
	//"[3]*(1+[0]*x[0]*x[0]+[1]*TMath::Sin(TMath::ACos(x[0]))*TMath::Sin(TMath::ACos(x[0]))*TMath::Cos(2*x[1])+[2]*TMath::Sin(2*TMath::ACos(x[0]))*TMath::Cos(x[1]))"
	Double_t result = par[3]*(1+par[0]*x[1]*x[1]+par[1]*(1-x[1]*x[1])*TMath::Cos(2*x[0])+par[2]*TMath::Sin(2*TMath::ACos(x[1]))*TMath::Cos(x[0]));
	return result;
}

void correcteddata(int sys2 = 0){
	int selectedfile,file;
	if(sys2>=0 && sys2<=NFILE) file=sys2,selectedfile=sys2+1;
	else if(sys2==100) file=0,selectedfile=NFILE;
	else return;

	//	TF2 *f2 = new TF2("f2","[3]*(1+[0]*x[0]*x[0]+[1]*TMath::Sin(TMath::ACos(x[0]))*TMath::Sin(TMath::ACos(x[0]))*TMath::Cos(2*x[1])+[2]*TMath::Sin(2*TMath::ACos(x[0]))*TMath::Cos(x[1]))",-1,1,-TMath::Pi(),TMath::Pi(),4);

	TF2 *f2 = new TF2("f2",angular,-TMath::Pi(),TMath::Pi(),-1,1,4);
	TF1 *fx = new TF1("fx",func1,-1,1,2);
	TF1 *fy = new TF1("fy",func2,-TMath::Pi(),TMath::Pi(),4);

	for(;file<selectedfile;file++){
		efficiencyfile[file] = new TFile(Form("~/jpsi/test20160707/rootfile/OutFile_cent_0_9_%d.root",file),"read"); 
		for(int trig=0;trig<NTRIG;trig++){//marker
			if((sys2>=20 && sys2<=23) || (sys2>=29 && sys2<=30) || sys2==35) rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160707/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,sys2),"read");//marker
			else rawdatafile[file][trig] = new TFile(Form("~/jpsi/test20160707/rootfile/%s_trg%d_1TrkPid_0_%d.ana.root",trigSet[trig].Data(),trig+3,0),"read");//marker
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

					costheta[file][trig][pt][frame][0] = (TH1F*)rawdata2D[file][trig][pt][frame][0]->ProjectionY();		
					phi[file][trig][pt][frame][0] = (TH1F*)rawdata2D[file][trig][pt][frame][0]->ProjectionX();		

					eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
					eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
					//					eff2D[file][trig][pt][frame][0]->GetXaxis()->SetRangeUser(-1,1);
					//					eff2D[file][trig][pt][frame][0]->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
					eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);

					TH2F* eff_pass_2dhist;
					eff_pass_2dhist = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
					eff_pass_2dhist->SetName("eff_pass_2dhist_costheta");

					costheta_eff[file][trig][pt][frame][0] =(TH1F*)eff_pass_2dhist->ProjectionY();
					costheta_eff[file][trig][pt][frame][0]->SetName(Form("costheta_pass_%d_%d_%d_%d",file,trig,pt,frame));

					eff_pass_2dhist = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
					eff_pass_2dhist->SetName("eff_pass_2dhist_phi");

					phi_eff[file][trig][pt][frame][0] =(TH1F*)eff_pass_2dhist->ProjectionX();
					phi_eff[file][trig][pt][frame][0]->SetName(Form("phi_pass_%d_%d_%d_%d",file,trig,pt,frame));

					eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
					eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
					//					eff2D[file][trig][pt][frame][1]->GetXaxis()->SetRangeUser(-1,1);
					//					eff2D[file][trig][pt][frame][1]->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
					eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);
					eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);

					TH2F* eff_total_2dhist;
					eff_total_2dhist = (TH2F*)eff2D[file][trig][pt][frame][1]->Clone();
					eff_total_2dhist->SetName("eff_total_2dhist_costheta");

					costheta_eff[file][trig][pt][frame][1] =(TH1F*)eff_total_2dhist->ProjectionY();
					costheta_eff[file][trig][pt][frame][1]->SetName(Form("costheta_total_%d_%d_%d_%d",file,trig,pt,frame));

					eff_total_2dhist = (TH2F*)eff2D[file][trig][pt][frame][1]->Clone();
					eff_total_2dhist->SetName("eff_total_2dhist_phi");

					phi_eff[file][trig][pt][frame][1] =(TH1F*)eff_total_2dhist->ProjectionX();
					phi_eff[file][trig][pt][frame][1]->SetName(Form("phi_total_%d_%d_%d_%d",file,trig,pt,frame));

					canvas[trig][pt][frame]->cd(5);
					costheta[file][trig][pt][frame][0]->Draw();

					canvas[trig][pt][frame]->cd(9);
					phi[file][trig][pt][frame][0]->Draw();

					canvas[trig][pt][frame]->cd(6);
					efficiency1D[file][trig][pt][frame][0] = new TGraphAsymmErrors(costheta_eff[file][trig][pt][frame][0],costheta_eff[file][trig][pt][frame][1],"N"); 
					efficiency1D[file][trig][pt][frame][0]->SetName(Form("costheta_eff_%d_%d_%d_%d",file,trig,pt,frame));
					efficiency1D[file][trig][pt][frame][0]->Draw("ap");

					canvas[trig][pt][frame]->cd(10);
					efficiency1D[file][trig][pt][frame][1] = new TGraphAsymmErrors(phi_eff[file][trig][pt][frame][0],phi_eff[file][trig][pt][frame][1],"N");
					efficiency1D[file][trig][pt][frame][1]->SetName(Form("phi_eff_%d_%d_%d_%d",file,trig,pt,frame));
					efficiency1D[file][trig][pt][frame][1]->Draw("ap");

					eff1D[file][trig][pt][frame][0] = new TH1F(Form("eff_file%d_trig%d_pt%d_frame%d_phase%d",file,trig,pt,frame,0),Form("eff_file%d_trig%d_pt%d_frame%d_phase%d",file,trig,pt,frame,0),10,-1,1); 
					eff1D[file][trig][pt][frame][1] = new TH1F(Form("eff_file%d_trig%d_pt%d_frame%d_phase%d",file,trig,pt,frame,1),Form("eff_file%d_trig%d_pt%d_frame%d_phase%d",file,trig,pt,frame,1),10,-TMath::Pi(),TMath::Pi());
					for(int phase=0;phase<2;phase++){
						for(int npoint=0;npoint<10;npoint++){
							Double_t x,y;
							efficiency1D[file][trig][pt][frame][phase]->GetPoint(npoint,x,y);
							eff1D[file][trig][pt][frame][phase]->SetBinContent(npoint+1,y);
							eff1D[file][trig][pt][frame][phase]->SetBinError(npoint+1,efficiency1D[file][trig][pt][frame][phase]->GetErrorY(npoint));
						}
					}	

					costheta[file][trig][pt][frame][1] = (TH1F*)costheta[file][trig][pt][frame][0]->Clone(); 
					costheta[file][trig][pt][frame][1]->SetName(Form("corrected_costheta_file%d_trig%d_pt%d_frame%d",file,trig,pt,frame));
					costheta[file][trig][pt][frame][1]->Divide(eff1D[file][trig][pt][frame][0]);//correct raw data 
					phi[file][trig][pt][frame][1] = (TH1F*)phi[file][trig][pt][frame][0]->Clone(); 
					phi[file][trig][pt][frame][1]->SetName(Form("corrected_phi_file%d_trig%d_pt%d_frame%d",file,trig,pt,frame));
					phi[file][trig][pt][frame][1]->Divide(eff1D[file][trig][pt][frame][1]);//correct raw data 

					canvas[trig][pt][frame]->cd(7);
					costheta[file][trig][pt][frame][1]->Draw();
					costheta_eff[file][trig][pt][frame][1]->SetMinimum(costheta_eff[file][trig][pt][frame][0]->GetMinimum());
					//		costheta_eff[file][trig][pt][frame][1]->Draw();
					//		costheta_eff[file][trig][pt][frame][0]->Draw("same");
					canvas[trig][pt][frame]->cd(11);
					phi[file][trig][pt][frame][1]->Draw();
					phi_eff[file][trig][pt][frame][1]->SetMinimum(phi_eff[file][trig][pt][frame][0]->GetMinimum());
					//		phi_eff[file][trig][pt][frame][1]->Draw();
					//		phi_eff[file][trig][pt][frame][0]->Draw("same");

					if(pt==0) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf(",trigName[trig].Data(),frame));			
					else if(pt==NPT-1) canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf)",trigName[trig].Data(),frame));			
					else canvas[trig][pt][frame]->Print(Form("figures/%s_frame%d.pdf",trigName[trig].Data(),frame));			
				}
			}
		}
	}
}

void correctedpass(int drawpass){
	if(drawpass==1){
		TH2F* pass = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
		pass->SetName("pass");
		pass->Divide(eff2D[file][trig][pt][frame][2]);
		pass->Draw("colz");	
	}
}

void compare(){
	gStyle->SetPadRightMargin();

	TCanvas* comparison = new TCanvas("comparison","comparison",1200,800);
	TLegend* comparisonleg;

	comparison->Divide(3,2);
	TFile *outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");
	TGraphErrors* output[NTRIG][NPLOT];

	comparison->cd(1);
	output[0][0] = (TGraphErrors*)outputfile->Get("plot_HT0_theta_0_1eID");	
	output[2][0] = (TGraphErrors*)outputfile->Get("plot_HT2_theta_0_1eID");	
	output[0][0]->SetMarkerStyle(20);
	output[2][0]->SetMarkerStyle(22);
	output[0][0]->SetMarkerSize(1);
	output[2][0]->SetMarkerSize(1);
	output[0][0]->Draw("ap");
	output[2][0]->Draw("samep");
	lambda_parameters[2][0][0]->SetMarkerStyle(26);
	lambda_parameters[2][0][0]->SetMarkerColor(kBlue);
	lambda_parameters[2][0][0]->SetLineColor(kBlue);
	lambda_theta_hx->SetMarkerStyle(24);
	lambda_parameters[2][0][0]->Draw("samep");
	lambda_theta_hx->Draw("samep");
	comparisonleg = new TLegend(0.7,0.7,0.89,0.89);
	comparisonleg->AddEntry(output[0][0],"HT0","p");
	comparisonleg->AddEntry(output[2][0],"HT2","p");
	comparisonleg->AddEntry(lambda_theta_hx,"new HT0","p");
	comparisonleg->AddEntry(lambda_parameters[2][0][0],"new HT2","p");
	comparisonleg->Draw("same");

	comparison->cd(2);
	output[0][1] = (TGraphErrors*)outputfile->Get("plot_HT0_phi_0_1eID");
	output[2][1] = (TGraphErrors*)outputfile->Get("plot_HT2_phi_0_1eID");
	output[0][1]->SetMarkerStyle(20);
	output[2][1]->SetMarkerStyle(22);
	output[0][1]->SetMarkerSize(1);
	output[2][1]->SetMarkerSize(1);
	output[0][1]->Draw("ap");
	output[2][1]->Draw("samep");
	lambda_phi_hx->SetMarkerStyle(24);
	lambda_parameters[2][0][1]->SetMarkerStyle(26);
	lambda_parameters[2][0][1]->SetMarkerColor(kBlue);
	lambda_parameters[2][0][1]->SetLineColor(kBlue);
	lambda_phi_hx->Draw("samep");
	lambda_parameters[2][0][1]->Draw("samep");
	comparisonleg = new TLegend(0.7,0.7,0.89,0.89);
	comparisonleg->AddEntry(output[0][1],"HT0","p");
	comparisonleg->AddEntry(output[2][1],"HT2","p");
	comparisonleg->AddEntry(lambda_phi_hx,"new HT0","p");
	comparisonleg->AddEntry(lambda_parameters[2][0][1],"new HT2","p");
	comparisonleg->Draw("same");

	comparison->SaveAs("figures/comparison.pdf");
}

void compare_lambda(int trig,int frame, int variable){
	if(variable==0)	comparisonlambda[trig][variable+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
	if(variable==1) comparisonlambda[trig][variable+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));
	if(variable==2) comparisonlambda[trig][variable+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_inv_%d_1eID",trigName[trig].Data(),frame));

	if(trig==0){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kRed);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kRed);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(24);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==1){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kBlue);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kBlue);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(26);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==2){
		comparisonlambda[trig][variable+3*frame]->SetMarkerColor(kBlack);
		comparisonlambda[trig][variable+3*frame]->SetLineColor(kBlack);
		comparisonlambda[trig][variable+3*frame]->SetMarkerStyle(27);
		comparisonlambda[trig][variable+3*frame]->SetMarkerSize(1.2);
	}
	comparisonlambda[trig][variable+3*frame]->Draw("samep");
}

void crosscheck(){
	TFile *oldfile = new TFile("~/jpsi/test20160210_Barbara/efficiency1eID.root","read");
	TH1F *oldhist[NTRIG][NPT];
	TH1F *h1[6],*h2[NPT][3];

	TFile *defaultfile[3];// embedding ht0_1eID ht2_1eID
	TH1F *defaulthist[2],*defaultrawdata[3];

	TCanvas *canvas2 = new TCanvas("canvas2","canvas2",1200,800);
	canvas2->Divide(6,2);

	gStyle->SetOptStat(0);
	defaultfile[0] = new TFile("~/jpsi/test20160210_0/rootfile/OutFile_cent_0_9.root","read");
	defaultfile[1] = new TFile("~/jpsi/test20160210_0/rootfile/ht0_trg3_1TrkPid_0.ana.root","read");
	defaultfile[2] = new TFile("~/jpsi/test20160210_0/rootfile/ht2_trg5_1TrkPid_0.ana.root","read");

	int xmin,xmax,ymin,ymax,zmin,zmax;

	for(int trig=0;trig<1;trig++){
		for(int pt=0;pt<NPT;pt++){
			if((TH1F*)oldfile->Get(Form("hEff%sHXThetaPt%d",trigName[trig].Data(),pt))!=0x0) oldhist[trig][pt] = (TH1F*)oldfile->Get(Form("hEff%sHXThetaPt%d",trigName[trig].Data(),pt));	
			ymin = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(-TMath::Pi());
			ymax = eff3D[0][1][1][0][0]->GetYaxis()->FindBin(TMath::Pi());
			zmin = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt]);//marker
			zmax = eff3D[0][0][1][0][0]->GetZaxis()->FindBin(PtEdge[pt+1])-1;

			canvas2->cd(pt+1);
			defaulthist[0] = (TH1F*)((TH2F*)defaultfile[0]->Get(Form("hHt%dJpsiCosThetaPt1",trig)))->ProjectionX(Form("defaulthist_%d_%d",trig,pt),zmin,zmax);
			defaulthist[0]->SetTitle("comparison of pass between 3D projection and 2D projection");
			defaulthist[0]->GetXaxis()->SetTitle("cos#theta");

			defaulthist[1] = (TH1F*)((TH2F*)defaultfile[0]->Get("hMcJpsiCosThetaPt"))->ProjectionX(Form("Mcdefaulthist_%d_%d",trig,pt),zmin,zmax);
			defaulthist[1]->SetTitle("comparison of total between 3D projection and 2D projection");
			defaulthist[1]->GetXaxis()->SetTitle("cos#theta");
			h1[0] = (TH1F*)eff3D[0][0][1][0][0]->ProjectionX("h1",ymin,ymax,zmin,zmax);
			h1[0]->RebinX(4);
			h1[0]->SetMarkerStyle(20);

			h1[1] = (TH1F*)eff3D[0][1][1][0][1]->ProjectionX("h2",ymin,ymax,zmin,zmax);
			h1[1]->RebinX(4);
			h1[1]->SetMarkerStyle(20);

			//	h1[1]->Draw("p");
			//	defaulthist[0]->Draw("samep");

			TGraphAsymmErrors* test = new TGraphAsymmErrors(h1[0],h1[1],"N");
			//	test->Print("all");	
			//	h1[0]->Draw("p");
			oldhist[trig][pt]->SetTitle("efficiency comparison");
			oldhist[trig][pt]->GetXaxis()->SetTitle("cos#theta");
			oldhist[trig][pt]->Draw("p");
			//			oldhist[trig][pt]->Print("all");
			test->Draw("samep");
			//			test->Print("all");

			//	oldhist[0]->Print("all");
			//			canvas->cd(2);
			//			defaulthist[1]->Draw("p");
			//			h1[1]->Draw("samep");

			canvas2->cd(7+pt);
			defaultrawdata[0] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPt"))->ProjectionX("rawdata",zmin,zmax);
			defaultrawdata[0]->SetTitle("Unlike Sign Comparison");
			//	defaultrawdata[0]->Draw("p");	
			h2[pt][0] = (TH1F*)rawdata3D[0][0][2][0][0]->ProjectionX("h4",ymin,ymax,zmin,zmax);
			h2[pt][0]->RebinX(4);
			h2[pt][0]->SetMarkerStyle(24);
			//	h2[0]->Draw("same");

			//canvas->cd(5);
			defaultrawdata[1] = (TH1F*)((TH2F*)defaultfile[1]->Get("hJpsiCosThetaPtBG"))->ProjectionX("rawdatabkg",1,8);
			defaultrawdata[1]->SetTitle("Like Sign Comparison");
			//	defaultrawdata[1]->Draw("p");	
			h2[pt][1] = (TH1F*)rawdata3D[0][0][2][0][1]->ProjectionX("h5",1,41,1,8);
			h2[pt][1]->RebinX(4);
			h2[pt][1]->SetMarkerStyle(24);
			//	h2[1]->Draw("same");
			h2[pt][2] = (TH1F*)h2[pt][0]->Clone();
			h2[pt][2]->Add(h2[pt][1],-1);
			h2[pt][2]->SetTitle("raw data");
			h2[pt][2]->Draw();

			defaultrawdata[2] = (TH1F*)defaultrawdata[0]->Clone();
			defaultrawdata[2]->Add(defaultrawdata[1],-1);	
			defaultrawdata[2]->Draw("same");

		}
	}
	canvas2->SaveAs("figures/2D3Dcrosscheck.pdf");
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const Int_t nbins = 10;
	Int_t i;
	Double_t xx,yy;
	Double_t chisq=0,chisq1=0,chisq2=0;
	Double_t delta[2];
	Int_t nbin=0;
	for(i=0;i<nbins;i++){
		xx = x[i];
		yy = y[i];
		if(fabs(errorz1[i]>1e-6)) delta[0] = (z1[i]-func1(&xx,par))/errorz1[i];
		else delta[0]=0;
		if(fabs(errorz2[i]>1e-6)) delta[1] = (z2[i]-func2(&yy,par))/errorz2[i];
		else delta[1]=0;
		chisq += (delta[0]*delta[0])+(delta[1]*delta[1]);
		chisq1 += delta[0]*delta[0];
		chisq2 += delta[1]*delta[1];
	}
	f=chisq;
	chisquare1=chisq1;
	chisquare2=chisq2;
	chisquare=chisq;
}


Double_t func1(Double_t *x,Double_t *par)
{
	Double_t value1=par[0]*(1+par[1]*x[0]*x[0]);
	return value1;
}

Double_t func2(Double_t *y,Double_t *par)
{
	Double_t value2=par[3]*(1+2*par[2]*TMath::Cos(2*y[0])/(3+par[1]));
	return value2;
}

double getaverage(std::vector<double> inputvector){
	double sum = 0.;
	for(int i=0;i<inputvector.size();i++) sum += inputvector[i];
	sum = sum/inputvector.size();
	return sum;
}

double getmaximum(std::vector<double> inputvector){
	double max;
	max = *std::max_element(&inputvector[0],&inputvector[0]+(inputvector.size()));// marker
	return max;
}

void lambdaparameters(int sys3){// plot and write lambda parameters 

	lambdas = new TFile("lambdas.root","recreate");
	lambdas->cd();

	int selectedfile,file;
	if(sys3>=0 && sys3<=NFILE) file=sys3,selectedfile=sys3+1;
	else if(sys3==100) {
		file=0,selectedfile=NFILE;
	}
	else return;

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT],y_sys[NPT],theta[NPT],thetaerr[NPT],phi[NPT],phierr[NPT];

	for(;file<selectedfile;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){// marker add cs frame 
				for(int phase=0;phase<NPHASE+1;phase++){
					for(int pt=0;pt<NPT;pt++) {
						x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.;
						x_err[pt] = (PtEdge[pt+1]-PtEdge[pt])/2.;
						theta[pt] = fitparameters[file][trig][pt][frame][0];
						thetaerr[pt] = fitparameters[file][trig][pt][frame][1];
						phi[pt] = fitparameters[file][trig][pt][frame][2];
						phierr[pt] = fitparameters[file][trig][pt][frame][3];
						if(phase!=2){
							if(fitdata[trig][pt][frame][phase]->GetEntries()!=0){
								y[pt] = fitparameters[file][trig][pt][frame][phase*2];
								y_err[pt] = fitparameters[file][trig][pt][frame][phase*2+1];
								y_sys[pt] = systematic_error[trig][pt][frame][phase];
							}
							else{
								y[pt] = -100.;
								y_err[pt] = 0.;
								y_sys[pt] = 0.;
							}
						}
						else{//calculate lambda_invariant and its statistic error
							if(theta[pt]!=-100 && phi[pt]!=-100){		
								y[pt] = (theta[pt]+3*phi[pt])/(1-phi[pt]);
								y_err[pt] = TMath::Sqrt(theta[pt]*theta[pt]/((1-phi[pt])*(1-phi[pt]))+phi[pt]*phi[pt]*TMath::Power((3+theta[pt]/((1-phi[pt])*(1-phi[pt]))),2)+2*(3+theta[pt])/TMath::Power((1-phi[pt]),3)*covariant[trig][pt][frame]);
								y_sys[pt] = y[pt]-(fitparameters[0][trig][pt][frame][0]+3*fitparameters[0][trig][pt][frame][2])/(1-fitparameters[0][trig][pt][frame][2]);
							}
							else{
								y[pt] = -100.;
								y_err[pt] = 0.;
								y_sys[pt] = 0.;
							}
						}
					}
					lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
					lambda_parameters[trig][frame][phase]->SetName(Form("lambdas_%d_%d_%d",trig,frame,phase));
					lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
					lambda_parameters_sys[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_sys);
					lambda_parameters_sys[trig][frame][phase]->SetName(Form("lambdas_sys_%d_%d_%d",trig,frame,phase));
					lambda_parameters_sys[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
					lambda_parameters[trig][frame][phase]->Write();		
					lambda_parameters_sys[trig][frame][phase]->Write();		
				}
			}
		}
	}

	TCanvas* lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
	lambdacanvas->Divide(3,2);

	for(int frame=0;frame<NFRAME;frame++){// marker add cs frame 
		for(int phase=0;phase<NPHASE+1;phase++){
			for(int trig=0;trig<NTRIG;trig++){
				lambdacanvas->cd(frame*3+phase+1);
				legend = new TLegend(0.7,0.7,0.89,0.89);
				lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
				if(trig==0){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
				}
				if(trig==1){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
				}
				plotsetting();
				legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
				if(trig==0)	lambda_parameters[trig][frame][phase]->Draw("ap");
				else lambda_parameters[trig][frame][phase]->Draw("samep");
				//				if(phase!=2)
				compare_lambda(trig,frame,phase);
				//				lambda_parameters_sys[trig][frame][phase]->Draw("samep[]");
			}
		}
	}
	lambdacanvas->SaveAs("figures/lambdas.pdf");
}

void plotsetting(){
	lambda_parameters[0][0][0]->SetTitle("#lambda_{#theta} in HX frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][1][0]->SetTitle("#lambda_{#theta} in CS frame;J/#psi p_{T};#lambda_{#theta}");
	lambda_parameters[0][0][1]->SetTitle("#lambda_{#phi} in HX frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][1][1]->SetTitle("#lambda_{#phi} in CS frame;J/#psi p_{T};#lambda_{#phi}");
	lambda_parameters[0][0][2]->SetTitle("#lambda_{inv} in HX frame;J/#psi p_{T};#lambda_{inv}");
	lambda_parameters[0][1][2]->SetTitle("#lambda_{inv} in CS frame;J/#psi p_{T};#lambda_{inv}");
}


void printscheme(){
	cout<<" **************************************************************************************"<<"\n"
		<<" *  do2Dcorrection function                                                           *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  3D_unlike - 3D_like = 3D_raw_data projection to xy plane -> 2D_raw_data           *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  2D_raw_data project to costheta(phi) axis -> 1D_raw_costheta(phi)                 *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  3D_pass project to xy plane -> 2D_pass; 3D_total project to xy plane -> 2D_total  *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  2D_pass project to costheta and phi axis -> 1D_costheta_pass and 1D_phi_pass      *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  2D_total project to costheta and phi axis -> 1D_costheta_total and 1D_phi_total   *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  1D_costheta(phi)_pass / 1D_costheta(phi)_total = costheta(phi)_efficiency;        *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  1D_raw_costheta(phi) / costheta(phi)_efficiency = 1D_corrected_costheta(phi)      *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  simultaneously fit costheta and phi distribution with corresponding function      *"<<"\n"
		<<" *                                                                                    *"<<"\n"
		<<" *  extract lambda_theta, lambda_phi and lambda_invariant                             *"<<"\n"
		<<" **************************************************************************************"<<endl;
}

void simultaneousfit(int sys1 = 0){// marker add frame 
	static Double_t vstart[NTRIG][NPT][5];
	static Double_t step[NTRIG][NPT][5];

	int file,selectedfile;
	if(sys1>=0 && sys1<=NFILE) file=sys1,selectedfile=sys1+1;
	else if(sys1==100) file=0,selectedfile=NFILE;
	for(;file<selectedfile;file++){ //marker
		TCanvas *canvas1[NTRIG];
		for(int trig=0;trig<NTRIG;trig++){
			canvas1[trig]	= new TCanvas(Form("canvas1%d",trig),Form("canvas1%d",trig),1200,800);
			canvas1[trig]->Divide(NPT,4);
			for(int frame=0;frame<NFRAME;frame++){
				for(int pt=0;pt<NPT;pt++){
					//	int frame = 0;

					fitdata[trig][pt][frame][0] = (TH1F*)costheta[file][trig][pt][frame][1]->Clone();
					fitdata[trig][pt][frame][0]->SetName("costheta_corrected_data");
					fitdata[trig][pt][frame][0]->Scale(1./fitdata[trig][pt][frame][0]->Integral());

					fitdata[trig][pt][frame][1] = (TH1F*)phi[file][trig][pt][frame][1]->Clone();
					fitdata[trig][pt][frame][1]->SetName("phi_corrected_data");
					fitdata[trig][pt][frame][1]->Scale(1./fitdata[trig][pt][frame][1]->Integral());

					for(int nbin=1;nbin<40/NREBIN+1;nbin++){
						z1[nbin-1] = fitdata[trig][pt][frame][0]->GetBinContent(nbin);
						errorz1[nbin-1] = fitdata[trig][pt][frame][0]->GetBinError(nbin);
						x[nbin-1] = fitdata[trig][pt][frame][0]->GetBinCenter(nbin);
						//								if(fabs(z2[nbin-1])>1e-6)empty[trig][pt][frame][eID-1][1]--;
					}

					for(int nbin=1;nbin<40/NREBIN+1;nbin++){
						z2[nbin-1] = fitdata[trig][pt][frame][1]->GetBinContent(nbin);
						errorz2[nbin-1] = fitdata[trig][pt][frame][1]->GetBinError(nbin);
						y[nbin-1] = fitdata[trig][pt][frame][1]->GetBinCenter(nbin);
						//								if(fabs(z2[nbin-1])>1e-6)empty[trig][pt][frame][eID-1][1]--;
					}

					vstart[trig][pt][0]=0.1;
					vstart[trig][pt][1]=0.1;
					vstart[trig][pt][2]=0.1;
					vstart[trig][pt][3]=-0.1;
					vstart[trig][pt][4]=0.001;

					step[trig][pt][0]=0.001;
					step[trig][pt][1]=0.001;
					step[trig][pt][2]=0.001;
					step[trig][pt][3]=0.001;
					step[trig][pt][4]=0.001;

					arglist[0] = 1;
					TMinuit *gMinuit = new TMinuit(5);
					gMinuit->SetFCN(fcn);
					gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

					gMinuit->mnparm(0,"norm1",vstart[trig][pt][0],step[trig][pt][0],0,0,ierflg);
					gMinuit->mnparm(1,"lambda1",vstart[trig][pt][1],step[trig][pt][1],0,0,ierflg);
					fitflag[trig][pt][frame][1]=ierflg;
					gMinuit->mnparm(2,"lambda2",vstart[trig][pt][2],step[trig][pt][2],0,0,ierflg);
					fitflag[trig][pt][frame][2]=ierflg;
					gMinuit->mnparm(3,"norm2",vstart[trig][pt][3],step[trig][pt][3],0,0,ierflg);

					arglist[0] = 500.;
					arglist[1] = 1.;
					gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
					Double_t amin,edm,errdef;
					Int_t nvpar,nparx,icstat;
					gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

					for(int var=0;var<4;var++)gMinuit->GetParameter(var,fParamVal[file][trig][pt][frame][var],fParamErr[file][trig][pt][frame][var]);
					lambda_theta[file][trig][pt][0] = fParamVal[file][trig][pt][frame][1];
					lambda_theta_err[file][trig][pt][0] = fParamErr[file][trig][pt][frame][1];
					lambda_phi[file][trig][pt][0] = fParamVal[file][trig][pt][frame][2];
					lambda_phi_err[file][trig][pt][0] = fParamErr[file][trig][pt][frame][2];

					gMinuit->mnemat(&matrix[0][0],4);
					covariant[trig][pt][frame] = matrix[1][2];

					fitparameters[file][trig][pt][frame][0] = fParamVal[file][trig][pt][frame][1];
					fitparameters[file][trig][pt][frame][1] = fParamErr[file][trig][pt][frame][1];
					fitparameters[file][trig][pt][frame][2] = fParamVal[file][trig][pt][frame][2];
					fitparameters[file][trig][pt][frame][3] = fParamErr[file][trig][pt][frame][2];

					fit1[pt][frame][0] = new TF1(Form("fit1_%d_%d",pt,0),func1,-1,1,2);
					fit2[pt][frame][0] = new TF1(Form("fit2_%d_%d",pt,0),func2,-TMath::Pi(),TMath::Pi(),4);
					fit1[pt][frame][1] = new TF1(Form("fit1_%d_%d",pt,1),func1,-1,1,2);
					fit2[pt][frame][1] = new TF1(Form("fit2_%d_%d",pt,1),func2,-TMath::Pi(),TMath::Pi(),4);

					fit1[pt][frame][0]->SetParameters(fParamVal[0][trig][pt][frame][0],fParamVal[0][trig][pt][frame][1]);
					fit2[pt][frame][0]->SetParameters(fParamVal[0][trig][pt][frame][0],fParamVal[0][trig][pt][frame][1],fParamVal[0][trig][pt][frame][2],fParamVal[0][trig][pt][frame][3]);
					canvas1[trig]->cd(pt+1);
					if(fitdata[trig][pt][frame][0]->GetEntries()!=0){
						fitdata[trig][pt][frame][0]->Draw();
						fit1[pt][frame][0]->Draw("same");	
					}
					canvas1[trig]->cd(pt+NPT+1);
					if(fitdata[trig][pt][frame][1]->GetEntries()!=0){
						fitdata[trig][pt][frame][1]->Draw();
						fit2[pt][frame][0]->Draw("same");
					}
					canvas1[trig]->SaveAs(Form("figures/simultaneousfit%d.pdf",trig));
				}
			}
		}
	}
}

void systematic(int sys = 0){

	correcteddata(0);
	simultaneousfit(0);

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT];

	systematic[0] = new TCanvas("systematic","systematic",1200,800);
	systematic[0]->Divide(4,2);
	for(int i=0;i<NPT;i++){
		x[i] = (PtEdge[i]+PtEdge[i+1])/2.;
		x_err[i] = (PtEdge[i+1]-PtEdge[i])/2.;
		y[i] = fit2dtheta[0][0][i][0];
		y_err[i] = fit2dtheta_err[0][0][i][0];
		cout<<"lambda_theta======"<<y[i]<<endl;	
	}

	lambda_theta_hx =  new TGraphErrors(NPT,x,y,x_err,y_err);
	systematic[0]->cd(1);
	lambda_theta_hx->GetYaxis()->SetRangeUser(-2,2);
	lambda_theta_hx->Draw("ape");
	//	systematic[0]->SaveAs("systematic1.pdf");	

	for(int i=0;i<NPT;i++){
		x[i] = (PtEdge[i]+PtEdge[i+1])/2.;
		x_err[i] = (PtEdge[i+1]-PtEdge[i])/2.;
		y[i] = fit2dphi[0][0][i][0];
		y_err[i] = fit2dphi_err[0][0][i][0];

		//		y[i] = lambda_phi[0][0][i][0];
		//		y_err[i] = lambda_phi_err[0][0][i][0];	
		cout<<"lambda_phi======"<<y[i]<<endl;	
	}

	lambda_phi_hx =  new TGraphErrors(NPT,x,y,x_err,y_err);
	systematic[0]->cd(2);
	lambda_phi_hx->GetYaxis()->SetRangeUser(-2,2);
	lambda_phi_hx->Draw("ape");
	systematic[0]->SaveAs("figures/systematic1.pdf");	

	for(int file=0;file<NFILE;file++){
		for(int trig=0;trig<NTRIG;trig++){
			for(int frame=0;frame<NFRAME;frame++){
				for(int phase=0;phase<NPHASE;phase++){
					for(int pt=0;pt<NPT;pt++) {
						y[pt] = fitparameters[file][trig][pt][frame][phase*2];
						y_err[pt] = fitparameters[file][trig][pt][frame][phase*2+1];
					}
					lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);
				}
			}
		}
	}

	if(sys>=1 && sys<=NFILE){
		correcteddata(sys);
		simultaneousfit(sys);
		int trig=0;
		for(int pt=0;pt<5;pt++) cout<<"lambda_theta in "<<sys<<"   "<<trig<<"  pt "<<pt<<"    "<<lambda_theta[sys][trig][pt][0]<<"&&&&&&&&&"<<lambda_theta[0][trig][pt][0]<<endl;
		//		cout<<"systematic uncertainty from "<<sys<<endl;//marker
	}
	else if(sys==100){
		//		for(sys=1;sys<NFILE;sys++){
		correcteddata(sys);
		simultaneousfit(sys);
		//		}
		std::vector<double> sysvector_theta_hx,sysvector_theta_cs,sysvector_phi_hx,sysvector_phi_cs;

		int file,trig,pt,frame;

		for(trig=0;trig<NTRIG;trig++){//marker add the definition of systematic 
			for(pt=0;pt<NPT;pt++){
				for(file=1;file<12;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);// revise the frame 
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);// revise the frame 
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				ptsmearing[trig][pt][0][0] = getaverage(sysvector_theta_hx);// marker getmaximum()
				ptsmearing[trig][pt][1][0] = getaverage(sysvector_theta_cs);// marker getmaximum()
				ptsmearing[trig][pt][0][1] = getaverage(sysvector_phi_hx);// marker getmaximum()
				ptsmearing[trig][pt][1][1] = getaverage(sysvector_phi_cs);// marker getmaximum()
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"ptsmearing = "<<ptsmearing[trig][pt][0][0]<<endl;

				for(file=12;file<20;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				weight[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				weight[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				weight[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				weight[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"weight = "<<weight[trig][pt][0][0]<<endl;

				for(file=20;file<24;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				nhitsfit[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
				nhitsfit[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
				nhitsfit[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
				nhitsfit[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"nhitsfit = "<<nhitsfit[trig][pt][0][0]<<endl;

				for(file=24;file<29;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				nsigma[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				nsigma[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				nsigma[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				nsigma[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"nsigma = "<<nsigma[trig][pt][0][0]<<endl;

				for(file=29;file<31;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				dsmadc[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				dsmadc[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				dsmadc[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				dsmadc[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"dsmadc = "<<dsmadc[trig][pt][0][0]<<endl;

				for(file=31;file<35;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				beta[trig][pt][0][0] = getmaximum(sysvector_theta_hx);
				beta[trig][pt][1][0] = getmaximum(sysvector_theta_cs);
				beta[trig][pt][0][1] = getmaximum(sysvector_phi_hx);
				beta[trig][pt][1][1] = getmaximum(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"beta = "<<beta[trig][pt][0][0]<<endl;

				for(file=35;file<36;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				poe[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				poe[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				poe[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				poe[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"poe = "<<poe[trig][pt][0][0]<<endl;

				for(file=36;file<40;file++){
					cout<<"lamba1 "<<lambda_theta[file][trig][pt][0]<<" - lambda2 "<<lambda_theta[0][trig][pt][0]<<"====="<<lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]<<endl;
					sysvector_theta_hx.push_back(lambda_theta[file][trig][pt][0]-lambda_theta[0][trig][pt][0]);
					sysvector_theta_cs.push_back(lambda_theta[file][trig][pt][1]-lambda_theta[0][trig][pt][1]);
					sysvector_phi_hx.push_back(lambda_phi[file][trig][pt][0]-lambda_phi[0][trig][pt][0]);
					sysvector_phi_cs.push_back(lambda_phi[file][trig][pt][1]-lambda_phi[0][trig][pt][1]);
				}
				pol[trig][pt][0][0] = getaverage(sysvector_theta_hx);
				pol[trig][pt][1][0] = getaverage(sysvector_theta_cs);
				pol[trig][pt][0][1] = getaverage(sysvector_phi_hx);
				pol[trig][pt][1][1] = getaverage(sysvector_phi_cs);
				sysvector_theta_hx.clear();
				sysvector_theta_cs.clear();
				sysvector_phi_hx.clear();
				sysvector_phi_cs.clear();
				cout<<"pol = "<<pol[trig][pt][0][0]<<endl;

				cout<<"combining systematic uncertainty"<<endl;//marker

				for(int frame=0;frame<NFRAME;frame++){		
					for(int phase=0;phase<NPHASE;phase++){
						systematic_error[trig][pt][frame][phase] = TMath::Sqrt(fabs(ptsmearing[trig][pt][frame][phase])*fabs(ptsmearing[trig][pt][frame][phase])+fabs(weight[trig][pt][frame][phase])*fabs(weight[trig][pt][frame][phase])+fabs(nhitsfit[trig][pt][frame][phase])*fabs(nhitsfit[trig][pt][frame][phase])+fabs(nsigma[trig][pt][frame][phase])*fabs(nsigma[trig][pt][frame][phase])+fabs(dsmadc[trig][pt][frame][phase])*fabs(dsmadc[trig][pt][frame][phase])+fabs(beta[trig][pt][frame][phase])*fabs(beta[trig][pt][frame][phase])+fabs(poe[trig][pt][frame][phase])*fabs(poe[trig][pt][frame][phase])+fabs(pol[trig][pt][frame][phase])*fabs(pol[trig][pt][frame][phase]));
					}
				}
				cout<<"combining systematic uncertainty is done"<<endl;
			}
		}
	}
}

