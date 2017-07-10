#include "TCanvas.h"

TCanvas* canvas;

TGraphErrors* staterr[2][2][3];
//TGraphErrors* staterr[2][3][2];
TGraphErrors* systerr[2][2][3];
//TGraphErrors* ht2systerr[2][3][2];

void finalize(){

	TFile* statistic = new TFile("~/polresults/20160707/functional/splot_3D_20170510/lambdasvspt.root","read");
	TFile* systematic = new TFile("outtest.root","read");

	canvas = new TCanvas("canvas","canvas",2400,1600);
	canvas->Divide(3,2,0.01,0.01);

	//	canvas->cd(1);
	staterr[0][0][0] = (TGraphErrors*)statistic->Get("lambdas_0_0_0");
	staterr[0][0][0]->SetName("staterr_0_0_0");
	//	staterr[0][0][0]->GetYaxis()->SetRangeUser(-2,2);
	//	staterr[0][0][0]->Draw("ap");

	systerr[0][0][0] = (TGraphErrors*)systematic->Get("lambda_0_0_0");
	systerr[0][0][0]->SetName("systerr_0_0_0");
	systerr[0][0][0]->SetMarkerColor(kRed);
	systerr[0][0][0]->SetLineColor(kRed);
	//	systerr[0][0][0]->Draw("psame");

	staterr[1][0][0] = (TGraphErrors*)statistic->Get("lambdas_2_0_0");
	staterr[1][0][0]->SetName("staterr_2_0_0");
	staterr[1][0][0]->SetMarkerColor(kRed);
	staterr[1][0][0]->SetLineColor(kRed);
	//	staterr[1][0][0]->Draw("psame");

	systerr[1][0][0] = (TGraphErrors*)systematic->Get("lambda_2_0_0");
	systerr[1][0][0]->SetName("systerr_2_0_0");
	systerr[1][0][0]->SetMarkerColor(kRed);
	systerr[1][0][0]->SetLineColor(kRed);
	//	systerr[1][0][0]->SetMarkerSize(1.2);
	//	systerr[1][0][0]->Draw("psame");

	staterr[0][0][1] = (TGraphErrors*)statistic->Get("lambdas_0_1_0");
	staterr[0][0][1]->SetName("staterr_0_0_1");
	staterr[0][0][1]->SetMarkerColor(kRed);
	staterr[0][0][1]->SetLineColor(kRed);

	staterr[1][0][1] = (TGraphErrors*)statistic->Get("lambdas_2_1_0");
	staterr[1][0][1]->SetName("staterr_2_0_1");
	staterr[1][0][1]->SetMarkerColor(kRed);
	staterr[1][0][1]->SetLineColor(kRed);

	systerr[0][0][1] = (TGraphErrors*)systematic->Get("lambda_0_0_1");
	systerr[0][0][1]->SetName("systerr_0_0_1");
	systerr[0][0][1]->SetMarkerColor(kRed);
	systerr[0][0][1]->SetLineColor(kRed);

	systerr[1][0][1] = (TGraphErrors*)systematic->Get("lambda_2_0_1");
	systerr[1][0][1]->SetName("systerr_2_0_1");
	systerr[1][0][1]->SetMarkerColor(kRed);
	systerr[1][0][1]->SetLineColor(kRed);


	staterr[0][0][2] = (TGraphErrors*)statistic->Get("lambdas_0_2_0");
	staterr[0][0][2]->SetName("staterr_0_2_0");
	staterr[0][0][2]->SetMarkerColor(kRed);
	staterr[0][0][2]->SetLineColor(kRed);

	staterr[1][0][2] = (TGraphErrors*)statistic->Get("lambdas_2_2_0");
	staterr[1][0][2]->SetName("staterr_2_2_0");
	staterr[1][0][2]->SetMarkerColor(kRed);
	staterr[1][0][2]->SetLineColor(kRed);


	systerr[0][0][2] = (TGraphErrors*)systematic->Get("lambda_0_0_2");
	systerr[0][0][2]->SetName("systerr_0_0_2");
	systerr[0][0][2]->SetMarkerColor(kRed);
	systerr[0][0][2]->SetLineColor(kRed);


	systerr[1][0][2] = (TGraphErrors*)systematic->Get("lambda_2_0_2");
	systerr[1][0][2]->SetName("systerr_2_0_2");
	systerr[1][0][2]->SetMarkerColor(kRed);
	systerr[1][0][2]->SetLineColor(kRed);

	staterr[0][1][0] = (TGraphErrors*)statistic->Get("lambdas_0_0_1");
	staterr[0][1][0]->SetName("systerr_0_0_1");
	staterr[0][1][0]->SetMarkerColor(kRed);
	staterr[0][1][0]->SetLineColor(kRed);

	staterr[1][1][0] = (TGraphErrors*)statistic->Get("lambdas_2_0_1");
	staterr[1][1][0]->SetName("systerr_2_0_1");
	staterr[1][1][0]->SetMarkerColor(kRed);
	staterr[1][1][0]->SetLineColor(kRed);

	systerr[0][1][0] = (TGraphErrors*)systematic->Get("lambda_0_1_0");
	systerr[0][1][0]->SetName("systerr_0_1_0");
	systerr[0][1][0]->SetMarkerColor(kRed);
	systerr[0][1][0]->SetLineColor(kRed);

	systerr[1][1][0] = (TGraphErrors*)systematic->Get("lambda_2_1_0");
	systerr[1][1][0]->SetName("systerr_2_1_0");
	systerr[1][1][0]->SetMarkerColor(kRed);
	systerr[1][1][0]->SetLineColor(kRed);

	staterr[0][1][1] = (TGraphErrors*)statistic->Get("lambdas_0_1_1");
	staterr[0][1][1]->SetName("staterr_0_1_1");
	staterr[0][1][1]->SetMarkerColor(kRed);
	staterr[0][1][1]->SetLineColor(kRed);

	staterr[1][1][1] = (TGraphErrors*)statistic->Get("lambdas_2_1_1");
	staterr[1][1][1]->SetName("staterr_2_1_1");
	staterr[1][1][1]->SetMarkerColor(kRed);
	staterr[1][1][1]->SetLineColor(kRed);

	systerr[0][1][1] = (TGraphErrors*)systematic->Get("lambda_0_1_1");
	systerr[0][1][1]->SetName("systerr_0_1_1");
	systerr[0][1][1]->SetMarkerColor(kRed);
	systerr[0][1][1]->SetLineColor(kRed);

	systerr[1][1][1] = (TGraphErrors*)systematic->Get("lambda_2_1_1");
	systerr[1][1][1]->SetName("systerr_2_1_1");
	systerr[1][1][1]->SetMarkerColor(kRed);
	systerr[1][1][1]->SetLineColor(kRed);

	staterr[0][1][2] = (TGraphErrors*)statistic->Get("lambdas_0_2_1");
	staterr[0][1][2]->SetName("staterr_0_2_1");
	staterr[0][1][2]->SetMarkerColor(kRed);
	staterr[0][1][2]->SetLineColor(kRed);

	staterr[1][1][2] = (TGraphErrors*)statistic->Get("lambdas_2_2_1");
	staterr[1][1][2]->SetName("staterr_2_2_1");
	staterr[1][1][2]->SetMarkerColor(kRed);
	staterr[1][1][2]->SetLineColor(kRed);

	systerr[0][1][2] = (TGraphErrors*)systematic->Get("lambda_0_1_2");
	systerr[0][1][2]->SetName("systerr_0_1_2");
	systerr[0][1][2]->SetMarkerColor(kRed);
	systerr[0][1][2]->SetLineColor(kRed);

	systerr[1][1][2] = (TGraphErrors*)systematic->Get("lambda_2_1_2");
	systerr[1][1][2]->SetName("systerr_2_1_2");
	systerr[1][1][2]->SetMarkerColor(kRed);
	systerr[1][1][2]->SetLineColor(kRed);

	title();
	removepoints();
	draw();

	canvas->SaveAs("final_20170709.pdf");
}


void title(){
	systerr[0][0][0]->SetTitle(";p_{T} GeV/c;#lambda_{#theta} in HX frame");
	systerr[0][0][1]->SetTitle(";p_{T} GeV/c;#lambda_{#phi} in HX frame");
	systerr[0][0][2]->SetTitle(";p_{T} GeV/c;#lambda_{#theta#phi} in HX frame");
	systerr[0][1][0]->SetTitle(";p_{T} GeV/c;#lambda_{#theta} in CS frame");
	systerr[0][1][1]->SetTitle(";p_{T} GeV/c;#lambda_{#phi} in CS frame");
	systerr[0][1][2]->SetTitle(";p_{T} GeV/c;#lambda_{#theta#phi} in CS frame");
}

void removepoints(){

	staterr[0][0][0]->RemovePoint(5);
	staterr[0][0][0]->RemovePoint(4);
	staterr[0][0][0]->RemovePoint(3);
	staterr[0][0][0]->RemovePoint(0);

	systerr[0][0][0]->RemovePoint(5);
	systerr[0][0][0]->RemovePoint(4);
	systerr[0][0][0]->RemovePoint(3);
	systerr[0][0][0]->RemovePoint(0);

	staterr[1][0][0]->RemovePoint(5);
	staterr[1][0][0]->RemovePoint(2);
	staterr[1][0][0]->RemovePoint(1);
	staterr[1][0][0]->RemovePoint(0);

	systerr[1][0][0]->RemovePoint(5);
	systerr[1][0][0]->RemovePoint(2);
	systerr[1][0][0]->RemovePoint(1);
	systerr[1][0][0]->RemovePoint(0);


	staterr[0][0][1]->RemovePoint(5);
	staterr[0][0][1]->RemovePoint(4);
	staterr[0][0][1]->RemovePoint(3);
	staterr[0][0][1]->RemovePoint(0);

	systerr[0][0][1]->RemovePoint(5);
	systerr[0][0][1]->RemovePoint(4);
	systerr[0][0][1]->RemovePoint(3);
	systerr[0][0][1]->RemovePoint(0);

	staterr[1][0][1]->RemovePoint(5);
	staterr[1][0][1]->RemovePoint(2);
	staterr[1][0][1]->RemovePoint(1);
	staterr[1][0][1]->RemovePoint(0);

	systerr[1][0][1]->RemovePoint(5);
	systerr[1][0][1]->RemovePoint(2);
	systerr[1][0][1]->RemovePoint(1);
	systerr[1][0][1]->RemovePoint(0);

	staterr[0][0][2]->RemovePoint(5);
	staterr[0][0][2]->RemovePoint(4);
	staterr[0][0][2]->RemovePoint(3);
	staterr[0][0][2]->RemovePoint(0);

	systerr[0][0][2]->RemovePoint(5);
	systerr[0][0][2]->RemovePoint(4);
	systerr[0][0][2]->RemovePoint(3);
	systerr[0][0][2]->RemovePoint(0);

	staterr[1][0][2]->RemovePoint(5);
	staterr[1][0][2]->RemovePoint(2);
	staterr[1][0][2]->RemovePoint(1);
	staterr[1][0][2]->RemovePoint(0);

	systerr[1][0][2]->RemovePoint(5);
	systerr[1][0][2]->RemovePoint(2);
	systerr[1][0][2]->RemovePoint(1);
	systerr[1][0][2]->RemovePoint(0);

	staterr[0][1][0]->RemovePoint(5);
	staterr[0][1][0]->RemovePoint(4);
	staterr[0][1][0]->RemovePoint(3);
	staterr[0][1][0]->RemovePoint(0);

	systerr[0][1][0]->RemovePoint(5);
	systerr[0][1][0]->RemovePoint(4);
	systerr[0][1][0]->RemovePoint(3);
	systerr[0][1][0]->RemovePoint(0);

	staterr[1][1][0]->RemovePoint(5);
	staterr[1][1][0]->RemovePoint(2);
	staterr[1][1][0]->RemovePoint(1);
	staterr[1][1][0]->RemovePoint(0);

	systerr[1][1][0]->RemovePoint(5);
	systerr[1][1][0]->RemovePoint(2);
	systerr[1][1][0]->RemovePoint(1);
	systerr[1][1][0]->RemovePoint(0);


	staterr[0][1][1]->RemovePoint(5);
	staterr[0][1][1]->RemovePoint(4);
	staterr[0][1][1]->RemovePoint(3);
	staterr[0][1][1]->RemovePoint(0);

	systerr[0][1][1]->RemovePoint(5);
	systerr[0][1][1]->RemovePoint(4);
	systerr[0][1][1]->RemovePoint(3);
	systerr[0][1][1]->RemovePoint(0);

	staterr[1][1][1]->RemovePoint(5);
	staterr[1][1][1]->RemovePoint(2);
	staterr[1][1][1]->RemovePoint(1);
	staterr[1][1][1]->RemovePoint(0);

	systerr[1][1][1]->RemovePoint(5);
	systerr[1][1][1]->RemovePoint(2);
	systerr[1][1][1]->RemovePoint(1);
	systerr[1][1][1]->RemovePoint(0);


	staterr[0][1][2]->RemovePoint(5);
	staterr[0][1][2]->RemovePoint(4);
	staterr[0][1][2]->RemovePoint(3);
	staterr[0][1][2]->RemovePoint(0);

	systerr[0][1][2]->RemovePoint(5);
	systerr[0][1][2]->RemovePoint(4);
	systerr[0][1][2]->RemovePoint(3);
	systerr[0][1][2]->RemovePoint(0);

	staterr[1][1][2]->RemovePoint(5);
	staterr[1][1][2]->RemovePoint(2);
	staterr[1][1][2]->RemovePoint(1);
	staterr[1][1][2]->RemovePoint(0);

	systerr[1][1][2]->RemovePoint(5);
	systerr[1][1][2]->RemovePoint(2);
	systerr[1][1][2]->RemovePoint(1);
	systerr[1][1][2]->RemovePoint(0);



}

void draw(){
	canvas->cd(1);
	systerr[0][0][0]->SetMarkerStyle(29);
	systerr[1][0][0]->SetMarkerStyle(29);
	staterr[0][0][0]->SetMarkerStyle(29);
	staterr[1][0][0]->SetMarkerStyle(29);

	systerr[0][0][0]->SetMarkerSize(1);
	systerr[1][0][0]->SetMarkerSize(1);
	staterr[0][0][0]->SetMarkerSize(1);
	staterr[1][0][0]->SetMarkerSize(1);
	/*
	   systerr[0][0][0]->SetMarkerStyle(29);
	   systerr[1][0][0]->SetMarkerStyle(29);
	   systerr[0][0][0]->SetMarkerSize(1);
	   systerr[1][0][0]->SetMarkerSize(1);
	   */
	systerr[0][0][0]->GetYaxis()->SetRangeUser(-1.,1.);
	systerr[0][0][0]->GetXaxis()->SetLimits(0,9);

	systerr[0][0][0]->Draw("ap[]");
	systerr[1][0][0]->Draw("[]psame");
	staterr[0][0][0]->Draw("psame");
	staterr[1][0][0]->Draw("psame");

	canvas->cd(2);
	systerr[0][0][1]->GetYaxis()->SetRangeUser(-1.,1.);
	systerr[0][0][1]->GetXaxis()->SetLimits(0,9);
	systerr[0][0][1]->SetMarkerStyle(29);
	systerr[1][0][1]->SetMarkerStyle(29);
	staterr[0][0][1]->SetMarkerStyle(29);
	staterr[1][0][1]->SetMarkerStyle(29);

	systerr[0][0][1]->SetMarkerSize(1);
	systerr[1][0][1]->SetMarkerSize(1);
	staterr[0][0][1]->SetMarkerSize(1);
	staterr[1][0][1]->SetMarkerSize(1);

	systerr[0][0][1]->Draw("ap[]");
	systerr[1][0][1]->Draw("psame[]");
	staterr[0][0][1]->Draw("psame");
	staterr[1][0][1]->Draw("psame");

	canvas->cd(3);
	systerr[0][0][2]->GetYaxis()->SetRangeUser(-1.,1.);
	systerr[0][0][2]->GetXaxis()->SetLimits(0,9);
	systerr[0][0][2]->SetMarkerStyle(29);
	systerr[1][0][2]->SetMarkerStyle(29);
	staterr[0][0][2]->SetMarkerStyle(29);
	staterr[1][0][2]->SetMarkerStyle(29);

	systerr[0][0][2]->SetMarkerSize(1);
	systerr[1][0][2]->SetMarkerSize(1);
	staterr[0][0][2]->SetMarkerSize(1);
	staterr[1][0][2]->SetMarkerSize(1);

	systerr[0][0][2]->Draw("ap[]");
	systerr[1][0][2]->Draw("psame[]");
	staterr[0][0][2]->Draw("psame");
	staterr[1][0][2]->Draw("psame");

	systerr[0][0][2]->Print();
	systerr[1][0][2]->Print();

	canvas->cd(4);
	systerr[0][1][0]->GetYaxis()->SetRangeUser(-1.5,1.5);
	systerr[0][1][0]->GetXaxis()->SetLimits(0,9);
	systerr[0][1][0]->SetMarkerStyle(29);
	systerr[1][1][0]->SetMarkerStyle(29);

	staterr[0][1][0]->SetMarkerStyle(29);
	staterr[1][1][0]->SetMarkerStyle(29);


	systerr[0][1][0]->SetMarkerSize(1);
	systerr[1][1][0]->SetMarkerSize(1);
	staterr[0][1][0]->SetMarkerSize(1);
	staterr[1][1][0]->SetMarkerSize(1);

	systerr[0][1][0]->Draw("ap[]");
	systerr[1][1][0]->Draw("psame[]");
	staterr[0][1][0]->Draw("psame");
	staterr[1][1][0]->Draw("psame");

	systerr[0][1][0]->Print();
	systerr[1][1][0]->Print();

	canvas->cd(5);
	systerr[0][1][1]->GetYaxis()->SetRangeUser(-1.,1.);
	systerr[0][1][1]->GetXaxis()->SetLimits(0,9);
	systerr[0][1][1]->SetMarkerStyle(29);
	systerr[1][1][1]->SetMarkerStyle(29);
	staterr[0][1][1]->SetMarkerStyle(29);
	staterr[1][1][1]->SetMarkerStyle(29);


	systerr[0][1][1]->SetMarkerSize(1);
	systerr[1][1][1]->SetMarkerSize(1);
	staterr[0][1][1]->SetMarkerSize(1);
	staterr[1][1][1]->SetMarkerSize(1);

	systerr[0][1][1]->Draw("ap[]");
	systerr[1][1][1]->Draw("psame[]");
	staterr[0][1][1]->Draw("psame");
	staterr[1][1][1]->Draw("psame");

	canvas->cd(6);
	systerr[0][1][2]->GetYaxis()->SetRangeUser(-1.,1.);
	systerr[0][1][2]->GetXaxis()->SetLimits(0,9);
	systerr[0][1][2]->SetMarkerStyle(29);
	systerr[1][1][2]->SetMarkerStyle(29);
	staterr[0][1][2]->SetMarkerStyle(29);
	staterr[1][1][2]->SetMarkerStyle(29);


	systerr[0][1][2]->SetMarkerSize(1);
	systerr[1][1][2]->SetMarkerSize(1);
	staterr[0][1][2]->SetMarkerSize(1);
	staterr[1][1][2]->SetMarkerSize(1);

	systerr[0][1][2]->Draw("ap[]");
	systerr[1][1][2]->Draw("psame[]");
	staterr[0][1][2]->Draw("psame");
	staterr[1][1][2]->Draw("psame");
}

/*systematic(){
}*/


