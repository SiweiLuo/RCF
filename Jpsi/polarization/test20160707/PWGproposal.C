
TString trigSet[3] = {"MB","HT0","HT2"};
TString histType[3] = {"unlikeSign","likeSign","eff"};

void PWGproposal(){

	TFile* infile[3];
/*   	infile[0] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_0_0_0_0_20171025.root","read");
	infile[1] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_1_1_0_0_20171025.root","read");
	infile[2] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_3_4_0_0_20171025.root","read");
*/

	infile[0] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht-1_sys0.ana.root","read");
	infile[1] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht0_sys0.ana.root","read");
	infile[2] = new TFile("/star/u/siwei/polarization/invmass_2eid_inv30/rootfile20170922/ht2_sys0.ana.root","read");

	TH2F *unlike[3],*like[3],*efficiency[3];
	//unlike[0] = infile[0]->Get("tothist");
	//unlike[0]->SetName("MB_unlike_sign");

	for(int i=0;i<3;i++){
		unlike[i] = (TH2F*)((TH3F*)infile[i]->Get("hJpsiCosThetaPhiPt"))->Project3D("xy");
		unlike[i]->SetName(Form("%s_unlike_sign",trigSet[i].Data()));
		like[i] = (TH2F*)((TH3F*)infile[i]->Get("hJpsiCosThetaPhiPtBG"))->Project3D("xy");
		like[i]->SetName(Form("%s_like_sign",trigSet[i].Data()));
//		efficiency[i] = (TH2F*)infile[i]->Get("effhist");
//		efficiency[i]->SetName(Form("%s_eff",trigSet[i].Data()));
	}

	TCanvas *dataeff[3];
	for(int i=0;i<3;i++){
		dataeff[i] = new TCanvas(Form("dataeff_%d",i),Form("dataeff_%d",i),1200,400);
		dataeff[i]->Divide(3,1);
		for(int j=0;j<3;j++){
			dataeff[i]->cd(j+1);	
			if(i==0) unlike[j]->Draw("colz");
			if(i==1) like[j]->Draw("colz");
//			if(i==2) efficiency[j]->Draw("colz");
		}
		dataeff[i]->SaveAs(Form("~/WWW/Proposal/figure/%s.pdf",histType[i].Data()));
	}	

	TFile *infilept[3];

	infilept[0] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_0_0_0_0_20171025.root","read");
	infilept[1] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_1_1_0_0_20171025.root","read");
	infilept[2] = new TFile("~/polresults/20160707/functional/splot_3D_20171025_2eid_inv30/splot_3D_3_4_0_0_20171025.root","read");




}
