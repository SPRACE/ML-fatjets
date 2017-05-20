#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

TGraph* makeROCCurvesTau21(std::string fsignalName, std::string fbackgroundName){

	TCanvas* cv = new TCanvas("cv","cv",600,600);
	TFile* fsignal = TFile::Open(fsignalName.c_str());
	TFile* fbackground = TFile::Open(fbackgroundName.c_str());
	TH1* hsig = (TH1*)fsignal->Get("tau21");
	TH1* hbkg = (TH1*)fbackground->Get("tau21");
	hsig->Smooth();
	hbkg->Smooth();
	double totalSig = hsig->Integral();
	double totalBkg = hbkg->Integral();

	int nBinsX = hsig->GetNbinsX();
	TGraph* ROCGraph = new TGraph(nBinsX); 

	for(int i=1; i<=nBinsX; ++i) {
		double effSig = hsig->Integral(1,i)/totalSig;
		double effBkg = hbkg->Integral(1,i)/totalBkg;
		ROCGraph->SetPoint(i-1,effBkg,effSig);
	}

	ROCGraph->SetLineColor(kRed);
	ROCGraph->SetMarkerColor(kRed);
	ROCGraph->SetMarkerStyle(kFullSquare);
	ROCGraph->SetMarkerSize(0.5);
	ROCGraph->Draw("ALP");
	std::cout << "Area under the curve = " << ROCGraph->Integral(0,1) << std::endl;
	ROCGraph->GetYaxis()->SetTitle("Signal efficiency");
	ROCGraph->GetXaxis()->SetTitle("Background efficiency");
	TLine* line = new TLine(0,0,1,1);
	line->SetLineWidth(2);
	line->SetLineStyle(kDashed);
	line->Draw();
	double integral = 0.5+ROCGraph->Integral();
	TText *th1 = new TText(0.05,1,TString::Format("Area = %g",integral).Data());
    th1->SetTextAlign(11); th1->SetTextSize(0.04);
   	th1->Draw();
   	cv->SetGridx();
   	cv->SetGridy();
	cv->SaveAs("canvas.pdf");
	return ROCGraph;
}