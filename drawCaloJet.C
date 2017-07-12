#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

void drawCaloJet(int jet, int nEtaBins, int nPhiBins, double jetRadius, const char* fileName,
                const char* imageName = "jet.png") {

    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                nEtaBins,-jetRadius,jetRadius,
                                nPhiBins,-jetRadius,jetRadius);

    std::ifstream inFile;
    std::string line;
    inFile.open(fileName);
    for(int i=0; i!=jet; ++i) {
        if(not std::getline(inFile, line)) {return;}
    }

    double content;
    std::istringstream iss(line);
    std::cout << line <<  std::endl;
    for(int iEta = 1; iEta <= nEtaBins; ++iEta) {
        for(int iPhi = 1; iPhi <= nPhiBins; ++iPhi) {
            iss >> content;
            calorimeter->SetBinContent(iEta,iPhi,content);
        }
    }
    TCanvas* cv = new TCanvas("cv","cv",600,600);
    TStyle st;
    st.SetPalette(kInvertedDarkBodyRadiator);
    st.SetOptStat(0);
    cv->SetRightMargin(0.15);
    st.cd();
    calorimeter->Draw("COLZ");
    calorimeter->GetXaxis()->SetTitle("x");
    calorimeter->GetYaxis()->SetTitle("y");
    calorimeter->GetZaxis()->SetRangeUser(5E-4,1);
    cv->SetLogz(true);
    cv->SaveAs(imageName);
}
