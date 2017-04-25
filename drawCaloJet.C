#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"


void drawCaloJet(int jet, int nEtaBins, int nPhiBins) {

    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                nEtaBins,0,0.1*nEtaBins,
                                nPhiBins,0,0.1*nPhiBins); 

    std::ifstream inFile;
    std::string line;
    inFile.open("background.txt");
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
    st.SetPalette(kLightTemperature);            
    st.SetOptStat(0);
    cv->SetRightMargin(0.15);
    st.cd();            
    calorimeter->Draw("COLZ");
    calorimeter->GetZaxis()->SetRangeUser(0.1,1000);
    cv->SetLogz(true);
    cv->SaveAs("jet.png");
}

