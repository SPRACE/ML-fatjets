#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

double calculateTau21(int nEtaBins, int nPhiBins, double jetRadius,
						const char *fileName, TH1* histo)
{
	using namespace std;

    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                 nEtaBins, -jetRadius, jetRadius,
                                 nPhiBins, -jetRadius, jetRadius);

    double etaBinSize = jetRadius / ((nEtaBins - 1) / 2.0);
    double phiBinSize = jetRadius / ((nPhiBins - 1) / 2.0);
    double startingEta = -jetRadius + etaBinSize / 2;
    double startingPhi = -jetRadius + phiBinSize / 2;

    std::ifstream inFile;
    std::string line;
    inFile.open(fileName);
    while (std::getline(inFile, line)) {
        
        std::istringstream iss(line);
        //std::cout << line <<  std::endl;
        fastjet::PseudoJet fji;
        std::vector <fastjet::PseudoJet> fjInputs;
        for (int iEta = 1; iEta <= nEtaBins; ++iEta) {
            for (int iPhi = 1; iPhi <= nPhiBins; ++iPhi) {
                Double_t content = 0;
                Double_t etaCenter = calorimeter->GetXaxis()->GetBinCenter(iEta);
                Double_t phiCenter = calorimeter->GetYaxis()->GetBinCenter(iPhi);
                iss >> content;
                fji.reset_PtYPhiM(content, etaCenter, phiCenter);
                fjInputs.push_back(fji);
            }
        }

        /// Run Fastjet algorithm and sort jets in pT order.
        double pTMin = 0.0;
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetRadius);
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);
        fastjet::PseudoJet this_jet = sortedJets.at(0);

        /// N-subjettiness
        double beta = 1;
        fastjet::contrib::NsubjettinessRatio nSub21_beta1(2, 1, fastjet::contrib::OnePass_KT_Axes(),
                fastjet::contrib::NormalizedMeasure(beta, jetRadius));
        double tau21_beta1 = nSub21_beta1(this_jet);

        histo->Fill(tau21_beta1);
    }
    return 1;
}


int main()
{
	std::string outRootFileName = "output.root";
   	TFile* outRootFile = TFile::Open(outRootFileName.c_str(),"RECREATE");
   	TH1D* histoTau21 = new TH1D("tau21", "tau21", 50, 0, 1);
    calculateTau21(25, 25, 1.2, "background_PU0_13TeV_MJ-65-95_PTJ-250-300_ext.txt", histoTau21);
	histoTau21->Write();
	outRootFile->Close();
	
    return 0;
}
