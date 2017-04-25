// Based on main72.cc program that is a part of the PYTHIA event generator.
// Author: Thiago Tomei

#include <iostream>
#include <fstream>
#include <string>
#include "Pythia8/Pythia.h"

// The FastJet3.h header enables automatic initialisation of
// fastjet::PseudoJet objects from Pythia8 Particle and Vec4 objects,
// as well as advanced features such as access to (a copy of)
// the original Pythia 8 Particle directly from the PseudoJet,
// and fastjet selectors that make use of the Particle properties.
// See the extensive comments in the header file for further details
// and examples.
#include "Pythia8Plugins/FastJet3.h"
//#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

using namespace Pythia8;

void fillCalorimeterWithEvent(TH2* calorimeter, const Pythia8::Event& event, int select, double etaMax)
{
    // Begin FastJet analysis: extract particles from event record.
        Vec4   pTemp;
        double mTemp;
        int nAnalyze = 0;
        for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

                // Require visible/charged particles inside detector.
                if (select > 2 &&  event[i].isNeutral()) continue;
                else if (select == 2 && !event[i].isVisible()) continue;
                if (etaMax < 20. && abs(event[i].eta()) > etaMax) continue;

                // Optionally modify mass and energy.
                pTemp = event[i].p();
                mTemp = event[i].m();

                /// Add the particle to the calorimeter
                calorimeter->Fill(pTemp.eta(),pTemp.phi(),pTemp.eT());
                ++nAnalyze;
        }
}

int main()
{

    // Separator
    const std::string sep = "\t";

    // Number of events, generated and listed ones (for jets).
    int nEvent    = 1;
    int nEventPU  = 0;
    int nListJets = 3;

    // LHC parameters
    double sqrtsInGeV = 13000.0;
    double meanPU = 30.0;
    
    // Select common parameters for SlowJet and FastJet analyses.
    int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 0.8;    // Jet size.
    double pTMin   = 30.0;    // Min jet pT.
    double etaMax  = 5.0;    // Pseudorapidity range of detector.
    int    select  = 2;      // Which particles are included?
    int    massSet = 2;      // Which mass are they assumed to have?

    // Generator. Shorthand for event.
    Pythia pythia("../share/Pythia8/xmldoc",false);
    Pythia pythiaPU("../share/Pythia8/xmldoc",false);

    Event &event = pythia.event;
    Event &eventPU = pythiaPU.event;

    TRandom3 rng(5713919);
    
    // Process selection.
    pythia.readString("HardQCD:all = on");
    pythiaPU.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("PhaseSpace:pTHatMin = 300.");

    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythiaPU.readString("Next:numberShowInfo = 0");
    pythiaPU.readString("Next:numberShowProcess = 0");
    pythiaPU.readString("Next:numberShowEvent = 0");

    // Very very quiet.
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythiaPU.readString("Init:showProcesses = off");
    pythiaPU.readString("Init:showMultipartonInteractions = off");
    std::ofstream* nullStream;
    fastjet::ClusterSequence::set_fastjet_banner_stream(nullStream);
    
    // LHC initialization.
    TString LHCinit = TString("Beams:eCM = ")+Form("%g",sqrtsInGeV);
    pythia.readString(LHCinit.Data());
    pythia.init();
    pythiaPU.readString(LHCinit.Data());
    pythiaPU.init();

    // Set up FastJet jet finder.
    //   one can use either explicitly use antikt, cambridge, etc., or
    //   just use genkt_algorithm with specification of power
    //fastjet::JetAlgorithm algorithm;
    //if (power == -1)      algorithm = fastjet::antikt_algorithm;
    //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
    //if (power ==  1)      algorithm = fastjet::kt_algorithm;
    //fastjet::JetDefinition jetDef(algorithm, R);
    // there's no need for a pointer to the jetDef (it's a fairly small object)
    fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
    std::vector <fastjet::PseudoJet> fjInputs;

    // Histograms.
    TH1F *nJets  = new TH1F("nJets", "number of jets", 15, -0.5, 14.5);
    TH1F *nParts = new TH1F("nParts", "number of particles", 200, -0.5, 199.5);
    int nEtaBins = 50;
    double maxCaloEta = 2.5;
    double minCaloEta = -2.5;
    int nPhiBins = 63;
    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                nEtaBins,minCaloEta,maxCaloEta,
                                nPhiBins,-3.14159,3.14159);
    int RinTowers = TMath::Nint(R/((maxCaloEta - minCaloEta)/nEtaBins));
    cout << "R in towers = " << RinTowers << endl;
    TCanvas* cv = new TCanvas("cv","cv",600,600);

    // File
    std::ofstream outFile;
    outFile.open("background.txt");
    outFile.precision(6);
    outFile << std::scientific;

    /// Factoring out "fill calorimeter with event" function.

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
        calorimeter->Reset();
        if (!pythia.next()) continue;

        /// Generate main event
        fillCalorimeterWithEvent(calorimeter,event,select,etaMax);
        
        nEventPU = rng.Poisson(meanPU);
        cout << "Adding PU = " << nEventPU << endl;
        /// Generate PU events
        for (int iEventPU = 0; iEventPU != nEventPU; ++iEventPU) {
            if (!pythiaPU.next()) continue;
            fillCalorimeterWithEvent(calorimeter,eventPU,select,etaMax);
        }   
                    
        /// Transform calo towers back into particles
        std::vector <fastjet::PseudoJet> fjInputs;
        for(int nBinsX=1 ; nBinsX<=calorimeter->GetNbinsX() ; ++nBinsX ) {
            for(int nBinsY=1 ; nBinsY<=calorimeter->GetNbinsY() ; ++nBinsY ) {
                fastjet::PseudoJet fji; 
                Double_t etaCenter = calorimeter->GetXaxis()->GetBinCenter(nBinsX);
                Double_t phiCenter = calorimeter->GetYaxis()->GetBinCenter(nBinsY);
                Double_t sumEt     = calorimeter->GetBinContent(nBinsX,nBinsY);
                fji.reset_PtYPhiM(sumEt,etaCenter,phiCenter);
                fjInputs.push_back(fji);
            }
        }             
        //cout << fjInputs.size() << endl;

        TStyle st;
        st.SetPalette(kLightTemperature);            
        st.SetOptStat(0);
        st.cd();            
        calorimeter->Draw("COLZ");
        calorimeter->GetZaxis()->SetRangeUser(0.1,1000);
        cv->SetLogz(true);
        //cv->SaveAs("calorimeter.png");

        /// Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);

        nJets->Fill(sortedJets.size());
        
        /// Find leading jet
        double jetPt = sortedJets[0].perp();
        if (jetPt < 300) continue;
        double etaJet =  sortedJets[0].eta();
        double phiJet =  sortedJets[0].phi_std();
        cout << etaJet << " " << phiJet << endl;
        double globalJetBin = calorimeter->FindFixBin(etaJet,phiJet);
        int jetEtaBin; int jetPhiBin; int jetZBin;
        calorimeter->GetBinXYZ(globalJetBin, jetEtaBin, jetPhiBin, jetZBin);
        cout << jetEtaBin << " " << jetPhiBin << endl;
        

        /// Print only the jet towers around the jet
        /// Now here we have to be careful
        /// If we are at a border in eta, we have to stop
        /// If we are at a border in phi, we have to loop
        int startingEta = jetEtaBin - RinTowers;
        int endingEta = jetEtaBin + RinTowers;
        int startingPhi = jetPhiBin - RinTowers;
        int endingPhi = jetPhiBin + RinTowers;
        for(int iEta = startingEta; iEta <= endingEta; ++ iEta) {
            for(int iPhi = startingPhi; iPhi <= endingPhi; ++ iPhi) {
                double content = 0;
                int targetPhi = iPhi;
                // Loop on eta
                if (targetPhi < 1) {targetPhi += nPhiBins;}
                if (targetPhi > nPhiBins) {targetPhi -= nPhiBins;}
                content = calorimeter->GetBinContent(iEta,targetPhi,0);
                // Stop on phi
                if (iEta < 1 or iEta > nEtaBins) {content = 0;}

                //cout << iEta << sep << iPhi << sep << content << endl;
                outFile << content << sep;
            }
        }
    outFile << endl;
    } //Close event loop 
       

    // Statistics. Histograms.
    //pythia.stat();
    outFile.close();

    TFile *f = TFile::Open("file.root", "RECREATE");
    nJets->Write();
    nParts->Write();
    f->Close();

    // Done.
    return 0;
}
