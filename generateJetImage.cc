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

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace Pythia8;

int main()
{

    // Separator
    const std::string sep = "\t";

    // Number of events, generated and listed ones (for jets).
    int nEvent    = 1;
    int nListJets = 3;

    // Select common parameters for SlowJet and FastJet analyses.
    int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 0.4;    // Jet size.
    double pTMin   = 30.0;    // Min jet pT.
    double etaMax  = 5.0;    // Pseudorapidity range of detector.
    int    select  = 2;      // Which particles are included?
    int    massSet = 2;      // Which mass are they assumed to have?

    // Generator. Shorthand for event.
    Pythia pythia;
    Event &event = pythia.event;

    // Process selection.
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 300.");

    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    // LHC initialization.
    pythia.readString("Beams:eCM = 13000.");
    pythia.init();

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
    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                50,-2.5,2.5,63,-3.14159,3.14159);
    TCanvas* cv = new TCanvas("cv","cv",600,600);

    // File
    std::ofstream outFile;
    outFile.open("background.txt");
    outFile.precision(6);
    outFile << std::scientific;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;


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
        cv->SaveAs("calorimeter.png");

        /// Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);

        nJets->Fill(sortedJets.size());
        
        double jetPt = sortedJets[0].perp();
        if (jetPt < 300) continue;
        double etaJet =  sortedJets[0].eta();
        double phiJet =  sortedJets[0].phi_std();
        cout << etaJet << " " << phiJet << endl;
        double globalJetBin = calorimeter->FindFixBin(etaJet,phiJet);
        int jetEtaBin; int jetPhiBin; int jetZBin;
        calorimeter->GetBinXYZ(globalJetBin, jetEtaBin, jetPhiBin, jetZBin);
        cout << jetEtaBin << " " << jetPhiBin << endl;
        
        /// Now here we have to be careful
        /// If we are at a border in eta, we have to stop
        /// If we are at a border in phi, we have to loop

/*            

        for (size_t ii = 0; ii != sortedJets.size(); ++ii) {
            nParts->Fill(sortedJets[ii].constituents().size());

            double jetPt = sortedJets[ii].perp();
            if (jetPt > 300) {
                outFile << "0" << sep << jetPt << sep;
                size_t nParticles = sortedJets[ii].constituents().size();
                // Loop over particles, but not more than 200.
                for (size_t kk = 0; kk != nParticles and kk != 200; ++kk) {
                    double px =  sortedJets[ii].constituents().at(kk).px();
                    double py =  sortedJets[ii].constituents().at(kk).py();
                    double pz =  sortedJets[ii].constituents().at(kk).pz();
                    outFile << px / jetPt << sep << py / jetPt << sep << pz / jetPt << sep;
                }
                // Pad with zeros
                size_t nPads = (nParticles < 200) ? (200 - nParticles) : 0;
                for (size_t kk = 0; kk != nPads; ++kk) {
                    outFile << 0 << sep << 0 << sep << 0 << sep;
                }
                outFile << std::endl;
            }
*/
        }
        // List first few FastJet jets and some info about them.
        // Note: the final few columns are illustrative of what information
        // can be extracted, but does not exhaust the possibilities.
        /*if (iEvent < nListJets) {
            cout << "\n --------  FastJet jets, p = " << setw(2) << power
                 << "  --------------------------------------------------\n\n "
                 << "  i         pT        y      phi  mult chgmult photons"
                 << "      hardest  pT in neutral " << endl
                 << "                                                       "
                 << "  constituent        hadrons " << endl;
            for (int i = 0; i < int(sortedJets.size()); ++i) {
                vector<fastjet::PseudoJet> constituents
                    = sortedJets[i].constituents();
                fastjet::PseudoJet hardest
                    = fastjet::SelectorNHardest(1)(constituents)[0];
                vector<fastjet::PseudoJet> neutral_hadrons
                    = (fastjet::SelectorIsHadron()
                       && fastjet::SelectorIsNeutral())(constituents);
                double neutral_hadrons_pt = join(neutral_hadrons).perp();
                cout << setw(4) << i << fixed << setprecision(3) << setw(11)
                     << sortedJets[i].perp() << setw(9)  << sortedJets[i].rap()
                     << setw(9) << sortedJets[i].phi_std()
                     << setw(6) << constituents.size()
                     << setw(8) << fastjet::SelectorIsCharged().count(constituents)
                     << setw(8) << fastjet::SelectorId(22).count(constituents)
                     << setw(13) << hardest.user_info<Particle>().name()
                     << "     " << setw(10) << neutral_hadrons_pt << endl;
            }
            cout << "\n --------  End FastJet Listing  ------------------"
                 << "---------------------------------" << endl;
        }*/

        // End of event loop.
    

    // Statistics. Histograms.
    pythia.stat();

    TFile *f = TFile::Open("file.root", "RECREATE");
    nJets->Write();
    nParts->Write();
    f->Close();

    // Done.
    return 0;
}
