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

using namespace Pythia8;

int main()
{

    // Separator
    const std::string sep = "\t";

    // Number of events, generated and listed ones (for jets).
    int nEvent    = 49999;
    int nListJets = 3;

    // Select common parameters for SlowJet and FastJet analyses.
    int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 0.8;    // Jet size.
    double pTMin   = 30.0;    // Min jet pT.
    double etaMax  = 5.0;    // Pseudorapidity range of detector.
    int    select  = 2;      // Which particles are included?
    int    massSet = 2;      // Which mass are they assumed to have?

    // Generator. Shorthand for event.
    Pythia pythia;
    Event &event = pythia.event;
    Event &process = pythia.process;

    // Process selection.
    //pythia.readString("HardQCD:all = on");
    //pythia.readString("PhaseSpace:pTHatMin = 300.");

    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = cmsgrid_final.lhe.gz");
    pythia.init();

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

    // File
    std::ofstream outFile;
    outFile.open("signal.txt");
    outFile.precision(6);
    outFile << std::scientific;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;

        if (iEvent < 10) process.list();

        // Find the hadronic Z, else continue.
        int ZhadIndex = -1;
        for (int p = 0; p != process.size(); ++p) {
            if (process.at(p).id() == 23 and process.at(process.at(p).daughter1()).isQuark())
                ZhadIndex = p;
        }
        if (ZhadIndex < 0) continue;
        Particle Zhad = process.at(ZhadIndex);
        fastjet::PseudoJet fjZhad(Zhad.px(), Zhad.py(), Zhad.pz(), Zhad.e());

        // Begin FastJet analysis: extract particles from event record.
        fjInputs.resize(0);
        Vec4   pTemp;
        double mTemp;
        int nAnalyze = 0;
        for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

                // Require visible/charged particles inside detector.
                if (select > 2 &&  event[i].isNeutral()) continue;
                else if (select == 2 && !event[i].isVisible()) continue;
                if (etaMax < 20. && abs(event[i].eta()) > etaMax) continue;

                // Create a PseudoJet from the complete Pythia particle.
                fastjet::PseudoJet particleTemp = event[i];

                // Optionally modify mass and energy.
                pTemp = event[i].p();
                mTemp = event[i].m();
                if (massSet < 2) {
                    mTemp = (massSet == 0 || event[i].id() == 22) ? 0. : 0.13957;
                    pTemp.e(sqrt(pTemp.pAbs2() + mTemp * mTemp));
                    particleTemp.reset_momentum(pTemp.px(), pTemp.py(),
                                                pTemp.pz(), pTemp.e());
                }

                // Store acceptable particles as input to Fastjet.
                // Conversion to PseudoJet is performed automatically
                // with the help of the code in FastJet3.h.
                fjInputs.push_back(particleTemp);
                ++nAnalyze;
            }

        // Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);

        nJets->Fill(sortedJets.size());

        // Loop over jets, find the matched one
        double minDR = 9999.9;
        int matchedJet = -1;
        for (size_t ii = 0; ii != sortedJets.size(); ++ii) {
            double thisDR = sqrt(fjZhad.squared_distance(sortedJets[ii]));
            if (thisDR < minDR) {
                minDR = thisDR;
                matchedJet = ii;
            }
        }
        if (minDR > 0.4 or matchedJet < 0) continue;

        nParts->Fill(sortedJets[matchedJet].constituents().size());

        double jetPt = sortedJets[matchedJet].perp();
        if (jetPt > 300) {
            outFile << "1" << sep << jetPt << sep;
            size_t nParticles = sortedJets[matchedJet].constituents().size();
            // Loop over particles, but not more than 200.
            for (size_t kk = 0; kk != nParticles and kk != 200; ++kk) {
                double px =  sortedJets[matchedJet].constituents().at(kk).px();
                double py =  sortedJets[matchedJet].constituents().at(kk).py();
                double pz =  sortedJets[matchedJet].constituents().at(kk).pz();
                outFile << px / jetPt << sep << py / jetPt << sep << pz / jetPt << sep;
            }
            // Pad with zeros
            size_t nPads = (nParticles < 200) ? (200 - nParticles) : 0;
            for (size_t kk = 0; kk != nPads; ++kk) {
                outFile << 0 << sep << 0 << sep << 0 << sep;
            }
            outFile << std::endl;
        }


        // List first few FastJet jets and some info about them.
        // Note: the final few columns are illustrative of what information
        // can be extracted, but does not exhaust the possibilities.
        if (iEvent < nListJets) {
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
        }

        // End of event loop.
    }

    // Statistics. Histograms.
    outFile.close();
    pythia.stat();

    TFile *f = TFile::Open("file.root", "RECREATE");
    nJets->Write();
    nParts->Write();
    f->Close();

    // Done.
    return 0;
}
