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
#include "fastjet/tools/Filter.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"

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

void zeroCalorimeterAroundJet(TH2* calorimeter, int jetEtaBin, int jetPhiBin, int RinTowers, int nEtaBins, int nPhiBins)
{
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
                // Stop on eta
                if (iEta < 1 or iEta > nEtaBins) {content = 0;}
                // Loop on phi
                if (targetPhi < 1) {targetPhi += nPhiBins;}
                if (targetPhi > nPhiBins) {targetPhi -= nPhiBins;}

                /// Zero the calorimeter
                calorimeter->SetBinContent(iEta,targetPhi,0);
             }
        }
}

void fillCalorimeterWithJet(TH2* calorimeter, const fastjet::PseudoJet& jet)
{
    for (int i=0; i!=jet.constituents().size(); ++i) {
        calorimeter->Fill(jet.constituents().at(i).eta(),
                            jet.constituents().at(i).phi_std(),
                            jet.constituents().at(i).Et());
    }
} 

void fillCalorimeterWithParticles(TH2* calorimeter, const std::vector<TVector3>& particles)
{
    for (int i=0; i!=particles.size(); ++i) {
        calorimeter->Fill(particles.at(i).X(),
                            particles.at(i).Y(),
                            particles.at(i).Z());
    }
} 

void printConstituents(const fastjet::PseudoJet& jet)
{
    for (int i=0; i!=jet.constituents().size(); ++i) {
        cout << jet.constituents().at(i).eta() << " "
        << jet.constituents().at(i).phi_std() << " " 
        << jet.constituents().at(i).Et() << " " << endl;
    }
}

void saveCalorimeterImage(TH2* calorimeter, const char* name, double rangeMin = 0.1, double rangeMax = 1000)
{
        TCanvas cv("cv","cv",600,600);
        TStyle st;
        st.SetPalette(kLightTemperature);            
        st.SetOptStat(0);
        st.cd();            
        calorimeter->Draw("COLZ");
        if(rangeMin > 0 and rangeMax > 0)
            calorimeter->GetZaxis()->SetRangeUser(rangeMin,rangeMax);
        cv.SetLogz(true);
        cv.SetRightMargin(0.15);
        cv.SaveAs(name);
}

std::vector<TVector3> convertJetToParticles (const fastjet::PseudoJet& jet) 
{
    std::vector<TVector3> result;
    for(int i=0; i!= jet.constituents().size(); ++i) {
        result.push_back(TVector3(
            jet.constituents().at(i).eta(),
            jet.constituents().at(i).phi_std(),
            jet.constituents().at(i).Et()
        ));
    }
    return result;
}

void translateParticles(std::vector<TVector3>& particles, double x, double y) 
{
    /// Careful, still for this translation we need to make sure that phi is cyclical.
    /// After this we don't need to care anymore!
    for (int i=0; i!=particles.size(); ++i) {
        particles.at(i).SetXYZ(particles.at(i).X()+x,
                                TVector2::Phi_mpi_pi(particles.at(i).Y()+y),
                                particles.at(i).Z());
    }
    return;
}

void rotateParticles(std::vector<TVector3>& particles, double phi) 
{
    for (int i=0; i!=particles.size(); ++i) {
        particles.at(i).RotateZ(phi);
    }
    return;
}

void reflectParticlesIfNeeded(std::vector<TVector3>& particles) 
{
    double rightSideEt = 0;
    double leftSideEt = 0;
    for (int i=0; i!=particles.size(); ++i) {
        if(particles.at(i).X() >= 0) rightSideEt += particles.at(i).Z();
        else leftSideEt += particles.at(i).Z();
    }
    if(leftSideEt > rightSideEt) {
        for (int i=0; i!=particles.size(); ++i) {
            particles.at(i).SetX(-particles.at(i).X());
        }
    }
    return;
}

void normalizeParticles(std::vector<TVector3>& particles)
{
    double sumEtSquared = 0;
    double norm = 0;
    for (int i=0; i!=particles.size(); ++i) {
        sumEtSquared += particles.at(i).Z()*particles.at(i).Z();
    }
    //cout << "sumEtSquared = " << sumEtSquared << endl;
    for (int i=0; i!=particles.size(); ++i) {
        particles.at(i).SetZ(particles.at(i).Z()/sqrt(sumEtSquared));
        norm += particles.at(i).Z()*particles.at(i).Z();
    }
    //cout << "sumEtSquared = " << norm << endl;
}

int main()
{

    // Separator
    const std::string sep = "\t";

    // Number of events, generated and listed ones (for jets).
    int nEvent    = 22000;
    int nEventPU  = 30;
    int nListJets = 3;

    // LHC parameters
    double sqrtsInGeV = 13000.0;
    double meanPU = 30.0;
    
    // Select common parameters for SlowJet and FastJet analyses.
    int    power   = 0;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 1.2;    // Jet size.
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
    pythia.readString("PhaseSpace:pTHatMin = 150.");
    pythia.readString("PhaseSpace:pTHatMax = 300.");

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
    std::ofstream* nullStream = 0;
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
    //cout << "R in towers = " << RinTowers << endl;
    TH2D *caloJet = new TH2D("caloJet", "Calorimeter jet representation",
                                2*RinTowers+1,-0.05-R,R+0.05,
                                2*RinTowers+1,-0.05-R,R+0.05);

    // File
    std::ofstream outFile;
    outFile.open("background.txt");
    outFile.precision(6);
    outFile << std::scientific;

    /// Factoring out "fill calorimeter with event" function.

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
        calorimeter->Reset();
        caloJet->Reset();
        if (!pythia.next()) continue;

        /// Generate main event
        fillCalorimeterWithEvent(calorimeter,event,select,etaMax);
        
        nEventPU = rng.Poisson(meanPU);
        //cout << "Adding PU = " << nEventPU << endl;
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
                //cout << etaCenter << " " << phiCenter << endl;
            }
        }             
        //cout << fjInputs.size() << endl;
        
        /*
        fjInputs.clear();
        fastjet::PseudoJet jjj;
        jjj.reset_PtYPhiM(400,0.3,3.1+0.3-2*M_PI,0);
        fjInputs.push_back(jjj);
        jjj.reset_PtYPhiM(200,-0.3,3.1-0.3,0);
        fjInputs.push_back(jjj);
        */
        
        /// Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);

        nJets->Fill(sortedJets.size());

      
        /// Find leading jet
        double jetPt = sortedJets[0].perp();
        if (jetPt < 200) continue;
        if (jetPt > 250) continue;

        //cout << "jet mass = " << jetMass << endl;
        //cout << "Found pretrim constituents:" << sortedJets[0].constituents().size() << endl;   
             
        /// 1) Do the noise reduction with TRIMMING
        fastjet::Filter trimmer (fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), 
                                 fastjet::SelectorPtFractionMin(0.05));
        fastjet::PseudoJet trimmedJet = trimmer(sortedJets[0]);
        assert(trimmedJet.has_structure_of<fastjet::Filter>());

        double jetMass = trimmedJet.m();
        if(jetMass < 65) continue;  
        if(jetMass > 95) continue;  

        const fastjet::Filter::StructureType & fj_struct = trimmedJet.structure_of<fastjet::Filter>();
        
        /// 2) Doing the trimming automatically defines the points of interest as the subjets
        int nSubjets = trimmedJet.pieces().size();
        //cout << "Found subjets: " << nSubjets << endl;
        //cout << "Found posttrim constituents:" << trimmedJet.constituents().size() << endl;
        //printConstituents(sortedJets[0]);

        /// Found the position of the trimmed jet in the calorimeter, in order to zero it
        /// out and substitute for the trimmed constituents later
        double etaJet =  trimmedJet.eta();
        double phiJet =  trimmedJet.phi_std();
        double globalJetBin = calorimeter->FindFixBin(etaJet,phiJet);
        int jetEtaBin; int jetPhiBin; int jetZBin;
        calorimeter->GetBinXYZ(globalJetBin, jetEtaBin, jetPhiBin, jetZBin);
        //cout << "Jet centroid: " << etaJet << " " << phiJet << endl;
        for(int i=0; i!= nSubjets; ++i) {
            //cout << "Points of interest: " << trimmedJet.pieces().at(i).eta() << " " << trimmedJet.pieces().at(i).phi_std() << " " << trimmedJet.pieces().at(i).Et() << endl;
                                        }

        
        /// Convert this to a vector of TVector3, with X = eta, Y = phi_std, Z = Et
        std::vector<TVector3> particles = convertJetToParticles(trimmedJet);
  
        /// 3) Alignment

        /// Translation, such that the jet is centered in 0,0 in eta,phi
        translateParticles(particles,-etaJet,-phiJet);

        /// Rotation, such that if there are at least two subjets, the two more energetic ones are vertical
        double phiRotation = 0;
        if(nSubjets > 1) {
            /// Translate the two subjets as well to not make mistakes
            TVector2 v1(trimmedJet.pieces().at(0).eta()-etaJet,
                        TVector2::Phi_mpi_pi(trimmedJet.pieces().at(0).phi_std()-phiJet));
            TVector2 v2(trimmedJet.pieces().at(1).eta()-etaJet,
                        TVector2::Phi_mpi_pi(trimmedJet.pieces().at(1).phi_std()-phiJet));
            phiRotation = TVector2::Phi_mpi_pi((v2-v1).Phi()+M_PI_2);
        }
        rotateParticles(particles,-phiRotation);
        /// Reflection, such that the eta > 0 side has more sumEt than the eta < 0 side
        reflectParticlesIfNeeded(particles);
        
        /// 4) Normalization, such that sum of square of calo towers is one.
        normalizeParticles(particles);

        /// Fill a mini-calorimeter (just the jet) with the particles
        fillCalorimeterWithParticles(caloJet,particles);
        
        /// And write it to the outFile
        for(int iEta = 1; iEta <= caloJet->GetNbinsX(); ++iEta) {
            for(int iPhi= 1; iPhi <= caloJet->GetNbinsY(); ++iPhi) {
                outFile << caloJet->GetBinContent(iEta,iPhi) << sep;
            }
        }
        outFile << endl;
        
        /// Save nice plots
        if(nEvent==16) {
            saveCalorimeterImage(caloJet,"caloJet.png",1E-4,1);
            saveCalorimeterImage(calorimeter,"calorimeter.png");
            zeroCalorimeterAroundJet(calorimeter, jetEtaBin, jetPhiBin, RinTowers, nEtaBins, nPhiBins);
            saveCalorimeterImage(calorimeter,"calorimeter_hole.png");
            fillCalorimeterWithJet(calorimeter,trimmedJet);
            saveCalorimeterImage(calorimeter,"calorimeter_trimmed.png");
        }
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
