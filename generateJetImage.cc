// Based on main72.cc program that is a part of the PYTHIA event generator.
// Author: Thiago Tomei

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "optionparser.h"

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

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TString.h"

// Delphes stuff
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "external/ExRootAnalysis/ExRootUtilities.h"

#include "CalorimeterFunctions.h"
#include "CalorimeterFillers.h"
 
int main(int argc, char *argv[])
{
    using namespace Pythia8;

    struct Arg: public option::Arg {
        static void printError(const char *msg1, const option::Option &opt, const char *msg2)
        {
            fprintf(stderr, "ERROR: %s", msg1);
            fwrite(opt.name, opt.namelen, 1, stderr);
            fprintf(stderr, "%s", msg2);
        }
        static option::ArgStatus Unknown(const option::Option &option, bool msg)
        {
            if (msg) printError("Unknown option '", option, "'\n");
            return option::ARG_ILLEGAL;
        }
        static option::ArgStatus Required(const option::Option &option, bool msg)
        {
            if (option.arg != 0)
                return option::ARG_OK;
            if (msg) printError("Option '", option, "' requires an argument\n");
            return option::ARG_ILLEGAL;
        }
        static option::ArgStatus NonEmpty(const option::Option &option, bool msg)
        {
            if (option.arg != 0 && option.arg[0] != 0)
                return option::ARG_OK;
            if (msg) printError("Option '", option, "' requires a non-empty argument\n");
            return option::ARG_ILLEGAL;
        }
        static option::ArgStatus Numeric(const option::Option &option, bool msg)
        {
            char *endptr = 0;
            if (option.arg != 0 && strtol(option.arg, &endptr, 10)) {};
            if (endptr != option.arg && *endptr == 0)
                return option::ARG_OK;
            if (msg) printError("Option '", option, "' requires a numeric argument\n");
            return option::ARG_ILLEGAL;
        }
    };

    enum  optionIndex { UNKNOWN, OUTFNAME, NUMEVENTS, AVERAGEPU, SQRTS,
                        MINJMASS, MAXJMASS, MINJETPT, MAXJETPT,
                        PROCESS
                      };
    const option::Descriptor usage[] = {
        {
            UNKNOWN,  0, "", "", option::Arg::None, "USAGE: generateJetImage [options]\n\n"
            "Options:"
        },
        {OUTFNAME, 0, "n", "outFileName", Arg::NonEmpty, " --outFileName \tOutput file name" },
        {NUMEVENTS, 0, "n", "numEvents", Arg::NonEmpty, " --numEvents, -n \tNumber of tried events" },
        {AVERAGEPU, 0, "", "averagePU", Arg::NonEmpty, " --averagePU \tAverage pileup" },
        {SQRTS,    0, "", "sqrts", Arg::NonEmpty, " --sqrts \tCenter of mass energy in GeV"},
        {MINJMASS, 0, "", "minJetMass", Arg::NonEmpty, " --minJetMass \tMin jet mass in GeV"},
        {MAXJMASS, 0, "", "maxJetMass", Arg::NonEmpty, " --maxJetMass \tMax jet mass in GeV"},
        {MINJETPT, 0, "", "minJetPt", Arg::NonEmpty, " --minJetPt\tMin jet mass in GeV"},
        {MAXJETPT, 0, "", "maxJetPt", Arg::NonEmpty, " --maxJetPt\tMax jet mass in GeV"},
        {PROCESS,  0, "", "process", Arg::NonEmpty, " --process\tProcess: 0=bkg, 1=sig"},
        {0, 0, 0, 0, 0, 0}
    };

    /// Parsing the options
    argc -= (argc > 0);
    argv += (argc > 0); // skip program name argv[0] if present
    option::Stats  stats(usage, argc, argv);
    option::Option options[stats.options_max], buffer[stats.buffer_max];
    option::Parser parse(usage, argc, argv, options, buffer);

    if (parse.error()) return 1;
    if (options[UNKNOWN]) {
        option::printUsage(std::cout, usage);
        return 0;
    }

    // Separator
    const std::string sep = "\t";

    // Basic run parameters
    int nEvent    = 1000;
    int nEventPU  = 0;
    int nListJets = 3;
    int processType =  0; //0 = background, 1 = signal
    string outFileName = "output.txt";

    // LHC parameters
    double sqrtsInGeV = 13000.0;
    double meanPU = 0.0;

    double _minJetMass = 0;
    double _maxJetMass = 9999;
    double _minJetPt = 0;
    double _maxJetPt = 9999;


    if (options[OUTFNAME]) outFileName = options[OUTFNAME].arg;
    if (options[NUMEVENTS]) nEvent = atof(options[NUMEVENTS].arg);
    if (options[AVERAGEPU]) meanPU = atof(options[AVERAGEPU].arg);
    if (options[SQRTS]) sqrtsInGeV = atof(options[SQRTS].arg);
    if (options[MINJMASS]) _minJetMass = atof(options[MINJMASS].arg);
    if (options[MAXJMASS]) _maxJetMass = atof(options[MAXJMASS].arg);
    if (options[MINJETPT]) _minJetPt = atof(options[MINJETPT].arg);
    if (options[MAXJETPT]) _maxJetPt = atof(options[MAXJETPT].arg);
    if (options[PROCESS]) processType = atoi(options[PROCESS].arg);

    // Select jet finding parameters.
    int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 0.8;    // Jet size.
    double pTMin   = 30.0;    // Min jet pT.
    double etaMax  = 5.0;    // Pseudorapidity range of detector.
    int    select  = 2;      // Which particles are included?
    int    massSet = 2;      // Which mass are they assumed to have?

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
    int nEtaBins = 34;
    double maxCaloEta = 1.479;
    double minCaloEta = -1.479;
    int nPhiBins = 72;
    TH2D *calorimeter = new TH2D("calorimeter", "Calorimeter representation",
                                 nEtaBins, minCaloEta, maxCaloEta,
                                 nPhiBins, -3.141593, 3.141593);
    double etaStep = calorimeter->GetXaxis()->GetBinWidth(1);
    double phiStep = calorimeter->GetYaxis()->GetBinWidth(1);
    int RinTowers = TMath::Nint(R / ((maxCaloEta - minCaloEta) / nEtaBins));
    if(RinTowers%2 == 0) RinTowers = (RinTowers - 1);
//     std::cout << "R in towers = " << RinTowers << std::endl;
    TH2D *caloJet = new TH2D("caloJet", "Calorimeter jet representation",
                             2 * RinTowers + 1, -etaStep/2 - R, R + etaStep/2,
                             2 * RinTowers + 1, -phiStep/2 - R, R + phiStep/2);

    // Files
    std::ofstream outFile;
    outFile.open(outFileName.c_str());
    outFile.precision(6);
    outFile << std::scientific;

    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add("~/Desktop/tag_1_delphes_events341.root");

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    
    // Get pointers to branches used in this analysis
    TClonesArray *branchTower = treeReader->UseBranch("Tower");

    /*TFile* inputFile = TFile::Open("~/Desktop/tag_1_delphes_events341.root");
    TTree* tree = (TTree*) inputFile->Get("Delphes");
    TTree          *fChain;
    Int_t           fCurrent; //!current Tree number in a TChain
    const static Int_t kMaxEvent = 9999;
    const static Int_t kMaxTower = 50000;
    Int_t           Tower_;
    Float_t         Tower_Eta[kMaxTower];
    Float_t         Tower_Phi[kMaxTower];
    Float_t         Tower_Eem[kMaxTower];
    Float_t         Tower_Ehad[kMaxTower];
    Float_t         Tower_Edges[kMaxTower][4];
    Int_t           Tower_size;
    TBranch        *b_Tower_;   //!
    TBranch        *b_Tower_Eta;   //!
    TBranch        *b_Tower_Phi;   //!
    TBranch        *b_Tower_Eem;   //!
    TBranch        *b_Tower_Ehad;   //!
    TBranch        *b_Tower_Edges;   //!
    TBranch        *b_Tower_size;   //!
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("Tower", &Tower_, &b_Tower_);
    fChain->SetBranchAddress("Tower.Eta", Tower_Eta, &b_Tower_Eta);
    fChain->SetBranchAddress("Tower.Phi", Tower_Phi, &b_Tower_Phi);
    fChain->SetBranchAddress("Tower.Eem", Tower_Eem, &b_Tower_Eem);
    fChain->SetBranchAddress("Tower.Ehad", Tower_Ehad, &b_Tower_Ehad);
    fChain->SetBranchAddress("Tower.Edges[4]", Tower_Edges, &b_Tower_Edges);
    fChain->SetBranchAddress("Tower_size", &Tower_size, &b_Tower_size);
    */
    
    /// Factoring out "fill calorimeter with event" function.

    // Begin event loop. Generate event. Skip if error.
    Long64_t nentries = chain.GetEntriesFast();
    nentries = 12;
    for (Long64_t jentry=8; jentry<nentries;jentry++) {

        treeReader->ReadEntry(jentry);        
        calorimeter->Reset();
        caloJet->Reset();

        /// Fill calorimeter with it
        fillCalorimeterWithDelphesCalorimeter(calorimeter, 
                                                branchTower,
                                                1.479);

        /// Transform calo towers back into particles
        std::vector <fastjet::PseudoJet> fjInputs;
        for (int nBinsX = 1 ; nBinsX <= calorimeter->GetNbinsX() ; ++nBinsX) {
            for (int nBinsY = 1 ; nBinsY <= calorimeter->GetNbinsY() ; ++nBinsY) {
                fastjet::PseudoJet fji;
                Double_t etaCenter = calorimeter->GetXaxis()->GetBinCenter(nBinsX);
                Double_t phiCenter = calorimeter->GetYaxis()->GetBinCenter(nBinsY);
                Double_t sumEt     = calorimeter->GetBinContent(nBinsX, nBinsY);
                fji.reset_PtYPhiM(sumEt, etaCenter, phiCenter);
                fjInputs.push_back(fji);
                //std::cout << etaCenter << " " << phiCenter << std::endl;
            }
        }

        /// Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);

        nJets->Fill(sortedJets.size());

        if (sortedJets.size() == 0) continue;

        /// Find leading jet
        double jetPt = sortedJets[0].perp();
        double jetEta = sortedJets[0].eta();
        double jetPhi = sortedJets[0].phi_std();
        if (jetPt < _minJetPt) continue;
        if (jetPt > _maxJetPt) continue;

        //std::cout << "jet pt = " << jetPt << std::endl;
        //std::cout << "jet eta = " << jetEta << std::endl;
        //std::cout << "jet phi = " << jetPhi << std::endl;
        //std::cout << "Found pretrim constituents:" << sortedJets[0].constituents().size() << std::endl;

        /// 1) Do the noise reduction with TRIMMING
        fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
                                fastjet::SelectorPtFractionMin(0.05));
        fastjet::PseudoJet trimmedJet = trimmer(sortedJets[0]);
        assert(trimmedJet.has_structure_of<fastjet::Filter>());

        double jetMass = trimmedJet.m();
        if (jetMass < _minJetMass) continue;
        if (jetMass > _maxJetMass) continue;

        const fastjet::Filter::StructureType &fj_struct = trimmedJet.structure_of<fastjet::Filter>();

        /// 2) Doing the trimming automatically defines the points of interest as the subjets
        int nSubjets = trimmedJet.pieces().size();
        //std::cout << "Found subjets: " << nSubjets << std::endl;
        //std::cout << "Found posttrim constituents:" << trimmedJet.constituents().size() << std::endl;
        //printConstituents(sortedJets[0]);

        /// Found the position of the trimmed jet in the calorimeter, in order to zero it
        /// out and substitute for the trimmed constituents later
        double etaJet =  trimmedJet.eta();
        double phiJet =  trimmedJet.phi_std();
        double globalJetBin = calorimeter->FindFixBin(etaJet, phiJet);
        int jetEtaBin;
        int jetPhiBin;
        int jetZBin;
        calorimeter->GetBinXYZ(globalJetBin, jetEtaBin, jetPhiBin, jetZBin);
        //std::cout << "Jet centroid: " << etaJet << " " << phiJet << std::endl;
        for (int i = 0; i != nSubjets; ++i) {
            //std::cout << "Points of interest: " << trimmedJet.pieces().at(i).eta() << " " << trimmedJet.pieces().at(i).phi_std() << " " << trimmedJet.pieces().at(i).Et() << std::endl;
        }


        /// Convert this to a vector of TVector3, with X = eta, Y = phi_std, Z = Et
        std::vector<TVector3> particles = convertJetToParticles(trimmedJet);

        /// 3) Alignment

        /// Translation, such that the jet is centered in 0,0 in eta,phi
        translateParticles(particles, -etaJet, -phiJet);

        /// Rotation, such that if there are at least two subjets, the two more energetic ones are vertical
        double phiRotation = 0;
        if (nSubjets > 1) {
            /// Translate the two subjets as well to not make mistakes
            TVector2 v1(trimmedJet.pieces().at(0).eta() - etaJet,
                        TVector2::Phi_mpi_pi(trimmedJet.pieces().at(0).phi_std() - phiJet));
            TVector2 v2(trimmedJet.pieces().at(1).eta() - etaJet,
                        TVector2::Phi_mpi_pi(trimmedJet.pieces().at(1).phi_std() - phiJet));
            phiRotation = TVector2::Phi_mpi_pi((v2 - v1).Phi() + M_PI_2);
        }
        rotateParticles(particles, -phiRotation);
        /// Reflection, such that the eta > 0 side has more sumEt than the eta < 0 side
        reflectParticlesIfNeeded(particles);

        /// 4) Normalization, such that sum of square of calo towers is one.
        normalizeParticles(particles);

        /// Fill a mini-calorimeter (just the jet) with the particles
        fillCalorimeterWithParticles(caloJet, particles);

        /// And write it to the outFile
        for (int iEta = 1; iEta <= caloJet->GetNbinsX(); ++iEta) {
            for (int iPhi = 1; iPhi <= caloJet->GetNbinsY(); ++iPhi) {
                outFile << caloJet->GetBinContent(iEta, iPhi) << sep;
            }
        }
        outFile << std::endl;

        //std::cout << "iEvent == " << iEvent << std::endl;
        /// Save nice plots
        if (jentry == 11) {
            std::cout << "Printing event" << std::endl;
            saveCalorimeterImage(caloJet, "caloJet.png", 1E-4, 1);
            saveCalorimeterImage(calorimeter, "calorimeter.png");
            zeroCalorimeterAroundJet(calorimeter, jetEtaBin, jetPhiBin, RinTowers, nEtaBins, nPhiBins);
            saveCalorimeterImage(calorimeter, "calorimeter_hole.png");
            fillCalorimeterWithJet(calorimeter, trimmedJet);
            saveCalorimeterImage(calorimeter, "calorimeter_trimmed.png");
        }
        
    } //Close event loop


    // Statistics. Histograms.
    //pythia.stat();
    outFile.close();

    TFile *f = TFile::Open("file.root", "RECREATE");
    nJets->Write();
    nParts->Write();
    calorimeter->Write();
    caloJet->Write();
    f->Close();

    // Done.
    return 0;
}
