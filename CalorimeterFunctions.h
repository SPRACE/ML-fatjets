/// File: header1.h
/// Author: Thiago Tomei (SPRACE)
/// Goal: To group all functions that deal only with the internal representation
///		  of our calorimeter.

void zeroCalorimeterAroundJet(TH2 *calorimeter, int jetEtaBin, int jetPhiBin, int RinTowers, int nEtaBins, int nPhiBins)
{
    /// Now here we have to be careful
    /// If we are at a border in eta, we have to stop
    /// If we are at a border in phi, we have to loop
    int startingEta = jetEtaBin - RinTowers;
    int endingEta = jetEtaBin + RinTowers;
    int startingPhi = jetPhiBin - RinTowers;
    int endingPhi = jetPhiBin + RinTowers;
    for (int iEta = startingEta; iEta <= endingEta; ++ iEta) {
        for (int iPhi = startingPhi; iPhi <= endingPhi; ++ iPhi) {
            double content = 0;
            int targetPhi = iPhi;
            // Stop on eta
            if (iEta < 1 or iEta > nEtaBins) {
                content = 0;
            }
            // Loop on phi
            if (targetPhi < 1) {
                targetPhi += nPhiBins;
            }
            if (targetPhi > nPhiBins) {
                targetPhi -= nPhiBins;
            }

            /// Zero the calorimeter
            calorimeter->SetBinContent(iEta, targetPhi, 0);
        }
    }
}

void fillCalorimeterWithJet(TH2 *calorimeter, const fastjet::PseudoJet &jet)
{
    for (int i = 0; i != jet.constituents().size(); ++i) {
        calorimeter->Fill(jet.constituents().at(i).eta(),
                          jet.constituents().at(i).phi_std(),
                          jet.constituents().at(i).Et());
    }
}

void fillCalorimeterWithParticles(TH2 *calorimeter, const std::vector<TVector3> &particles)
{
    for (int i = 0; i != particles.size(); ++i) {
        calorimeter->Fill(particles.at(i).X(),
                          particles.at(i).Y(),
                          particles.at(i).Z());
    }
}

void printConstituents(const fastjet::PseudoJet &jet)
{
    for (int i = 0; i != jet.constituents().size(); ++i) {
        std::cout << jet.constituents().at(i).eta() << " "
             << jet.constituents().at(i).phi_std() << " "
             << jet.constituents().at(i).Et() << " " << std::endl;
    }
}

void saveCalorimeterImage(TH2 *calorimeter, const char *name, double rangeMin = 0.1, double rangeMax = 1000)
{
    TCanvas cv("cv", "cv", 600, 600);
    TStyle st;
    st.SetPalette(kLightTemperature);
    st.SetOptStat(0);
    st.cd();
    calorimeter->Draw("COLZ");
    if (rangeMin > 0 and rangeMax > 0)
        calorimeter->GetZaxis()->SetRangeUser(rangeMin, rangeMax);
    cv.SetLogz(true);
    cv.SetRightMargin(0.15);
    cv.SaveAs(name);
}

std::vector<TVector3> convertJetToParticles(const fastjet::PseudoJet &jet)
{
    std::vector<TVector3> result;
    for (int i = 0; i != jet.constituents().size(); ++i) {
        result.push_back(TVector3(
                             jet.constituents().at(i).eta(),
                             jet.constituents().at(i).phi_std(),
                             jet.constituents().at(i).Et()
                         ));
    }
    return result;
}

void translateParticles(std::vector<TVector3> &particles, double x, double y)
{
    /// Careful, still for this translation we need to make sure that phi is cyclical.
    /// After this we don't need to care anymore!
    for (int i = 0; i != particles.size(); ++i) {
        particles.at(i).SetXYZ(particles.at(i).X() + x,
                               TVector2::Phi_mpi_pi(particles.at(i).Y() + y),
                               particles.at(i).Z());
    }
    return;
}

void rotateParticles(std::vector<TVector3> &particles, double phi)
{
    for (int i = 0; i != particles.size(); ++i) {
        particles.at(i).RotateZ(phi);
    }
    return;
}

void reflectParticlesIfNeeded(std::vector<TVector3> &particles)
{
    double rightSideEt = 0;
    double leftSideEt = 0;
    for (int i = 0; i != particles.size(); ++i) {
        if (particles.at(i).X() >= 0) rightSideEt += particles.at(i).Z();
        else leftSideEt += particles.at(i).Z();
    }
    if (leftSideEt > rightSideEt) {
        for (int i = 0; i != particles.size(); ++i) {
            particles.at(i).SetX(-particles.at(i).X());
        }
    }
    return;
}

void normalizeParticles(std::vector<TVector3> &particles)
{
    double sumEtSquared = 0;
    double norm = 0;
    for (int i = 0; i != particles.size(); ++i) {
        sumEtSquared += particles.at(i).Z() * particles.at(i).Z();
    }
    //std::cout << "sumEtSquared = " << sumEtSquared << std::endl;
    for (int i = 0; i != particles.size(); ++i) {
        particles.at(i).SetZ(particles.at(i).Z() / sqrt(sumEtSquared));
        norm += particles.at(i).Z() * particles.at(i).Z();
    }
    //std::cout << "sumEtSquared = " << norm << std::endl;
}