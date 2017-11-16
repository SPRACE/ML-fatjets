#include "Pythia8/Pythia.h"

void fillCalorimeterWithEvent(TH2 *calorimeter, const Pythia8::Event &event, int select, double etaMax)
{
    Pythia8::Vec4   pTemp;
    int nAnalyze = 0; // Not used for the time being

    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

            // Require visible/charged particles inside detector.
            if (select > 2 &&  event[i].isNeutral()) continue;
            else if (select == 2 && !event[i].isVisible()) continue;
            //if (select == 2 && abs(event[i].id())==13) continue; // No muons
            if (etaMax < 20. && std::abs(event[i].eta()) > etaMax) continue;

            // Optionally modify mass and energy.
            pTemp = event[i].p();
            //mTemp = event[i].m();

            /// Add the particle to the calorimeter
            calorimeter->Fill(pTemp.eta(), pTemp.phi(), pTemp.eT());
            ++nAnalyze;
        }
}

const static Int_t kMaxTower = 50000;
void fillCalorimeterWithDelphesCalorimeter(TH2 *calorimeter, 
	const TClonesArray* towers,
	double etaMax)
{
// 	std::cout << "Found towers: " << towers->GetEntries() << std::endl;
    for (int iTower=0; iTower!= towers->GetEntries(); ++iTower) {
    		Tower *theTower = (Tower*) towers->At(iTower);
//      		std::cout << "iTower == " << iTower << std::endl;
//   			std::cout << "eta,phi = " << theTower->Eta
//   										<< " " << theTower->Phi  << "\t\t";
//   			std::cout << "Edges = " << theTower->Edges[0] << " " 
//   									<< theTower->Edges[1] << " "
//   									<< theTower->Edges[2] << " "
//   									<< theTower->Edges[3] << std::endl;
//   			std::cout << "Eem = " << theTower->Eem 
//   						<< ", Ehad = " << theTower->Ehad << std::endl;
			double towerEta = theTower->Eta;
			if(std::abs(towerEta) > etaMax) continue;
			double towerPhi = theTower->Phi;
			double towerLoEta = theTower->Edges[0];
			double towerHiEta = theTower->Edges[1];
			double towerLoPhi = theTower->Edges[2];
			double towerHiPhi = theTower->Edges[3];
			double towerEem = theTower->Eem;
			double towerEhad = theTower->Ehad;
			double towerEt = theTower->ET;
 			double energy = 0;

 			/// Calculate centre-of-tower coordinates
 			double eta =  (towerLoEta+towerHiEta)/2.0;
 			TVector2 v1,v2; v1.SetMagPhi(1,towerLoPhi); v2.SetMagPhi(1,towerHiPhi);
 			double phi =  TVector2::Phi_mpi_pi((v1+v2).Phi());
 			if(towerEhad > 0 and towerEem > 0) std::cout << "WARNING::Tower eith both Ehad and Eem" << std::endl;
			if(towerEhad > 0) { // This is a hadronic tower
				energy = towerEhad;
// 				calorimeter->Fill(eta-2*0.0174,phi-2*0.0174,energy);
// 				calorimeter->Fill(eta-1*0.0174,phi-2*0.0174,energy);
// 				calorimeter->Fill(eta-0*0.0174,phi-2*0.0174,energy);
// 				calorimeter->Fill(eta+1*0.0174,phi-2*0.0174,energy);
// 				calorimeter->Fill(eta+2*0.0174,phi-2*0.0174,energy);
// 				
// 				calorimeter->Fill(eta-2*0.0174,phi-1*0.0174,energy);
// 				calorimeter->Fill(eta-1*0.0174,phi-1*0.0174,energy);
// 				calorimeter->Fill(eta-0*0.0174,phi-1*0.0174,energy);
// 				calorimeter->Fill(eta+1*0.0174,phi-1*0.0174,energy);
// 				calorimeter->Fill(eta+2*0.0174,phi-1*0.0174,energy);
// 				
// 				calorimeter->Fill(eta-2*0.0174,phi-0*0.0174,energy);
// 				calorimeter->Fill(eta-1*0.0174,phi-0*0.0174,energy);
 				calorimeter->Fill(eta-0*0.0174,phi-0*0.0174,towerEt);
// 				calorimeter->Fill(eta+1*0.0174,phi-0*0.0174,energy);
// 				calorimeter->Fill(eta+2*0.0174,phi-0*0.0174,energy);
// 				
// 				calorimeter->Fill(eta-2*0.0174,phi+1*0.0174,energy);
// 				calorimeter->Fill(eta-1*0.0174,phi+1*0.0174,energy);
// 				calorimeter->Fill(eta-0*0.0174,phi+1*0.0174,energy);
// 				calorimeter->Fill(eta+1*0.0174,phi+1*0.0174,energy);
// 				calorimeter->Fill(eta+2*0.0174,phi+1*0.0174,energy);
// 				
// 				calorimeter->Fill(eta-2*0.0174,phi+2*0.0174,energy);
// 				calorimeter->Fill(eta-1*0.0174,phi+2*0.0174,energy);
// 				calorimeter->Fill(eta-0*0.0174,phi+2*0.0174,energy);
// 				calorimeter->Fill(eta+1*0.0174,phi+2*0.0174,energy);
// 				calorimeter->Fill(eta+2*0.0174,phi+2*0.0174,energy);
			}
			else { // This is an electromagnetic tower tower
				energy = towerEem;
				calorimeter->Fill(eta,phi,towerEt);
			}	
		}
}


