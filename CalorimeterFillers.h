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
