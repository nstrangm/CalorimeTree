#include "TMath.h"
#include "Cuts.h"
#include "Utilities.h"
void applyJetPtUESubtraction(std::vector<DLJet> &Jets, Event &Event, DLJetCuts& cuts) {
  for (int iJet = 0; iJet < (int)Jets.size(); iJet++) {
    if(cuts.UEEstimationMethod == "JetArea")
      Jets.at(iJet).pt -= Event.Rho * Jets.at(iJet).Area;
    else if(cuts.UEEstimationMethod == "PerpCone")
      Jets.at(iJet).pt -= Jets.at(iJet).PerpConeRho * Jets.at(iJet).Area;
    else
      FATAL("Unknown UEEstimationMethod");
  }
}