#ifndef _UTILITIES_
#define _UTILITIES_

#include "Logging.h"
#include "PlottingClass.h"

class GlobalOptions
{
 public:
  int TreeFormat = 0;
  bool doQA = false;
  bool doIsoGamma = false;
  bool doJets = false;
  bool doPi0 = false;
  GlobalOptions(bool userWantsQA, bool userWantsIsoGamma, bool userWantsJets, bool userWantsPi0s);
  ~GlobalOptions();
};

GlobalOptions::GlobalOptions(bool userWantsQA, bool userWantsIsoGamma, bool userWantsJets, bool userWantsPi0s)
{
  doQA = userWantsQA;
  doIsoGamma = userWantsIsoGamma;
  doJets = userWantsJets;
  doPi0 = userWantsPi0s;
}

GlobalOptions::~GlobalOptions()
{
}

#endif