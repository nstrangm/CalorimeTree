#ifndef CUTS
#define CUTS

#include "PhysicsObjects.h"

// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// -------------------- Event Cuts ---------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
class EventCuts
{
 private:
  float zVtxMin = -10;
  float zVtxMax = 10;

 public:
  EventCuts(TString EventCutString);
  ~EventCuts(){};
  bool PassedCuts(Event Event);
};

EventCuts::EventCuts(TString EventCutString)
{
  // if(EventCutString.Contains("minbias")){

  // }
}

bool EventCuts::PassedCuts(Event Event)
{
  bool passed = true;
  if (Event.ZVtx < zVtxMin || Event.ZVtx > zVtxMax) {
    passed = false;
  }

  return passed;
}

// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// ------------------- IsoGamma Cuts -------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
class IsoGammaCuts
{
 private:
  float EMin = 0.;

 public:
  IsoGammaCuts(TString IsoGammaCutString);
  ~IsoGammaCuts(){};
  bool useRhoInsteadOfPerpCone = true;
  int NonLinMode = 1;
  bool PassedCuts(IsoGamma IsoGamma);
};

IsoGammaCuts::IsoGammaCuts(TString IsoGammaCutString)
{
  useRhoInsteadOfPerpCone = true;
  if (IsoGammaCutString.Contains("E700")) {
    EMin = 0.7;
  }
}

bool IsoGammaCuts::PassedCuts(IsoGamma IsoGamma)
{
  bool passed = true;
  if (IsoGamma.E < EMin) {
    passed = false;
  }

  return passed;
}

void doIsoGammaCuts(std::vector<IsoGamma> IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  for (unsigned long i = 0; i < IsoGammas.size(); i++) {
    if (!IsoGammaCuts.PassedCuts(IsoGammas.at(i))){
      IsoGammas.erase(IsoGammas.begin() + i);
    }
  }
}

// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// ----------------------- Jet Cuts --------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
class JetCuts
{
 private:
  float PtMin = 0.;

 public:
  JetCuts(TString JetCutString);
  ~JetCuts(){};
  bool PassedCuts(Jet Jet);
  float R = 0.4;
};

JetCuts::JetCuts(TString JetCutString)
{
  if (JetCutString.Contains("Pt700")) {
    PtMin = 0.7;
  }
}

bool JetCuts::PassedCuts(Jet Jet)
{
  bool passed = true;
  if (Jet.Pt() < PtMin) {
    passed = false;
  }

  return passed;
}

// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// ----------------------- Pi0 Cuts --------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
class Pi0Cuts
{
 private:
  float PtMin = 0.;

 public:
  Pi0Cuts(TString Pi0CutString);
  ~Pi0Cuts(){};
  bool PassedCuts(Pi0 Pi0);
};

Pi0Cuts::Pi0Cuts(TString Pi0CutString)
{
  if (Pi0CutString.Contains("Pt700")) {
    PtMin = 0.7;
  }
}

bool Pi0Cuts::PassedCuts(Pi0 Pi0)
{
  bool passed = true;
  if (Pi0.Pt() < PtMin) {
    passed = false;
  }

  return passed;
}

#endif