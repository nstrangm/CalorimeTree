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
  EventCuts(GlobalOptions optns);
  ~EventCuts(){};
  bool PassedCuts(Event Event);
};

EventCuts::EventCuts(GlobalOptions optns)
{
  // if(EventCutString.Contains("minbias")){

  // }
}

bool EventCuts::PassedCuts(Event Event)
{
  bool passed = true;
  if (Event.ZVtx < zVtxMin || Event.ZVtx > zVtxMax)
  {
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
  float EMax = 0.;

public:
  IsoGammaCuts(GlobalOptions optns);
  ~IsoGammaCuts(){};
  bool useRhoInsteadOfPerpCone = true;
  int NonLinMode = 1;
  bool PassedCuts(IsoGamma IsoGamma);
};

IsoGammaCuts::IsoGammaCuts(GlobalOptions optns)
{
  const char *YAMLFilePath = Form("%s/Cuts.yaml", optns.parentDir.Data());
  YAML::Node ycut = YAML::LoadFile(YAMLFilePath);

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file %s", optns.cutString.Data(), YAMLFilePath))

  YAML::Node cuts = ycut[(std::string)optns.cutString];

  EMin = cuts["cluster_min_E"].as<float>();
  EMax = cuts["cluster_max_E"].as<float>();
  useRhoInsteadOfPerpCone = cuts["useRhoInsteadOfPerpCone"].as<bool>();

  INFO(Form("EMin = %f, EMax = %f, useRhoInsteadOfPerpCone = %s", EMin, EMax, useRhoInsteadOfPerpCone ? "true" : "false"))
}

bool IsoGammaCuts::PassedCuts(IsoGamma IsoGamma)
{
  bool passed = true;
  if (IsoGamma.E < EMin)
  {
    passed = false;
  }

  return passed;
}

void doIsoGammaCuts(std::vector<IsoGamma> IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  for (unsigned long i = 0; i < IsoGammas.size(); i++)
  {
    if (!IsoGammaCuts.PassedCuts(IsoGammas.at(i)))
    {
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
  JetCuts(GlobalOptions optns);
  ~JetCuts(){};
  bool PassedCuts(Jet Jet);
  float R = 0.4;
};

JetCuts::JetCuts(GlobalOptions optns)
{
  // if (JetCutString.Contains("Pt700"))
  // {
  //   PtMin = 0.7;
  // }
}

bool JetCuts::PassedCuts(Jet Jet)
{
  bool passed = true;
  if (Jet.Pt() < PtMin)
  {
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
  Pi0Cuts(GlobalOptions optns);
  ~Pi0Cuts(){};
  bool PassedCuts(Pi0 Pi0);
};

Pi0Cuts::Pi0Cuts(GlobalOptions optns)
{
  // if (Pi0CutString.Contains("Pt700"))
  // {
  //   PtMin = 0.7;
  // }
}

bool Pi0Cuts::PassedCuts(Pi0 Pi0)
{
  bool passed = true;
  if (Pi0.Pt() < PtMin)
  {
    passed = false;
  }

  return passed;
}

#endif