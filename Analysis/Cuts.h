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
  bool PassedCuts(Event &Event);
};

EventCuts::EventCuts(GlobalOptions optns)
{
  // if(EventCutString.Contains("minbias")){

  // }
}

bool EventCuts::PassedCuts(Event &Event)
{
  bool passed = true;
  if (Event.ZVtx < zVtxMin || Event.ZVtx > zVtxMax)
  {
    passed = false;
    Event.ZVtxOK = false;
  }
  if (!Event.IsTriggered)
  {
    passed = false;
  }
  if (Event.Quality != 0)
  {
    passed = false;
    Event.QualityOK = false;
  }

  Event.Selected = true;
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
  float EMcalEtaMax = 0.;
  float EMcalEtaMin = 0.;
  float EMcalPhiMax = 0.;
  float EMcalPhiMin = 0.;
  float DcalEtaMax = 0.;
  float DcalEtaMin = 0.;
  float DcalPhiMax = 0.;
  float DcalPhiMin = 0.;
  float DcalHoleEtaMax = 0.;
  float DcalHoleEtaMin = 0.;
  float DcalHolePhiMax = 0.;
  float DcalHolePhiMin = 0.;
  float EMin = 0.;
  float EMax = 0.;
  unsigned short NcellsMin = 0;
  unsigned short NLMMax = 0;
  float DistanceToBadChannelMin = 0;
  float M02Min = 0.;
  float M02Max = 0.;
  float MatchDetaMin = 0;
  float MatchDetaMax = 0;
  float MatchDphiMin = 0;
  float MatchDphiMax = 0;
  float MatchVetoMax = 0;
  float FplusMax = 0;
  float pTisoCorrectedMax = 0;

public:
  IsoGammaCuts(GlobalOptions optns);
  ~IsoGammaCuts(){};
  bool useRhoInsteadOfPerpCone = true;
  int NonLinMode = 1;
  bool applyNonLin = false; // Apply NonLinearity in CalorimeTree only if not already done in AnalysisTask
  bool PassedCuts(IsoGamma IsoGamma);
};

IsoGammaCuts::IsoGammaCuts(GlobalOptions optns)
{
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];

  applyNonLin = chosencut["cluster_applyNonLin"].IsDefined() ? chosencut["cluster_applyNonLin"].as<bool>() : standardcut["cluster_applyNonLin"].as<bool>();
  DEBUG
  // Acceptance cut
  // EMcal
  EMcalEtaMax = chosencut["cluster_max_EMcal_eta"].IsDefined() ? chosencut["cluster_max_EMcal_eta"].as<float>() : standardcut["cluster_max_EMcal_eta"].as<float>();
  EMcalEtaMin = chosencut["cluster_min_EMcal_eta"].IsDefined() ? chosencut["cluster_min_EMcal_eta"].as<float>() : standardcut["cluster_min_EMcal_eta"].as<float>();
  EMcalPhiMax = chosencut["cluster_max_EMcal_phi"].IsDefined() ? chosencut["cluster_max_EMcal_phi"].as<float>() : standardcut["cluster_max_EMcal_phi"].as<float>();
  EMcalPhiMin = chosencut["cluster_min_EMcal_phi"].IsDefined() ? chosencut["cluster_min_EMcal_phi"].as<float>() : standardcut["cluster_min_EMcal_phi"].as<float>();
  // DCal
  DcalEtaMax = chosencut["cluster_max_Dcal_eta"].IsDefined() ? chosencut["cluster_max_Dcal_eta"].as<float>() : standardcut["cluster_max_Dcal_eta"].as<float>();
  DcalEtaMin = chosencut["cluster_min_Dcal_eta"].IsDefined() ? chosencut["cluster_min_Dcal_eta"].as<float>() : standardcut["cluster_min_Dcal_eta"].as<float>();
  DcalPhiMax = chosencut["cluster_max_Dcal_phi"].IsDefined() ? chosencut["cluster_max_Dcal_phi"].as<float>() : standardcut["cluster_max_Dcal_phi"].as<float>();
  DcalPhiMin = chosencut["cluster_min_Dcal_phi"].IsDefined() ? chosencut["cluster_min_Dcal_phi"].as<float>() : standardcut["cluster_min_Dcal_phi"].as<float>();
  // DCal hole
  DcalHoleEtaMax = chosencut["cluster_max_DcalHole_eta"].IsDefined() ? chosencut["cluster_max_DcalHole_eta"].as<float>() : standardcut["cluster_max_DcalHole_eta"].as<float>();
  DcalHoleEtaMin = chosencut["cluster_min_DcalHole_eta"].IsDefined() ? chosencut["cluster_min_DcalHole_eta"].as<float>() : standardcut["cluster_min_DcalHole_eta"].as<float>();
  DcalHolePhiMax = chosencut["cluster_max_DcalHole_phi"].IsDefined() ? chosencut["cluster_max_DcalHole_phi"].as<float>() : standardcut["cluster_max_DcalHole_phi"].as<float>();
  DcalHolePhiMin = chosencut["cluster_min_DcalHole_phi"].IsDefined() ? chosencut["cluster_min_DcalHole_phi"].as<float>() : standardcut["cluster_min_DcalHole_phi"].as<float>();
  // Load E cut
  EMin = chosencut["cluster_min_E"].IsDefined() ? chosencut["cluster_min_E"].as<float>() : standardcut["cluster_min_E"].as<float>();
  EMax = chosencut["cluster_max_E"].IsDefined() ? chosencut["cluster_max_E"].as<float>() : standardcut["cluster_max_E"].as<float>();

  // Load min Ncells
  NcellsMin = chosencut["cluster_min_Nc"].IsDefined() ? chosencut["cluster_min_Nc"].as<unsigned short>() : standardcut["cluster_min_Nc"].as<unsigned short>();

  // Load NLM cut
  NLMMax = chosencut["cluster_max_NLM"].IsDefined() ? chosencut["cluster_max_NLM"].as<unsigned short>() : standardcut["cluster_max_NLM"].as<unsigned short>();

  // Load min dist to bad channel cut
  DistanceToBadChannelMin = chosencut["cluster_min_distbadch"].IsDefined() ? chosencut["cluster_min_distbadch"].as<float>() : standardcut["cluster_min_distbadch"].as<float>();

  // Load Sigma long cut
  M02Min = chosencut["cluster_min_M02"].IsDefined() ? chosencut["cluster_min_M02"].as<float>() : standardcut["cluster_min_M02"].as<float>();
  M02Max = chosencut["cluster_max_M02"].IsDefined() ? chosencut["cluster_max_M02"].as<float>() : standardcut["cluster_max_M02"].as<float>();

  // Load track matching cut
  MatchDetaMin = chosencut["cluster_min_MatchDeta"].IsDefined() ? chosencut["cluster_min_MatchDeta"].as<float>() : standardcut["cluster_min_MatchDeta"].as<float>();
  MatchDetaMax = chosencut["cluster_max_MatchDeta"].IsDefined() ? chosencut["cluster_max_MatchDeta"].as<float>() : standardcut["cluster_max_MatchDeta"].as<float>();
  MatchDphiMin = chosencut["cluster_min_MatchDeta"].IsDefined() ? chosencut["cluster_min_MatchDeta"].as<float>() : standardcut["cluster_min_MatchDeta"].as<float>();
  MatchDphiMax = chosencut["cluster_max_MatchDeta"].IsDefined() ? chosencut["cluster_max_MatchDeta"].as<float>() : standardcut["cluster_max_MatchDeta"].as<float>();
  MatchVetoMax = chosencut["cluster_max_MatchVeto"].IsDefined() ? chosencut["cluster_max_MatchVeto"].as<float>() : standardcut["cluster_max_MatchVeto"].as<float>();

  // Load FM+
  FplusMax = chosencut["cluster_max_Fplus"].IsDefined() ? chosencut["cluster_max_Fplus"].as<float>() : standardcut["cluster_max_Fplus"].as<float>();

  // Load Isolation cut
  pTisoCorrectedMax = chosencut["cluster_max_ptiso"].IsDefined() ? chosencut["cluster_max_ptiso"].as<float>() : standardcut["cluster_max_ptiso"].as<float>();

  useRhoInsteadOfPerpCone = chosencut["useRhoInsteadOfPerpCone"].IsDefined() ? chosencut["useRhoInsteadOfPerpCone"].as<bool>() : standardcut["useRhoInsteadOfPerpCone"].as<bool>();

  INFO(Form("EMin = %f, EMax = %f, useRhoInsteadOfPerpCone = %s", EMin, EMax, useRhoInsteadOfPerpCone ? "true" : "false"))
}

bool IsoGammaCuts::PassedCuts(IsoGamma IsoGamma)
{
  bool passed = true;

  // Check EMCAL acceptance
  if (IsoGamma.Eta() < EMcalEtaMax && IsoGamma.Eta() > EMcalEtaMin && IsoGamma.Phi() > EMcalPhiMin && IsoGamma.Phi() < EMcalPhiMax)
  {
    // True if is in EMCAL.
  }else{
    if (IsoGamma.Eta() < DcalEtaMax && IsoGamma.Eta() > DcalEtaMin && IsoGamma.Phi() > DcalPhiMin && IsoGamma.Phi() < DcalPhiMax)
    {// True if is in DCal.
      if (IsoGamma.Eta() > DcalHoleEtaMin && IsoGamma.Eta() < DcalHoleEtaMax && IsoGamma.Phi() < DcalHolePhiMin && IsoGamma.Phi() > DcalHolePhiMax)
      {
        passed = false; // But false if is in Dcal and in Dcal hole.
      }
    }
    else
    {
      passed = false;
    }
  }

  // check cluster energy
  if (IsoGamma.E < EMin || IsoGamma.E > EMax)
  {
    passed = false;
  }
  // check cells pr. cluster
  if (IsoGamma.NCells < NcellsMin)
  {
    passed = false;
  }
  // Check number of local maxima
  if (IsoGamma.NLM > NLMMax)
  {
    passed = false;
  }
  // Check distance to bad channel
  if (IsoGamma.DistanceToBadChannel < DistanceToBadChannelMin)
  {
    passed = false;
  }
  // check cluster long-axis
  if (IsoGamma.M02 < M02Min || IsoGamma.M02 > M02Max)
  {
    passed = false;
  }
  // Check track matching
  if (IsoGamma.MatchedTrack.dEta > MatchDetaMax || IsoGamma.MatchedTrack.dEta < MatchDetaMin || IsoGamma.MatchedTrack.dPhi < MatchDphiMin || IsoGamma.MatchedTrack.dPhi > MatchDphiMax)
  {if ((IsoGamma.E / IsoGamma.MatchedTrack.P) > MatchVetoMax)//Check veto
    {
      passed = false;
    }
  }
  // Check Fplus
  if (IsoGamma.EFrac > FplusMax)
  {
    passed = false;
  }
  // Check ptiso (fully corrected quantity)
  if (IsoGamma.IsoChargedCorrected > pTisoCorrectedMax)
  {
    passed = false;
  }
  // Check number of local maxima
  if (IsoGamma.NLM > NLMMax)
  {
    passed = false;
  }
  return passed;
}

void doIsoGammaCuts(std::vector<IsoGamma> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<IsoGamma>::iterator iter;
  for (iter = IsoGammas.begin(); iter != IsoGammas.end();)
  {
    if (!IsoGammaCuts.PassedCuts(*iter))
    {
      iter = IsoGammas.erase(iter);
      // cout<<"Out"<<"\n";
    }
    else
    {
      ++iter;
      // cout<<"Passed"<<"\n";
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