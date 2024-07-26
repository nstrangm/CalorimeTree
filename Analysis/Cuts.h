#ifndef CUTS
#define CUTS

#include "PhysicsObjects.h"
#include "Logging.h"

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
// ------------------- GammaGen Cuts -------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
class GammaGenCuts
{
private:
  float EMCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
  float DCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
  float DCalHoleEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};

  float pTisoCorrectedMax = 0;
  float RapidityMax = 0;
  TDirectory *hQADir = nullptr;

public:
  GammaGenCuts(GlobalOptions optns);
  ~GammaGenCuts(){};

  bool PassedGammaGenCuts(GammaGen GammaGen);
  bool isSignal(GammaGen GammaGen);
  bool PassedGammaGenIsolationCuts(GammaGen GammaGen);
};

GammaGenCuts::GammaGenCuts(GlobalOptions optns)
{
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];
  // Acceptance cut
  // EMcal
  EMCalEtaPhiMinMax[0][0] = chosencut["cluster_min_EMcal_eta"].IsDefined() ? chosencut["cluster_min_EMcal_eta"].as<float>() : standardcut["cluster_min_EMcal_eta"].as<float>();
  EMCalEtaPhiMinMax[0][1] = chosencut["cluster_max_EMcal_eta"].IsDefined() ? chosencut["cluster_max_EMcal_eta"].as<float>() : standardcut["cluster_max_EMcal_eta"].as<float>();
  EMCalEtaPhiMinMax[1][0] = chosencut["cluster_min_EMcal_phi"].IsDefined() ? chosencut["cluster_min_EMcal_phi"].as<float>() : standardcut["cluster_min_EMcal_phi"].as<float>();
  EMCalEtaPhiMinMax[1][1] = chosencut["cluster_max_EMcal_phi"].IsDefined() ? chosencut["cluster_max_EMcal_phi"].as<float>() : standardcut["cluster_max_EMcal_phi"].as<float>();
  // DCal
  DCalEtaPhiMinMax[0][0] = chosencut["cluster_min_Dcal_eta"].IsDefined() ? chosencut["cluster_min_Dcal_eta"].as<float>() : standardcut["cluster_min_Dcal_eta"].as<float>();
  DCalEtaPhiMinMax[0][1] = chosencut["cluster_max_Dcal_eta"].IsDefined() ? chosencut["cluster_max_Dcal_eta"].as<float>() : standardcut["cluster_max_Dcal_eta"].as<float>();
  DCalEtaPhiMinMax[1][0] = chosencut["cluster_min_Dcal_phi"].IsDefined() ? chosencut["cluster_min_Dcal_phi"].as<float>() : standardcut["cluster_min_Dcal_phi"].as<float>();
  DCalEtaPhiMinMax[1][1] = chosencut["cluster_max_Dcal_phi"].IsDefined() ? chosencut["cluster_max_Dcal_phi"].as<float>() : standardcut["cluster_max_Dcal_phi"].as<float>();
  // DCal hole
  DCalHoleEtaPhiMinMax[0][0] = chosencut["cluster_min_DcalHole_eta"].IsDefined() ? chosencut["cluster_min_DcalHole_eta"].as<float>() : standardcut["cluster_min_DcalHole_eta"].as<float>();
  DCalHoleEtaPhiMinMax[0][1] = chosencut["cluster_max_DcalHole_eta"].IsDefined() ? chosencut["cluster_max_DcalHole_eta"].as<float>() : standardcut["cluster_max_DcalHole_eta"].as<float>();
  DCalHoleEtaPhiMinMax[1][0] = chosencut["cluster_min_DcalHole_phi"].IsDefined() ? chosencut["cluster_min_DcalHole_phi"].as<float>() : standardcut["cluster_min_DcalHole_phi"].as<float>();
  DCalHoleEtaPhiMinMax[1][1] = chosencut["cluster_max_DcalHole_phi"].IsDefined() ? chosencut["cluster_max_DcalHole_phi"].as<float>() : standardcut["cluster_max_DcalHole_phi"].as<float>();
  // Load Isolation cut
  pTisoCorrectedMax = chosencut["cluster_max_ptiso"].IsDefined() ? chosencut["cluster_max_ptiso"].as<float>() : standardcut["cluster_max_ptiso"].as<float>();
  RapidityMax = chosencut["cluster_max_rapidity"].IsDefined() ? chosencut["cluster_max_rapidity"].as<float>() : standardcut["cluster_max_rapidity"].as<float>();
  //useRhoInsteadOfPerpCone = chosencut["useRhoInsteadOfPerpCone"].IsDefined() ? chosencut["useRhoInsteadOfPerpCone"].as<bool>() : standardcut["useRhoInsteadOfPerpCone"].as<bool>();
}

bool GammaGenCuts::PassedGammaGenCuts(GammaGen GammaGen)
{
  bool passed = true;
  //Raise fatal if eta for MC is outside allowed range:
  if(TMath::Abs(GammaGen.Eta())>RapidityMax){
    FATAL(Form("Generated photon found with |eta|>%f  |%f|>%f",RapidityMax,GammaGen.Eta(),RapidityMax))
  }
  // Check GammaGen acceptance
  if (!GammaGen.isInEMCalAcceptance(EMCalEtaPhiMinMax) && !GammaGen.isInDCalAcceptance(DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax))
  {
    passed = false;
  }
  return passed;
}

bool GammaGenCuts::isSignal(GammaGen GammaGen)
{

  if (CheckTagBit(GammaGen.MCTag, kMCPhoton))
  {
    if ((CheckTagBit(GammaGen.MCTag, kMCPrompt) || CheckTagBit(GammaGen.MCTag, kMCFragmentation)) && GammaGen.IsoChargedCorrected < pTisoCorrectedMax)
    {
      return true;
    }
  }
  return false;
}

void doGammaGenCuts(std::vector<GammaGen> &GammaGens, GammaGenCuts GammaGenCuts)
{
  std::vector<GammaGen>::iterator iter;
  for (iter = GammaGens.begin(); iter != GammaGens.end();)
  {
    if (!GammaGenCuts.PassedGammaGenCuts(*iter))
    {
      iter = GammaGens.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

bool GammaGenCuts::PassedGammaGenIsolationCuts(GammaGen GammaGen)
{
  bool passed = true;
  // check cluster long-axis
  if (GammaGen.IsoChargedCorrected > pTisoCorrectedMax)
  {
    passed = false;
  }
  return passed;
}

void doGammaGenIsolationCuts(std::vector<GammaGen> &GammaGens, GammaGenCuts GammaGenCuts)
{
  std::vector<GammaGen>::iterator iter;
  for(iter = GammaGens.begin();iter != GammaGens.end();)
  {
    if (!GammaGenCuts.PassedGammaGenIsolationCuts(*iter))
    {
      iter = GammaGens.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
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
  float EMCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
  float DCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
  float DCalHoleEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
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
  TDirectory *hQADir = nullptr;

public:
  IsoGammaCuts(GlobalOptions optns, TDirectory *hQADirIsoGammas);
  ~IsoGammaCuts(){};

  bool useRhoInsteadOfPerpCone = true;
  int NonLinMode = 1;
  bool applyNonLin = false; // Apply NonLinearity in CalorimeTree only if not already done in AnalysisTask
  bool PassedIsoGammaCuts(IsoGamma IsoGamma);
  bool PassedClusterCuts(IsoGamma IsoGamma);
  bool isSignal(IsoGamma IsoGamma);
  bool isSignalClusterLevelIso(IsoGamma IsoGamma);//For troubleshooting, use "isSignal" for actual signal def.
};

IsoGammaCuts::IsoGammaCuts(GlobalOptions optns, TDirectory *hQADirIsoGammas)
{
  hQADir = hQADirIsoGammas;
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];

  applyNonLin = chosencut["cluster_applyNonLin"].IsDefined() ? chosencut["cluster_applyNonLin"].as<bool>() : standardcut["cluster_applyNonLin"].as<bool>();
  DEBUG
  // Acceptance cut
  // EMcal
  EMCalEtaPhiMinMax[0][0] = chosencut["cluster_min_EMcal_eta"].IsDefined() ? chosencut["cluster_min_EMcal_eta"].as<float>() : standardcut["cluster_min_EMcal_eta"].as<float>();
  EMCalEtaPhiMinMax[0][1] = chosencut["cluster_max_EMcal_eta"].IsDefined() ? chosencut["cluster_max_EMcal_eta"].as<float>() : standardcut["cluster_max_EMcal_eta"].as<float>();
  EMCalEtaPhiMinMax[1][0] = chosencut["cluster_min_EMcal_phi"].IsDefined() ? chosencut["cluster_min_EMcal_phi"].as<float>() : standardcut["cluster_min_EMcal_phi"].as<float>();
  EMCalEtaPhiMinMax[1][1] = chosencut["cluster_max_EMcal_phi"].IsDefined() ? chosencut["cluster_max_EMcal_phi"].as<float>() : standardcut["cluster_max_EMcal_phi"].as<float>();
  // DCal
  DCalEtaPhiMinMax[0][0] = chosencut["cluster_min_Dcal_eta"].IsDefined() ? chosencut["cluster_min_Dcal_eta"].as<float>() : standardcut["cluster_min_Dcal_eta"].as<float>();
  DCalEtaPhiMinMax[0][1] = chosencut["cluster_max_Dcal_eta"].IsDefined() ? chosencut["cluster_max_Dcal_eta"].as<float>() : standardcut["cluster_max_Dcal_eta"].as<float>();
  DCalEtaPhiMinMax[1][0] = chosencut["cluster_min_Dcal_phi"].IsDefined() ? chosencut["cluster_min_Dcal_phi"].as<float>() : standardcut["cluster_min_Dcal_phi"].as<float>();
  DCalEtaPhiMinMax[1][1] = chosencut["cluster_max_Dcal_phi"].IsDefined() ? chosencut["cluster_max_Dcal_phi"].as<float>() : standardcut["cluster_max_Dcal_phi"].as<float>();
  // DCal hole
  DCalHoleEtaPhiMinMax[0][0] = chosencut["cluster_min_DcalHole_eta"].IsDefined() ? chosencut["cluster_min_DcalHole_eta"].as<float>() : standardcut["cluster_min_DcalHole_eta"].as<float>();
  DCalHoleEtaPhiMinMax[0][1] = chosencut["cluster_max_DcalHole_eta"].IsDefined() ? chosencut["cluster_max_DcalHole_eta"].as<float>() : standardcut["cluster_max_DcalHole_eta"].as<float>();
  DCalHoleEtaPhiMinMax[1][0] = chosencut["cluster_min_DcalHole_phi"].IsDefined() ? chosencut["cluster_min_DcalHole_phi"].as<float>() : standardcut["cluster_min_DcalHole_phi"].as<float>();
  DCalHoleEtaPhiMinMax[1][1] = chosencut["cluster_max_DcalHole_phi"].IsDefined() ? chosencut["cluster_max_DcalHole_phi"].as<float>() : standardcut["cluster_max_DcalHole_phi"].as<float>();
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

bool IsoGammaCuts::PassedIsoGammaCuts(IsoGamma IsoGamma)
{
  bool passed = true;
  // check cluster long-axis
  if (IsoGamma.M02 < M02Min || IsoGamma.M02 > M02Max)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(8., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(8., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check ptiso (fully corrected quantity)
  if (IsoGamma.IsoChargedCorrected > pTisoCorrectedMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);
  return passed;
}

bool IsoGammaCuts::PassedClusterCuts(IsoGamma IsoGamma)
{
  bool passed = true;
  //Check if event was loaded.
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(0., IsoGamma.Pt(), IsoGamma.EventWeight);
  if (hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(0., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check cluster acceptance
  if (!IsoGamma.isInEMCalAcceptance(EMCalEtaPhiMinMax) && !IsoGamma.isInDCalAcceptance(DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax))
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(1., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(1., IsoGamma.Pt(), IsoGamma.EventWeight);
  // check cluster energy
  if (IsoGamma.E < EMin || IsoGamma.E > EMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(2., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(2., IsoGamma.Pt(), IsoGamma.EventWeight);
  // check cells pr. cluster
  if (IsoGamma.NCells < NcellsMin)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(3., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(3., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check number of local maxima
  if (IsoGamma.NLM > NLMMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(4., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(4., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check distance to bad channel
  if (IsoGamma.DistanceToBadChannel < DistanceToBadChannelMin)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(5., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(5., IsoGamma.Pt(), IsoGamma.EventWeight);
   //Check track matching
  if (IsoGamma.MatchedTrack.P>0 && IsoGamma.MatchedTrack.dEta < MatchDetaMax && IsoGamma.MatchedTrack.dEta > MatchDetaMin && IsoGamma.MatchedTrack.dPhi > MatchDphiMin && IsoGamma.MatchedTrack.dPhi < MatchDphiMax && (IsoGamma.E / IsoGamma.MatchedTrack.P) < MatchVetoMax)
  {
    if (hQADir != nullptr){
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(6., IsoGamma.Pt(), IsoGamma.EventWeight);
      ((TH2F *)hQADir->FindObject("hIsoGammadEtadphiCut"))->Fill(IsoGamma.MatchedTrack.dPhi, IsoGamma.MatchedTrack.dEta, IsoGamma.EventWeight);
    }
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(6., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check Fplus
  if (IsoGamma.EFrac > FplusMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(7., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(7., IsoGamma.Pt(), IsoGamma.EventWeight);
  return passed;
}

bool IsoGammaCuts::isSignal(IsoGamma IsoGamma)
{
  //Calculate generator level isolation:
  float pTisoGen= IsoGamma.TrueClusterIsoCharged / CalculateIsoCorrectionFactor(IsoGamma.Eta(), 0.8, 0.4) - IsoGamma.TrueClusterIsoBckPerp * TMath::Pi() * 0.4 * 0.4;
  if (CheckTagBit(IsoGamma.MCTag, kMCPhoton))
  {
    if (CheckTagBit(IsoGamma.MCTag, kMCPrompt) || CheckTagBit(IsoGamma.MCTag, kMCFragmentation) && pTisoGen < pTisoCorrectedMax)
    {
      return true;
    }
  }
  return false;
}

bool IsoGammaCuts::isSignalClusterLevelIso(IsoGamma IsoGamma)//For troubleshooting, use "isSignal" for actual signal def.
{
  //Calculate generator level isolation:
  float pTisoGen= IsoGamma.IsoCharged / CalculateIsoCorrectionFactor(IsoGamma.Eta(), 0.8, 0.4) - IsoGamma.IsoBckPerp * TMath::Pi() * 0.4 * 0.4;
  if (CheckTagBit(IsoGamma.MCTag, kMCPhoton))
  {
    if (CheckTagBit(IsoGamma.MCTag, kMCPrompt) || CheckTagBit(IsoGamma.MCTag, kMCFragmentation) && pTisoGen < pTisoCorrectedMax)
    {
      return true;
    }
  }
  return false;
}

void doIsoGammaClusterCuts(std::vector<IsoGamma> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<IsoGamma>::iterator iter;
  for (iter = IsoGammas.begin(); iter != IsoGammas.end();)
  {
    if (!IsoGammaCuts.PassedClusterCuts(*iter))
    {
      iter = IsoGammas.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

void doIsoGammaCuts(std::vector<IsoGamma> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<IsoGamma>::iterator iter;
  for (iter = IsoGammas.begin(); iter != IsoGammas.end();)
  {
    if (!IsoGammaCuts.PassedIsoGammaCuts(*iter))
    {
      iter = IsoGammas.erase(iter);
    }
    else
    {
      ++iter; 
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