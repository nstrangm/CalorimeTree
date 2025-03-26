#ifndef CUTS
#define CUTS

#include "PhysicsObjects.h"
#include "Logging.h"
#include "TRandom3.h"

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
  int OccupancyMin = 0;
  int OccupancyMax = 99999;

public:
  EventCuts(GlobalOptions optns);
  ~EventCuts() {};
  bool PassedCuts(Event &Event);
};

EventCuts::EventCuts(GlobalOptions optns)
{
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))
  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];
  OccupancyMin = chosencut["event_occupancy_min"].IsDefined() ? chosencut["event_occupancy_min"].as<int>() : standardcut["event_occupancy_min"].as<int>();
  OccupancyMax = chosencut["event_occupancy_max"].IsDefined() ? chosencut["event_occupancy_max"].as<int>() : standardcut["event_occupancy_max"].as<int>();
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
  if (Event.Occupancy < OccupancyMin || Event.Occupancy > OccupancyMax)
  {
    passed = false;
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
  ~GammaGenCuts() {};

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
  TString ClusterAcceptance = (TString)(chosencut["gamma_Acceptance"].IsDefined() ? chosencut["gamma_Acceptance"].as<std::string>() : standardcut["gamma_Acceptance"].as<std::string>()).c_str();
  SetAcceptance(ClusterAcceptance, EMCalEtaPhiMinMax, DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax);
  // Load Isolation cut
  pTisoCorrectedMax = chosencut["gamma_max_ptiso"].IsDefined() ? chosencut["gamma_max_ptiso"].as<float>() : standardcut["gamma_max_ptiso"].as<float>();
  RapidityMax = chosencut["gamma_max_rapidity"].IsDefined() ? chosencut["gamma_max_rapidity"].as<float>() : standardcut["gamma_max_rapidity"].as<float>();
  
  LOG(Form("Gen Gamma %s Acceptance: %.2f < EMCal_Eta < %.2f | %.2f < EMCal_Phi < %.2f | %.2f < DCal_Eta < %.2f | %.2f < DCal_Phi < %.2f | %.2f < DCalHole_Eta < %.2f | %.2f < DCalHole_Phi < %.2f", ClusterAcceptance.Data(), EMCalEtaPhiMinMax[0][0], EMCalEtaPhiMinMax[0][1], EMCalEtaPhiMinMax[1][0], EMCalEtaPhiMinMax[1][1], DCalEtaPhiMinMax[0][0], DCalEtaPhiMinMax[0][1], DCalEtaPhiMinMax[1][0], DCalEtaPhiMinMax[1][1], DCalHoleEtaPhiMinMax[0][0], DCalHoleEtaPhiMinMax[0][1], DCalHoleEtaPhiMinMax[1][0], DCalHoleEtaPhiMinMax[1][1]))
}

bool GammaGenCuts::PassedGammaGenCuts(GammaGen GammaGen)
{
  bool passed = true;
  // Raise fatal if eta for MC is outside allowed range:
  if (TMath::Abs(GammaGen.Eta()) > RapidityMax)
  {
    FATAL(Form("Generated photon found with |eta|>%f  |%f|>%f", RapidityMax, GammaGen.Eta(), RapidityMax))
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
  for (iter = GammaGens.begin(); iter != GammaGens.end();)
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
  int ClusterDefinition = 10;
  float TimeMin = -1000;
  float TimeMax = +1000;
  unsigned short NcellsMin = 0;
  unsigned short NLMMax = 0;
  float DistanceToBadChannelMin = 0;
  float M02Min = 0.;
  float M02Max = 0.;
  bool doTrackMatching = false;
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
  ~IsoGammaCuts() {};

  bool useRhoInsteadOfPerpCone = true;
  int NonLinMode = 1;
  bool applyNonLin = false; // Apply NonLinearity in CalorimeTree only if not already done in AnalysisTask
  bool PassedShowerShapeCuts(Cluster IsoGamma);
  bool PassedIsoGammaCuts(Cluster IsoGamma);
  bool PassedClusterCuts(Cluster IsoGamma);
  bool isSignal(Cluster IsoGamma);
  bool isSignalClusterLevelIso(Cluster IsoGamma); // For troubleshooting, use "isSignal" for actual signal def.
};

IsoGammaCuts::IsoGammaCuts(GlobalOptions optns, TDirectory *hQADirIsoGammas)
{
  hQADir = hQADirIsoGammas;
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  bool bmerged = ((TString)hQADirIsoGammas->GetTitle()).Contains("mPi0");
  bool bgammaForPi0 = ((TString)hQADirIsoGammas->GetTitle()).Contains("gammaForPi0");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];

  applyNonLin = chosencut["gamma_applyNonLin"].IsDefined() ? chosencut["gamma_applyNonLin"].as<bool>() : standardcut["gamma_applyNonLin"].as<bool>();

  // Acceptance cuts
  TString ClusterAcceptance = (TString)(chosencut["gamma_Acceptance"].IsDefined() ? chosencut["gamma_Acceptance"].as<std::string>() : standardcut["gamma_Acceptance"].as<std::string>()).c_str();
  SetAcceptance(ClusterAcceptance, EMCalEtaPhiMinMax, DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax);

  // Load E cut
  const char* minE = bgammaForPi0 ? "gPi0_min_E" : "gamma_min_E";
  const char* maxE = bgammaForPi0 ? "gPi0_max_E" : "gamma_max_E";
  EMin = chosencut[minE].IsDefined() ? chosencut[minE].as<float>() : standardcut[minE].as<float>();
  EMax = chosencut[maxE].IsDefined() ? chosencut[maxE].as<float>() : standardcut[maxE].as<float>();

  // Load Time cut
  TimeMin = chosencut["gamma_min_time"].IsDefined() ? chosencut["gamma_min_time"].as<float>() : standardcut["gamma_min_time"].as<float>();
  TimeMax = chosencut["gamma_max_time"].IsDefined() ? chosencut["gamma_max_time"].as<float>() : standardcut["gamma_max_time"].as<float>();

  // Load min Ncells
  NcellsMin = chosencut["gamma_min_Nc"].IsDefined() ? chosencut["gamma_min_Nc"].as<unsigned short>() : standardcut["gamma_min_Nc"].as<unsigned short>();

  // Load NLM cut
  NLMMax = chosencut["gamma_max_NLM"].IsDefined() ? chosencut["gamma_max_NLM"].as<unsigned short>() : standardcut["gamma_max_NLM"].as<unsigned short>();

  ClusterDefinition = chosencut["gamma_ClusterDefinition"].IsDefined() ? chosencut["gamma_ClusterDefinition"].as<int>() : standardcut["gamma_ClusterDefinition"].as<int>();
  // Load min dist to bad channel cut
  DistanceToBadChannelMin = chosencut["gamma_min_distbadch"].IsDefined() ? chosencut["gamma_min_distbadch"].as<float>() : standardcut["gamma_min_distbadch"].as<float>();

  // Load Sigma long cut
  const char* minM02 = bmerged ? "mPi0_min_M02" : "gamma_min_M02";
  const char* maxM02 = bmerged ? "mPi0_max_M02" : "gamma_max_M02";
  if(bgammaForPi0)
  {
    minM02 = "gPi0_min_M02";
    maxM02 = "gPi0_max_M02";
  }
  M02Min = chosencut[minM02].IsDefined() ? chosencut[minM02].as<float>() : standardcut[minM02].as<float>();
  M02Max = chosencut[maxM02].IsDefined() ? chosencut[maxM02].as<float>() : standardcut[maxM02].as<float>();

  // Load track matching cut
  doTrackMatching = chosencut["gamma_doTrackMatching"].IsDefined() ? chosencut["gamma_doTrackMatching"].as<bool>() : standardcut["gamma_doTrackMatching"].as<bool>();
  MatchDetaMin = chosencut["gamma_min_MatchDeta"].IsDefined() ? chosencut["gamma_min_MatchDeta"].as<float>() : standardcut["gamma_min_MatchDeta"].as<float>();
  MatchDetaMax = chosencut["gamma_max_MatchDeta"].IsDefined() ? chosencut["gamma_max_MatchDeta"].as<float>() : standardcut["gamma_max_MatchDeta"].as<float>();
  MatchDphiMin = chosencut["gamma_min_MatchDeta"].IsDefined() ? chosencut["gamma_min_MatchDeta"].as<float>() : standardcut["gamma_min_MatchDeta"].as<float>();
  MatchDphiMax = chosencut["gamma_max_MatchDeta"].IsDefined() ? chosencut["gamma_max_MatchDeta"].as<float>() : standardcut["gamma_max_MatchDeta"].as<float>();
  MatchVetoMax = chosencut["gamma_max_MatchVeto"].IsDefined() ? chosencut["gamma_max_MatchVeto"].as<float>() : standardcut["gamma_max_MatchVeto"].as<float>();

  // Load FM+
  FplusMax = chosencut["gamma_max_Fplus"].IsDefined() ? chosencut["gamma_max_Fplus"].as<float>() : standardcut["gamma_max_Fplus"].as<float>();

  // Load Isolation cut
  const char* maxptiso = bmerged ? "mPi0_max_ptiso" : "gamma_max_ptiso";
  pTisoCorrectedMax = chosencut[maxptiso].IsDefined() ? chosencut[maxptiso].as<float>() : standardcut[maxptiso].as<float>();

  useRhoInsteadOfPerpCone = chosencut["useRhoInsteadOfPerpCone"].IsDefined() ? chosencut["useRhoInsteadOfPerpCone"].as<bool>() : standardcut["useRhoInsteadOfPerpCone"].as<bool>();

  LOG(Form("%s IsoGammaCuts: applyNonLin = %s, %.1f < E < %.1f, %.1f < Time < %.1f, NCells >= %d, NLM <= %d, DistToBadChannel >= %.1f, %.1f < M02 < %.1f, F+ <= %.1f, pTIso < %.1f GeV/c, useRhoInsteadOfPerpCone = %s", optns.cutString.Data(), applyNonLin ? "true" : "false", EMin, EMax,TimeMin,TimeMax, NcellsMin, NLMMax, DistanceToBadChannelMin, M02Min, M02Max, FplusMax, pTisoCorrectedMax, useRhoInsteadOfPerpCone ? "true" : "false"))
  LOG(Form("Gamma %s Acceptance: %.2f < EMCal_Eta < %.2f | %.2f < EMCal_Phi < %.2f | %.2f < DCal_Eta < %.2f | %.2f < DCal_Phi < %.2f | %.2f < DCalHole_Eta < %.2f | %.2f < DCalHole_Phi < %.2f", ClusterAcceptance.Data(), EMCalEtaPhiMinMax[0][0], EMCalEtaPhiMinMax[0][1], EMCalEtaPhiMinMax[1][0], EMCalEtaPhiMinMax[1][1], DCalEtaPhiMinMax[0][0], DCalEtaPhiMinMax[0][1], DCalEtaPhiMinMax[1][0], DCalEtaPhiMinMax[1][1], DCalHoleEtaPhiMinMax[0][0], DCalHoleEtaPhiMinMax[0][1], DCalHoleEtaPhiMinMax[1][0], DCalHoleEtaPhiMinMax[1][1]))
}

// additional function in case one wants apply shower shape cut but not isolation cut
bool IsoGammaCuts::PassedShowerShapeCuts(Cluster IsoGamma)
{
  bool passed = true;
  // check cluster long-axis
  if (IsoGamma.M02 < M02Min || IsoGamma.M02 > M02Max)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);
  return passed;
}

bool IsoGammaCuts::PassedIsoGammaCuts(Cluster IsoGamma)
{
  bool passed = true;
  // check cluster long-axis
  if (!PassedShowerShapeCuts(IsoGamma)){
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(9., IsoGamma.Pt(), IsoGamma.EventWeight);

  // Check ptiso (fully corrected quantity)
  if (IsoGamma.IsoChargedCorrected > pTisoCorrectedMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(10., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(10., IsoGamma.Pt(), IsoGamma.EventWeight);
  return passed;
}

bool IsoGammaCuts::PassedClusterCuts(Cluster IsoGamma)
{
  bool passed = true;
  // Check if event was loaded.
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(0., IsoGamma.Pt(), IsoGamma.EventWeight);
  if (hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(0., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check cluster definition
  if (IsoGamma.Definition != ClusterDefinition)
  {
    // LOG(Form("Cluster definition %d does not match expected definition %d", IsoGamma.Definition, ClusterDefinition))
    passed = false;
  }
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
  // check cluster time
  if (IsoGamma.Time < TimeMin || IsoGamma.Time > TimeMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(3., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(3., IsoGamma.Pt(), IsoGamma.EventWeight);
  
  // check cells per cluster
  if (IsoGamma.NCells < NcellsMin)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(4., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(4., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check number of local maxima
  if (IsoGamma.NLM > NLMMax)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(5., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(5., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check distance to bad channel
  if (IsoGamma.DistanceToBadChannel < DistanceToBadChannelMin)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(6., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(6., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check track matching
  if (doTrackMatching && IsoGamma.MatchedTrack.P > 0 && IsoGamma.MatchedTrack.dEta < MatchDetaMax && IsoGamma.MatchedTrack.dEta > MatchDetaMin && IsoGamma.MatchedTrack.dPhi > MatchDphiMin && IsoGamma.MatchedTrack.dPhi < MatchDphiMax && (IsoGamma.E / IsoGamma.MatchedTrack.P) < MatchVetoMax)
  {
    passed = false;
    if (hQADir != nullptr)
    {
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(7., IsoGamma.Pt(), IsoGamma.EventWeight);
      ((TH2F *)hQADir->FindObject("hIsoGammadEtadphiCut"))->Fill(IsoGamma.MatchedTrack.dPhi, IsoGamma.MatchedTrack.dEta, IsoGamma.EventWeight);
    }
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(7., IsoGamma.Pt(), IsoGamma.EventWeight);
  // Check Fplus
  if (IsoGamma.EFrac > FplusMax || IsoGamma.IsExotic)
  {
    passed = false;
    if (hQADir != nullptr)
      ((TH2F *)hQADir->FindObject("hpTSpectrumLossFromIndividualCuts"))->Fill(8., IsoGamma.Pt(), IsoGamma.EventWeight);
  }
  if (passed && hQADir != nullptr)
    ((TH2F *)hQADir->FindObject("hpTSpectraAfterSubsequentCuts"))->Fill(8., IsoGamma.Pt(), IsoGamma.EventWeight);
  return passed;
}

bool IsoGammaCuts::isSignal(Cluster IsoGamma)
{
  // Calculate generator level isolation:
  float pTisoGen = IsoGamma.TrueClusterIsoCharged / CalculateIsoCorrectionFactor(IsoGamma.Eta(), 0.8, 0.4) - IsoGamma.TrueClusterIsoBckPerp * TMath::Pi() * 0.4 * 0.4;
  if (CheckTagBit(IsoGamma.MCTag, kMCPhoton))
  {
    if (CheckTagBit(IsoGamma.MCTag, kMCPrompt) || CheckTagBit(IsoGamma.MCTag, kMCFragmentation) && pTisoGen < pTisoCorrectedMax)
    {
      return true;
    }
  }
  return false;
}

bool IsoGammaCuts::isSignalClusterLevelIso(Cluster IsoGamma) // For troubleshooting, use "isSignal" for actual signal def.
{
  // Calculate generator level isolation:
  float pTisoGen = IsoGamma.IsoCharged / CalculateIsoCorrectionFactor(IsoGamma.Eta(), 0.8, 0.4) - IsoGamma.IsoBckPerp * TMath::Pi() * 0.4 * 0.4;
  if (CheckTagBit(IsoGamma.MCTag, kMCPhoton))
  {
    if (CheckTagBit(IsoGamma.MCTag, kMCPrompt) || CheckTagBit(IsoGamma.MCTag, kMCFragmentation) && pTisoGen < pTisoCorrectedMax)
    {
      return true;
    }
  }
  return false;
}

void doIsoGammaClusterCuts(std::vector<Cluster> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<Cluster>::iterator iter;
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

void doIsoGammaShowerShapeCuts(std::vector<Cluster> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<Cluster>::iterator iter;
  for (iter = IsoGammas.begin(); iter != IsoGammas.end();)
  {
    if (!IsoGammaCuts.PassedShowerShapeCuts(*iter))
    {
      iter = IsoGammas.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

void doIsoGammaCuts(std::vector<Cluster> &IsoGammas, IsoGammaCuts IsoGammaCuts)
{
  std::vector<Cluster>::iterator iter;
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
class DLJetCuts
{
private:
  float PtMin = 0.;
  float AreaMin = 0.;
public:
  DLJetCuts(GlobalOptions optns);
  ~DLJetCuts() {};
  bool PassedCuts(Jet Jet);
  float R = 0.4;
  bool doUEsubtraction = false;
  string UEEstimationMethod = "JetArea";
  float JetEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
};

DLJetCuts::DLJetCuts(GlobalOptions optns)
{
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];

  // Acceptance cuts
  R = chosencut["jet_Radius"].IsDefined() ? chosencut["jet_Radius"].as<float>() : standardcut["jet_Radius"].as<float>();

  TString JetAcceptance = (TString)(chosencut["jet_Acceptance"].IsDefined() ? chosencut["jet_Acceptance"].as<std::string>() : standardcut["jet_Acceptance"].as<std::string>()).c_str();
  SetAcceptance(JetAcceptance, JetEtaPhiMinMax, R);

  PtMin = chosencut["jetMinPt"].IsDefined() ? chosencut["jetMinPt"].as<float>() : standardcut["jetMinPt"].as<float>();

  AreaMin = chosencut["jet_MinArea"].IsDefined() ? chosencut["jet_MinArea"].as<float>() : standardcut["jet_MinArea"].as<float>();

  doUEsubtraction = chosencut["jet_applyUESubtraction"].IsDefined() ? chosencut["jet_applyUESubtraction"].as<bool>() : standardcut["jet_applyUESubtraction"].as<bool>();
  
  UEEstimationMethod = chosencut["jet_UEEstimationMethod"].IsDefined() ? chosencut["jet_UEEstimationMethod"].as<std::string>() : standardcut["jet_UEEstimationMethod"].as<std::string>();

  LOG(Form("%s JetCuts: pT > %.1f GeV/c", optns.cutString.Data(), PtMin))
  LOG(Form("Jet %s Acceptance: %.2f < Jet_Eta < %.2f | %.2f < Jet_Phi < %.2f", JetAcceptance.Data(), JetEtaPhiMinMax[0][0], JetEtaPhiMinMax[0][1], JetEtaPhiMinMax[1][0], JetEtaPhiMinMax[1][1]))
  LOG(Form("Jet %s Area > %.2f", JetAcceptance.Data(), AreaMin))
}

bool DLJetCuts::PassedCuts(Jet Jet)
{
  bool passed = true;
  if (Jet.Pt() < PtMin)
  {
    passed = false;
  }
  // Check GammaGen acceptance
  if (!Jet.isInJetAcceptance(JetEtaPhiMinMax))
  {
    passed = false;
  }
  if (Jet.Area < AreaMin)
  {
    passed = false;
  }

  return passed;
}

void doJetCuts(std::vector<DLJet> &Jets, DLJetCuts JetCuts)
{
  std::vector<DLJet>::iterator iter;
  for (iter = Jets.begin(); iter != Jets.end();)
  {
    if (!JetCuts.PassedCuts(*iter))
    {
      iter = Jets.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
// -----------Exclusive Trigger Particle Selection-----------
//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
// this class handles many things related to the selection
// of exclusive trigger particles, i.e. determination of signal
// and reference trigger particle (only one random per event in specified pt range)
class ExclusiveTriggerParticleSelection
{
private:
  float PtMinSignal = 0.;
  float PtMaxSignal = 0.;
  float PtMinReference = 0.;
  float PtMaxReference = 0.;
  bool isSigEvent = false;
  float SignalFraction = 0.9;
  TRandom3 rndm;
public:
  bool doExclusiveSelections = false;
  ExclusiveTriggerParticleSelection(GlobalOptions optns);
  ~ExclusiveTriggerParticleSelection() {};
  bool isSignalEvent();
  template <typename T>
  bool getTriggerParticle(std::vector<T> Objects, T &trigObject);
  float getSignalFraction() {return SignalFraction;};
  float getPtMinSignal() {return PtMinSignal;};
  float getPtMaxSignal() {return PtMaxSignal;};
  float getPtMinReference() {return PtMinReference;};
  float getPtMaxReference() {return PtMaxReference;};
};

ExclusiveTriggerParticleSelection::ExclusiveTriggerParticleSelection(GlobalOptions optns)
{
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");
  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];

  PtMinSignal = chosencut["trigger_sig_ptmin"].IsDefined() ? chosencut["trigger_sig_ptmin"].as<float>() : standardcut["trigger_sig_ptmin"].as<float>();
  PtMaxSignal = chosencut["trigger_sig_ptmax"].IsDefined() ? chosencut["trigger_sig_ptmax"].as<float>() : standardcut["trigger_sig_ptmax"].as<float>();
  PtMinReference = chosencut["trigger_ref_ptmin"].IsDefined() ? chosencut["trigger_ref_ptmin"].as<float>() : standardcut["trigger_ref_ptmin"].as<float>();
  PtMaxReference = chosencut["trigger_ref_ptmax"].IsDefined() ? chosencut["trigger_ref_ptmax"].as<float>() : standardcut["trigger_ref_ptmax"].as<float>();
  SignalFraction = chosencut["event_signal_fraction"].IsDefined() ? chosencut["event_signal_fraction"].as<float>() : standardcut["event_signal_fraction"].as<float>();
  doExclusiveSelections = chosencut["doExclusiveSelections"].IsDefined() ? chosencut["doExclusiveSelections"].as<bool>() : standardcut["doExclusiveSelections"].as<bool>();
}
bool ExclusiveTriggerParticleSelection::isSignalEvent()
{
  if (rndm.Rndm() < SignalFraction)
  {
    isSigEvent = true;
  } else {
    isSigEvent = false;
  }
  return isSigEvent;
}

// returns signal or reference trigger particle
template <typename T>
bool ExclusiveTriggerParticleSelection::getTriggerParticle(std::vector<T> Objects, T &trigObject)
{
  std::vector<T> triggerObjects;
  if(isSigEvent){ // signal trigger particle
     for (auto &obj : Objects)
     {
        if (obj.Pt() > PtMinSignal && obj.Pt() < PtMaxSignal)
        {
          triggerObjects.push_back(obj);
        }
     }
  } else { // reference trigger particle
    // find all objs in vector that pass PtMinReference < pt < PtMaxReference
    for (auto &obj : Objects)
    {
      if (obj.Pt() > PtMinReference && obj.Pt() < PtMaxReference)
      {
        triggerObjects.push_back(obj);
      }
    }
  }
  // check if triggerObjects vector is empty
  if (triggerObjects.empty())
  {
    return false;
  }

  // select one random object from the vector of trigger objects
  int index = rndm.Integer(triggerObjects.size());
  trigObject = triggerObjects[index];
  return true;
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
  float MinvMax = 0.17;
  float MinvMin = 0.1;

public:
  Pi0Cuts(GlobalOptions optns);
  ~Pi0Cuts() {};
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
  if (Pi0.Mass < MinvMin || Pi0.Mass > MinvMax)
  {
    passed = false;
  }

  return passed;
}

void doPi0Cuts(std::vector<Pi0> &Pi0s, Pi0Cuts pi0Cuts)
{
  std::vector<Pi0>::iterator iter;
  for (iter = Pi0s.begin(); iter != Pi0s.end();)
  {
    if (!pi0Cuts.PassedCuts(*iter))
    {
      iter = Pi0s.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}


#endif