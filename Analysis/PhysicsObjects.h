#ifndef _PHYSICSOBJECTS_
#define _PHYSICSOBJECTS_

#include "TreeUtilities.h"
#include "Logging.h"
#include "MCUtils.h"
#include "Geometry.h"

class PhysicsObject
{
public:
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float px = 0;
  float py = 0;
  float pz;
  float Px() { return px; }
  float Py() { return py; }
  float Pz() { return pz; }
  float P() { return TMath::Sqrt(px * px + py * py + pz * pz); }
  float Pt() { return pt; }
  float Eta() { return eta; }
  float Phi() { return phi; }
  PhysicsObject() {};
  // explicit PhysicsObject(float PxFromTree, float PyFromTree, float PzFromTree); // explicit means this constructor can be inherited by the Gammas, Jets and Pi0s
  explicit PhysicsObject(float PtFromTree, float EtaFromTree, float PhiFromTree); // explicit means this constructor can be inherited by the Gammas, Jets and Pi0s
  ~PhysicsObject() {};
};

PhysicsObject::PhysicsObject(float PtFromTree, float EtaFromTree, float PhiFromTree)
{
  pt = PtFromTree;
  eta = EtaFromTree;
  phi = PhiFromTree;
  px = PxFromPtPhi(pt, phi);
  py = PyFromPtPhi(pt, phi);
  pz = PzFromPtEta(pt, eta);
}

template <typename T1, typename T2>
float calculateAngleBetweenPhysicsObjects(T1 o1, T2 o2)
{
  return TMath::ACos((o1.Px() * o2.Px() + o1.Py() * o2.Py() + o1.Pz() * o2.Pz()) / (PFromPxPyPz(o1.Px(), o1.Py(), o1.Pz()) * PFromPxPyPz(o2.Px(), o2.Py(), o2.Pz())));
}

template <typename T1, typename T2>
float calculateDeltaPhiBetweenPhysicsObjects(T1 o1, T2 o2)
{
  float dPhi = fabs(o1.Phi() - o2.Phi());
  if (dPhi > TMath::Pi())
    dPhi = fabs(dPhi - 2 * TMath::Pi());
  return dPhi;
}

template <typename T1, typename T2>
float calculateDeltaEtaBetweenPhysicsObjects(T1 o1, T2 o2)
{
  return fabs(o1.Eta() - o2.Eta());
}

template <typename T1, typename T2>
float calculateRelativisticAngleBetweenPhysicsObjects(T1 o1, T2 o2)
{
  float dEta = o1.Eta() - o2.Eta();
  float dPhi = o1.Phi() - o2.Phi();
  if (dPhi > TMath::Pi())
    dPhi = fabs(dPhi - 2 * TMath::Pi());
  return TMath::Sqrt(dEta * dEta + dPhi * dPhi);
}

class Event
{
public:
  Event() {};
  Event(TreeBuffer tree, GlobalOptions optns);
  ~Event() {};

  float Rho = 0;

  int Multiplicity = 0;
  float Centrality = 0;
  uint16_t Selection = 0;
  uint32_t Alias = 0;

  unsigned short NPrimaryTracks = 0;
  bool IsTriggered = true;
  float ZVtx = 0;
  bool ZVtxOK = true;
  unsigned short Quality = 0;
  bool QualityOK = true;
  unsigned short NotAccepted = 0;

  bool Selected = false;

  // MC properties
  float weight = 1.;
};

Event::Event(TreeBuffer tree, GlobalOptions optns)
{
  Rho = tree.Event_Rho;

  switch (optns.TreeFormat)
  {
  case kRun2Tree:
    NPrimaryTracks = tree.Event_NPrimaryTracks;
    IsTriggered = tree.Event_IsTriggered;
    ZVtx = tree.Event_ZVertex;
    Quality = tree.Event_Quality;
    NotAccepted = tree.Event_NotAccepted;
    if (optns.isMC)
    {
      weight = tree.Event_Weight;
    }
    break;
  case kRun3Tree:
    Multiplicity = tree.Event_Multiplicity;
    Centrality = tree.Event_Centrality;
    Selection = tree.Event_Selection;
    Alias = tree.Event_Alias;
    break;
  case kBerkeleyTree:
  default:
    FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
  }
}

class TrackMatchedToIsoGamma
{
public:
  float P = 0;
  float Pt = 0;
  float dEta = 0;
  float dPhi = 0;
  bool IsConv = false;
  TrackMatchedToIsoGamma() {};
  void SetMatchedTrackProperties(float P_in, float Pt_in, float dEta_in, float dPhi_in, bool IsConv_in)
  {
    P = P_in;
    Pt = Pt_in;
    dEta = dEta_in;
    dPhi = dPhi_in;
    IsConv = IsConv_in;
  };
  ~TrackMatchedToIsoGamma() {};
};

class GammaGen : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~GammaGen() {};
  // bool isSignal();
  float E = 0;
  float IsoCharged = 0;
  float IsoBckPerp = 0;
  float IsoChargedCorrected = 0;
  int MCTag;
  bool isInEMCalAcceptance(float EMCalEtaPhiMinMax[2][2]);
  bool isInDCalAcceptance(float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2]);
};

bool GammaGen::isInEMCalAcceptance(float EMCalEtaPhiMinMax[2][2])
{
  return (Eta() < EMCalEtaPhiMinMax[0][1] && Eta() > EMCalEtaPhiMinMax[0][0] && Phi() > EMCalEtaPhiMinMax[1][0] && Phi() < EMCalEtaPhiMinMax[1][1]);
}

bool GammaGen::isInDCalAcceptance(float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2])
{
  if (Eta() < DCalEtaPhiMinMax[0][1] && Eta() > DCalEtaPhiMinMax[0][0] && Phi() > DCalEtaPhiMinMax[1][0] && Phi() < DCalEtaPhiMinMax[1][1])
  {                                                                                                                                                           // In DCal
    if (Eta() > DCalHoleEtaPhiMinMax[0][0] && Eta() < DCalHoleEtaPhiMinMax[0][1] && Phi() < DCalHoleEtaPhiMinMax[1][0] && Phi() > DCalHoleEtaPhiMinMax[1][1]) // In DCal hole
      return false;
    else // Not in DCal hole
      return true;
  }
  return false;
}

class Cluster : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~Cluster() {};
  bool isPhoton();
  // bool isSignal();
  bool isGammaFromPion();
  bool isGammaFromEta();
  bool isPion();
  bool isInEMCalAcceptance(float EMCalEtaPhiMinMax[2][2]);
  bool isInDCalAcceptance(float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2]);
  float E = 0;
  float EBeforeNL = 0; // To monitor the application of the NL
  float M02 = 0;
  float M02Recalc = 0; // Calculated in 5x5 window - not used for now
  float M20 = 0;
  unsigned short NCells = 0;
  float V1SplitMass = 0;
  float MinMassDiffToPi0 = 0; // ToDo: Problem: this is calculated in AliPhysics before applying the NL
  float MinMassDiffToEta = 0;
  unsigned short SM = 0; // EMCal/DCal super module number
  float EFrac = 0;
  float IsoCharged = 0;
  float IsoChargedCorrected = 0;
  float IsoBckPerp = 0;
  float TrueClusterIsoCharged = 0;
  float TrueClusterIsoBckPerp = 0;
  TrackMatchedToIsoGamma MatchedTrack;
  float DistanceToBadChannel = 1;
  float NLM = 0;
  float EventWeight = 1;

  int MCTag;
};

bool Cluster::isPhoton()
{
  if (CheckTagBit(MCTag, kMCPhoton))
  {
    return true;
  }
  return false;
}

bool Cluster::isGammaFromPion()
{
  if (CheckTagBit(MCTag, kMCPhoton))
  {
    if (CheckTagBit(MCTag, kMCPi0Decay))
    {
      return true;
    }
  }
  return false;
}

bool Cluster::isGammaFromEta()
{
  if (CheckTagBit(MCTag, kMCPhoton))
  {
    if (CheckTagBit(MCTag, kMCEtaDecay))
    {
      return true;
    }
  }
  return false;
}

bool Cluster::isPion()
{
  if (CheckTagBit(MCTag, kMCPion))
  {
    return true;
  }
  return false;
}

bool Cluster::isInEMCalAcceptance(float EMCalEtaPhiMinMax[2][2])
{
  return (Eta() < EMCalEtaPhiMinMax[0][1] && Eta() > EMCalEtaPhiMinMax[0][0] && Phi() > EMCalEtaPhiMinMax[1][0] && Phi() < EMCalEtaPhiMinMax[1][1]);
}

bool Cluster::isInDCalAcceptance(float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2])
{
  if (Eta() < DCalEtaPhiMinMax[0][1] && Eta() > DCalEtaPhiMinMax[0][0] && Phi() > DCalEtaPhiMinMax[1][0] && Phi() < DCalEtaPhiMinMax[1][1])
  {                                                                                                                                                           // In DCal
    if (Eta() > DCalHoleEtaPhiMinMax[0][0] && Eta() < DCalHoleEtaPhiMinMax[0][1] && Phi() < DCalHoleEtaPhiMinMax[1][0] && Phi() > DCalHoleEtaPhiMinMax[1][1]) // In DCal hole
      return false;
    else // Not in DCal hole
      return true;
  }
  return false;
}

void saveGenPhotonsFromEventInVector(TreeBuffer tree, std::vector<GammaGen> &GammaGens, GlobalOptions optns)
{
  for (int iGammaGen = 0; iGammaGen < (int)tree.GenPhoton_Px->size(); iGammaGen++)
  {
    if (optns.TreeFormat == kRun2Tree)
    {
      GammaGen gammaGen(PtFromPxPy(tree.GenPhoton_Px->at(iGammaGen), tree.GenPhoton_Py->at(iGammaGen)), EtaFromPxPyPz(tree.GenPhoton_Px->at(iGammaGen), tree.GenPhoton_Py->at(iGammaGen), tree.GenPhoton_Pz->at(iGammaGen)), PhiFromPxPy(tree.GenPhoton_Px->at(iGammaGen), tree.GenPhoton_Py->at(iGammaGen)));
      gammaGen.E = tree.GenPhoton_E->at(iGammaGen);
      gammaGen.IsoCharged = tree.GenPhoton_IsoCharged3->at(iGammaGen);
      gammaGen.IsoBckPerp = tree.GenPhoton_IsoBckPerp->at(iGammaGen);
      gammaGen.MCTag = tree.GenPhoton_MCTag->at(iGammaGen);
      GammaGens.push_back(gammaGen);
    }
    else if (optns.TreeFormat == kRun3Tree)
    {
      GammaGen gammaGen(0, 0, 0);
      // gammaGen.E = tree.GenPhoton_E->at(iGammaGen);
      // gammaGen.IsoCharged = tree.GenPhoton_IsoCharged3->at(iGammaGen);
      // gammaGen.IsoBckPerp = tree.GenPhoton_IsoBckPerp->at(iGammaGen);
      // gammaGen.MCTag = tree.GenPhoton_MCTag->at(iGammaGen);
      GammaGens.push_back(gammaGen);
    }
    else
    {
      FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
    }
  }
}

void saveClustersFromEventInVector(TreeBuffer tree, std::vector<Cluster> &IsoGammas, GlobalOptions optns)
{
  for (int iCluster = 0; iCluster < (int)tree.Cluster_E->size(); iCluster++)
  {
    if (optns.TreeFormat == kRun2Tree)
    {
      Cluster isoGamma(PtFromPxPy(tree.Cluster_Px->at(iCluster), tree.Cluster_Py->at(iCluster)), EtaFromPxPyPz(tree.Cluster_Px->at(iCluster), tree.Cluster_Py->at(iCluster), tree.Cluster_Pz->at(iCluster)), PhiFromPxPy(tree.Cluster_Px->at(iCluster), tree.Cluster_Py->at(iCluster)));
      isoGamma.EBeforeNL = tree.Cluster_E->at(iCluster);
      isoGamma.E = tree.Cluster_E->at(iCluster);
      isoGamma.M02 = tree.Cluster_M02->at(iCluster);
      isoGamma.M20 = tree.Cluster_M20->at(iCluster);
      isoGamma.NCells = tree.Cluster_NCells->at(iCluster);
      isoGamma.V1SplitMass = tree.Cluster_V1SplitMass->at(iCluster);
      isoGamma.MinMassDiffToPi0 = tree.Cluster_MinMassDiffToPi0->at(iCluster);
      isoGamma.MinMassDiffToEta = tree.Cluster_MinMassDiffToEta->at(iCluster);
      isoGamma.SM = tree.Cluster_SM->at(iCluster);
      isoGamma.EFrac = tree.Cluster_EFrac->at(iCluster);
      isoGamma.IsoCharged = tree.Cluster_IsoCharged3->at(iCluster);
      isoGamma.IsoBckPerp = tree.Cluster_IsoBckPerp->at(iCluster);
      isoGamma.MatchedTrack.SetMatchedTrackProperties(tree.Cluster_MatchTrackP->at(iCluster), tree.Cluster_MatchTrackPt->at(iCluster), tree.Cluster_MatchTrackdEta->at(iCluster), tree.Cluster_MatchTrackdPhi->at(iCluster), tree.Cluster_MatchTrackIsConv->at(iCluster));
      isoGamma.DistanceToBadChannel = tree.Cluster_DistanceToBadChannel->at(iCluster);
      isoGamma.NLM = tree.Cluster_NLM->at(iCluster);

      if (optns.isMC)
      {
        isoGamma.MCTag = tree.TrueCluster_MCTag->at(iCluster);
        isoGamma.EventWeight = tree.Event_Weight;
        isoGamma.TrueClusterIsoCharged = tree.TrueCluster_IsoCharged3->at(iCluster);
        isoGamma.TrueClusterIsoBckPerp = tree.TrueCluster_IsoBckPerp->at(iCluster);
      }

      IsoGammas.push_back(isoGamma);
    }
    else if (optns.TreeFormat == kRun3Tree)
    {
      Cluster isoGamma(PtFromPEta(tree.Cluster_E->at(iCluster), tree.Cluster_Eta->at(iCluster)), tree.Cluster_Eta->at(iCluster), tree.Cluster_Phi->at(iCluster));
      isoGamma.EBeforeNL = tree.Cluster_E->at(iCluster);
      isoGamma.E = tree.Cluster_E->at(iCluster);
      isoGamma.M02 = tree.Cluster_M02->at(iCluster);
      isoGamma.M20 = tree.Cluster_M20->at(iCluster);
      isoGamma.NCells = tree.Cluster_NCells->at(iCluster);
      isoGamma.IsoCharged = tree.Cluster_IsoCharged3->at(iCluster);
      isoGamma.IsoBckPerp = tree.Cluster_IsoBckPerp->at(iCluster);
      isoGamma.MatchedTrack.SetMatchedTrackProperties(0., 0., 0., 0., 0.);
      isoGamma.NLM = tree.Cluster_NLM->at(iCluster);
      if (optns.isMC)
      {
        FATAL("MC information for run3 clusters is not yet implemented!")
      }

      IsoGammas.push_back(isoGamma);
    }
    else
    {
      FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
    }
  }
}

// ########################################################################################################################
// ###### ToDo: Check whether this is run3 and if so do not apply a correction since this is done online already! #########
// ########################################################################################################################
void calculateIsolation(std::vector<Cluster> &IsoGammas, Event &Event, bool useRhoInsteadOfPerpCone)
{
  for (int iGamma = 0; iGamma < (int)IsoGammas.size(); iGamma++)
  {
    Cluster *isoGamma = &IsoGammas.at(iGamma);

    // Calculate corrected isolation pT subtracting backperp mult by cone area.
    float IsoChargedAcceptanceCorrected = isoGamma->IsoCharged / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4) - isoGamma->IsoBckPerp * TMath::Pi() * 0.4 * 0.4;

    if (useRhoInsteadOfPerpCone)
    {
      float IsoBckPerpAcceptanceCorrected = isoGamma->IsoCharged / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4) - Event.Rho * TMath::Pi() * 0.4 * 0.4;
    }
    else
    {
    }
    isoGamma->IsoChargedCorrected = IsoChargedAcceptanceCorrected;
  }
}

void GammaGencalculateIsolation(std::vector<GammaGen> &GammaGens, Event &Event)
{
  for (int iGamma = 0; iGamma < (int)GammaGens.size(); iGamma++)
  {
    GammaGen *GammaGen = &GammaGens.at(iGamma);
    float IsoChargedCorrected = GammaGen->IsoCharged - GammaGen->IsoBckPerp;
    GammaGen->IsoChargedCorrected = IsoChargedCorrected;
  }
}

class Jet : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~Jet() {};
  float Area = 0;
  unsigned short Nch = 0;
  unsigned short Nclus = 0;
  unsigned short Nconstits = 0;
  float M = 0;
  bool isInJetAcceptance(float JetEtaPhiMinMax[2][2]);
};

bool Jet::isInJetAcceptance(float JetEtaPhiMinMax[2][2])
{
  return (Eta() < JetEtaPhiMinMax[0][1] && Eta() > JetEtaPhiMinMax[0][0] && Phi() > JetEtaPhiMinMax[1][0] && Phi() < JetEtaPhiMinMax[1][1]);
}

class DLJet : public Jet
{
public:
  using Jet::Jet;
  ~DLJet() {};
  float Area = 0;
  unsigned short Nch = 0;
  unsigned short Nclus = 0;
  unsigned short Nconstits = 0;
  float M = 0;
  bool isInJetAcceptance(float JetEtaPhiMinMax[2][2]);
};

void saveJetsFromEventInVector(TreeBuffer tree, std::vector<DLJet> &Jets, GlobalOptions optns)
{
  for (int iJet = 0; iJet < (int)tree.Jet_Area->size(); iJet++)
  {
    if (optns.TreeFormat == kRun2Tree)
    {
      float px = tree.Jet_Px->at(iJet);
      float py = tree.Jet_Py->at(iJet);
      float pz = tree.Jet_Pz->at(iJet);
      DLJet jet(PtFromPxPy(px, py), EtaFromPxPyPz(px, py, pz), PhiFromPxPy(px, py));
      jet.Area = tree.Jet_Area->at(iJet);
      jet.Nch = tree.Jet_Nch->at(iJet);
      jet.Nclus = tree.Jet_Nclus->at(iJet);
      Jets.push_back(jet);
    }
    else if (optns.TreeFormat == kRun3Tree)
    {
      DLJet jet(tree.Jet_Pt->at(iJet), tree.Jet_Eta->at(iJet), tree.Jet_Phi->at(iJet));
      jet.Area = tree.Jet_Area->at(iJet);
      jet.M = tree.Jet_M->at(iJet);
      jet.Nconstits = tree.Jet_Nconstits->at(iJet);
      Jets.push_back(jet);
    }
    else
    {
      FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
    }
  }
}

class PLJet : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~PLJet() {};
  float Area = 0;
  unsigned short NPart = 0;
  Jet *ClosestDLJet = nullptr;
};

void savePLJetsFromEventInVector(TreeBuffer tree, std::vector<PLJet> &PLJets, GlobalOptions optns)
{
  for (int iJet = 0; iJet < (int)tree.PLJet_Px->size(); iJet++)
  {
    if (optns.TreeFormat == kRun2Tree)
    {
      float px = tree.PLJet_Px->at(iJet);
      float py = tree.PLJet_Py->at(iJet);
      float pz = tree.PLJet_Pz->at(iJet);
      PLJet jet(PtFromPxPy(px, py), EtaFromPxPyPz(px, py, pz), PhiFromPxPy(px, py));
      jet.Area = tree.PLJet_Area->at(iJet);
      jet.NPart = tree.PLJet_NPart->at(iJet);
      PLJets.push_back(jet);
      break;
    }
    else if (optns.TreeFormat == kRun3Tree)
    {
      // PLJet jet(tree.PLJet_Pt->at(iJet), tree.PLJet_Eta->at(iJet), tree.PLJet_Phi->at(iJet));
      // jet.Area = tree.PLJet_Area->at(iJet);
      // jet.NPart = tree.PLJet_NPart->at(iJet);
      // PLJets.push_back(jet);
    }
    else
    {
      FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
    }
  }
}

template <typename JT1, typename JT2>
short int getClosestJetNumber(JT1 Jet, std::vector<JT2> JetGroup, float R)
{
  float angleBetweenClosestJets = TMath::Pi();
  int closestJetNumber = -1;
  for (int iJetFromGroup = 0; iJetFromGroup < (int)JetGroup.size(); iJetFromGroup++)
  {
    float angleBetweenJets = calculateRelativisticAngleBetweenPhysicsObjects(Jet, JetGroup.at(iJetFromGroup));
    if (angleBetweenJets < angleBetweenClosestJets)
    {
      angleBetweenClosestJets = angleBetweenJets;
      closestJetNumber = iJetFromGroup;
    }
  }

  if (angleBetweenClosestJets > 2 * R) // If the distance is larger than 2xR, then no match can be made
    return -1;
  return closestJetNumber;
}

void mapPLtoDLjets(std::vector<DLJet> &DLJets, std::vector<PLJet> &PLJets, float R)
{
  std::vector<short int> closestPLJetFromDLJetNumber;
  std::vector<short int> closestDLJetFromPLJetNumber;
  for (int iDLJet = 0; iDLJet < (int)DLJets.size(); iDLJet++)
    closestPLJetFromDLJetNumber.push_back(getClosestJetNumber(DLJets.at(iDLJet), PLJets, R));
  for (int iPLJet = 0; iPLJet < (int)PLJets.size(); iPLJet++)
  {
    closestDLJetFromPLJetNumber.push_back(getClosestJetNumber(PLJets.at(iPLJet), DLJets, R));
    if (closestDLJetFromPLJetNumber.at(iPLJet) == -1) // No DL jet was found for this PL jet
      continue;
    if (closestPLJetFromDLJetNumber.at(closestDLJetFromPLJetNumber.at(iPLJet)) == iPLJet)
    {
      PLJets.at(iPLJet).ClosestDLJet = &DLJets.at(closestDLJetFromPLJetNumber.at(iPLJet));
    }
  }
}

class Pi0 : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~Pi0() {};
  float Mass;
};

class CorrelationPair
{
public:
  float DPhi = 0;
  float DEta = 0;
  float pTImbalance = 1; // pTjet/pTgamma
};

class GammaJetPair : public CorrelationPair
{
public:
  GammaJetPair(Cluster *isoGammaptr, Jet *jetptr);
  ~GammaJetPair() {};
  Cluster *isoGamma = nullptr;
  Jet *jet = nullptr;
};

GammaJetPair::GammaJetPair(Cluster *isoGammaptr, Jet *jetptr)
{
  isoGamma = isoGammaptr;
  jet = jetptr;
  DPhi = calculateDeltaPhiBetweenPhysicsObjects(*isoGamma, *jet);
  DEta = calculateDeltaEtaBetweenPhysicsObjects(*isoGamma, *jet);
  pTImbalance = jet->Pt() / isoGamma->Pt();
}

class GGPi0JetPair : public CorrelationPair
{
public:
  GGPi0JetPair(Pi0 *ggp0ptr, Jet *jetptr);
  ~GGPi0JetPair() {};
  Pi0 *ggPi0 = nullptr;
  Jet *jet = nullptr;
};

GGPi0JetPair::GGPi0JetPair(Pi0 *ggp0ptr, Jet *jetptr)
{
  ggPi0 = ggp0ptr;
  jet = jetptr;
  DPhi = calculateDeltaPhiBetweenPhysicsObjects(*ggPi0, *jet);
  DEta = calculateDeltaEtaBetweenPhysicsObjects(*ggPi0, *jet);
  pTImbalance = jet->Pt() / ggPi0->Pt();
}

template <typename X, typename Y, typename Z>
void pairXWithYIntoZ(std::vector<X> &x, std::vector<Y> &y, std::vector<Z> &z)
{
  for (int i = 0; i < (int)x.size(); i++)
  {
    for (int j = 0; j < (int)y.size(); j++)
    {
      Z pair(&x.at(i), &y.at(j));
      z.push_back(pair);
    }
  }
}

void pairGammasFromEventInVector(std::vector<Cluster> &IsoGammas, std::vector<Pi0> &Pi0s)
{
  for (int ig1 = 0; ig1 < (int)IsoGammas.size(); ig1++)
  {
    for (int ig2 = ig1 + 1; ig2 < (int)IsoGammas.size(); ig2++)
    {
      float pt = PtFromPxPy(IsoGammas.at(ig1).Px() + IsoGammas.at(ig2).Px(), IsoGammas.at(ig1).Py() + IsoGammas.at(ig2).Py());
      float eta = EtaFromPxPyPz(IsoGammas.at(ig1).Px() + IsoGammas.at(ig2).Px(), IsoGammas.at(ig1).Py() + IsoGammas.at(ig2).Py(), IsoGammas.at(ig1).Pz() + IsoGammas.at(ig2).Pz());
      float phi = PhiFromPxPy(IsoGammas.at(ig1).Px() + IsoGammas.at(ig2).Px(), IsoGammas.at(ig1).Py() + IsoGammas.at(ig2).Py());
      Pi0 pi0(pt, eta, phi);
      pi0.Mass = TMath::Sqrt(2 * IsoGammas.at(ig1).E * IsoGammas.at(ig2).E * (1 - TMath::Cos(calculateAngleBetweenPhysicsObjects(IsoGammas.at(ig1), IsoGammas.at(ig2)))));
      Pi0s.push_back(pi0);
    }
  }
}

#endif // _PHYSICSOBJECTS_