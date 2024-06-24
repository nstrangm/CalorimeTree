#ifndef _PHYSICSOBJECTS_
#define _PHYSICSOBJECTS_

#include "TreeUtilities.h"
#include "Logging.h"
#include "MCUtils.h"

class PhysicsObject
{
public:
  float Px = 0;
  float Py = 0;
  float Pz = 0;
  float P();
  float Pt();
  float Eta();
  float Phi();
  PhysicsObject(){};
  explicit PhysicsObject(float PxFromTree, float PyFromTree, float PzFromTree); // explicit means this constructor can be inherited by the Gammas, Jets and Pi0s
  ~PhysicsObject(){};
};

PhysicsObject::PhysicsObject(float PxFromTree, float PyFromTree, float PzFromTree)
{
  Px = PxFromTree;
  Py = PyFromTree;
  Pz = PzFromTree;
}

float PhysicsObject::P()
{
  return TMath::Sqrt(Px * Px + Py * Py + Pz * Pz);
}
float PhysicsObject::Pt()
{
  return TMath::Sqrt(Px * Px + Py * Py);
}
float PhysicsObject::Eta()
{
  return 0.5 * TMath::Log((P() + Pz) / (P() - Pz)); // https://en.wikipedia.org/wiki/Pseudorapidity
}
float PhysicsObject::Phi()
{
  float Phi = -1.;
  if (Px > 0)
  {
    if (Py > 0)
      Phi = TMath::ATan(Py / Px);
    else
      Phi = TMath::ATan(Py / Px) + 2 * TMath::Pi();
  }
  else
  {
    Phi = TMath::ATan(Py / Px) + TMath::Pi();
  }
  return Phi;
}

template <typename T1, typename T2>
float calculateAngleBetweenPhysicsObjects(T1 o1, T2 o2)
{
  return TMath::ACos((o1.Px * o2.Px + o1.Py * o2.Py + o1.Pz * o2.Pz) / (o1.P() * o2.P()));
}

template <typename T1, typename T2>
float calculateDeltaPhiBetweenPhysicsObjects(T1 o1, T2 o2)
{
  float dPhi = o1.Phi() - o2.Phi();
  if (dPhi > TMath::Pi())
    dPhi = fabs(dPhi - 2 * TMath::Pi());
  return dPhi;
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
  Event(){};
  Event(TreeBuffer tree, GlobalOptions optns);
  ~Event(){};

  float Rho = 0;
  unsigned short NPrimaryTracks = 0;
  bool IsTriggered = false;
  float ZVtx = 0;
  bool ZVtxOK = true;
  unsigned short Quality = 0;
  bool QualityOK = true;
  unsigned short NotAccepted = 0;

  bool Selected = false;

  // MC properties
  float weight = 1;
};

Event::Event(TreeBuffer tree, GlobalOptions optns)
{
  Rho = tree.Event_Rho;
  NPrimaryTracks = tree.Event_NPrimaryTracks;
  IsTriggered = tree.Event_IsTriggered;
  ZVtx = tree.Event_ZVertex;
  Quality = tree.Event_Quality;
  NotAccepted = tree.Event_NotAccepted;

  if (optns.isMC)
  {
    weight = tree.Event_Weight;
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
  TrackMatchedToIsoGamma(){};
  void SetMatchedTrackProperties(float P_in, float Pt_in, float dEta_in, float dPhi_in, bool IsConv_in)
  {
    P = P_in;
    Pt = Pt_in;
    dEta = dEta_in;
    dPhi = dPhi_in;
    IsConv = IsConv_in;
  };
  ~TrackMatchedToIsoGamma(){};
};

class GammaGen : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~GammaGen(){};
  float E =0;
  float IsoCharged=0;
  float IsoBckPerp=0;
  float MCTag=0;
};

class IsoGamma : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~IsoGamma(){};
  bool isPhoton();
  bool isSignal();
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
  TrackMatchedToIsoGamma MatchedTrack;
  float DistanceToBadChannel = 0;
  float NLM = 0;
  float EventWeight = 1;

  int MCTag;
};

bool IsoGamma::isPhoton()
{
  if (CheckTagBit(MCTag, kMCPhoton))
  {
    return true;
  }
  return false;
}

bool IsoGamma::isSignal()
{
  if (CheckTagBit(MCTag, kMCPhoton))
  {
    if (CheckTagBit(MCTag, kMCPrompt) || CheckTagBit(MCTag, kMCFragmentation))
    {
      return true;
    }
  }
  return false;
}

bool IsoGamma::isGammaFromPion()
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

bool IsoGamma::isGammaFromEta()
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

bool IsoGamma::isPion()
{
  if (CheckTagBit(MCTag, kMCPion))
  {
    return true;
  }
  return false;
}

bool IsoGamma::isInEMCalAcceptance(float EMCalEtaPhiMinMax[2][2]){
  return (Eta() < EMCalEtaPhiMinMax[0][1] && Eta() > EMCalEtaPhiMinMax[0][0] && Phi() > EMCalEtaPhiMinMax[1][0] && Phi() < EMCalEtaPhiMinMax[1][1]);
}

bool IsoGamma::isInDCalAcceptance(float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2]){
  if(Eta() < DCalEtaPhiMinMax[0][1] && Eta() > DCalEtaPhiMinMax[0][0] && Phi() > DCalEtaPhiMinMax[1][0] && Phi() < DCalEtaPhiMinMax[1][1]){ // In DCal
    if(Eta() > DCalHoleEtaPhiMinMax[0][0] && Eta() < DCalHoleEtaPhiMinMax[0][1] && Phi() < DCalHoleEtaPhiMinMax[1][0] && Phi() > DCalHoleEtaPhiMinMax[1][1]) // In DCal hole
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
    if (optns.isMC)
    {
      GammaGen gammaGen(tree.GenPhoton_Px->at(iGammaGen), tree.GenPhoton_Py->at(iGammaGen), tree.GenPhoton_Pz->at(iGammaGen));
      gammaGen.E = tree.GenPhoton_E->at(iGammaGen);
      gammaGen.IsoCharged = tree.GenPhoton_IsoCharged3->at(iGammaGen);
      gammaGen.IsoBckPerp = tree.GenPhoton_IsoBckPerp->at(iGammaGen);
      GammaGens.push_back(gammaGen);
    }
    
  }
}

void saveClustersFromEventInVector(TreeBuffer tree, std::vector<IsoGamma> &IsoGammas, GlobalOptions optns)
{
  for (int iCluster = 0; iCluster < (int)tree.Cluster_Px->size(); iCluster++)
  {
    IsoGamma isoGamma(tree.Cluster_Px->at(iCluster), tree.Cluster_Py->at(iCluster), tree.Cluster_Pz->at(iCluster));
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
    }

    IsoGammas.push_back(isoGamma);
  }
}

void calculateIsolation(std::vector<IsoGamma> &IsoGammas, Event &Event, bool useRhoInsteadOfPerpCone)
{
  for (int iGamma = 0; iGamma < (int)IsoGammas.size(); iGamma++)
  {
    IsoGamma *isoGamma = &IsoGammas.at(iGamma);
    // float IsoChargedAcceptanceCorrected = isoGamma->IsoCharged / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4);

    // Calculate corrected isolation pT subtracting backperp mult by cone area.
    float IsoChargedAcceptanceCorrected = isoGamma->IsoCharged / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4) - isoGamma->IsoBckPerp * TMath::Pi() * 0.4 * 0.4;

    if (useRhoInsteadOfPerpCone)
    {
      // cout<<"WRONG!"<<"\n";

      float IsoBckPerpAcceptanceCorrected = isoGamma->IsoCharged / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4) - Event.Rho * TMath::Pi() * 0.4 * 0.4;
      // isoGamma->IsoBckPerp / CalculateIsoCorrectionFactor(isoGamma->Eta(), 0.8, 0.4);
    }
    else
    {
    }
    isoGamma->IsoChargedCorrected = IsoChargedAcceptanceCorrected;
  }
}

class Jet : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~Jet(){};
  float Area = 0;
  unsigned short Nch = 0;
  unsigned short Nclus = 0;
};

void saveJetsFromEventInVector(TreeBuffer tree, std::vector<Jet> &Jets)
{
  for (int iJet = 0; iJet < (int)tree.Jet_Px->size(); iJet++)
  {
    Jet jet(tree.Jet_Px->at(iJet), tree.Jet_Py->at(iJet), tree.Jet_Pz->at(iJet));
    jet.Area = tree.Jet_Area->at(iJet);
    jet.Nch = tree.Jet_Nch->at(iJet);
    jet.Nclus = tree.Jet_Nclus->at(iJet);
    Jets.push_back(jet);
  }
}

class PLJet : public PhysicsObject
{
public:
  using PhysicsObject::PhysicsObject;
  ~PLJet(){};
  float Area = 0;
  unsigned short NPart = 0;
  Jet *ClosestDLJet = nullptr;
};

void savePLJetsFromEventInVector(TreeBuffer tree, std::vector<PLJet> &PLJets)
{
  for (int iJet = 0; iJet < (int)tree.PLJet_Px->size(); iJet++)
  {
    PLJet jet(tree.PLJet_Px->at(iJet), tree.PLJet_Py->at(iJet), tree.PLJet_Pz->at(iJet));
    jet.Area = tree.PLJet_Area->at(iJet);
    jet.NPart = tree.PLJet_NPart->at(iJet);
    PLJets.push_back(jet);
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

void mapPLtoDLjets(std::vector<Jet> &DLJets, std::vector<PLJet> &PLJets, float R)
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
  ~Pi0(){};
  float Mass;
};

class GammaJetPair
{
public:
  GammaJetPair(IsoGamma *isoGammaptr, Jet *jetptr);
  ~GammaJetPair(){};
  IsoGamma *isoGamma;
  Jet *jet;
  float DPhi = 0;
  float pTImbalance = 1; // pTjet/pTgamma
};

GammaJetPair::GammaJetPair(IsoGamma *isoGammaptr, Jet *jetptr)
{
  isoGamma = isoGammaptr;
  jet = jetptr;
  DPhi = calculateDeltaPhiBetweenPhysicsObjects(*isoGamma, *jet);
  pTImbalance = jet->Pt() / isoGamma->Pt();
}

void pairGammasWithJets(std::vector<IsoGamma> &IsoGammas, std::vector<Jet> &Jets, std::vector<GammaJetPair> &GammaJetPairs)
{
  for (int ig = 0; ig < (int)IsoGammas.size(); ig++)
  {
    for (int ij = 0; ij < (int)Jets.size(); ij++)
    {
      GammaJetPair gammaJetPair(&IsoGammas.at(ig), &Jets.at(ij));
      GammaJetPairs.push_back(gammaJetPair);
    }
  }
}

void pairIsoGammasFromEventInVector(std::vector<IsoGamma> &IsoGammas, std::vector<Pi0> &Pi0s)
{
  for (int ig1 = 0; ig1 < (int)IsoGammas.size(); ig1++)
  {
    for (int ig2 = ig1 + 1; ig2 < (int)IsoGammas.size(); ig2++)
    {
      Pi0 pi0(IsoGammas.at(ig1).Px + IsoGammas.at(ig2).Px, IsoGammas.at(ig1).Py + IsoGammas.at(ig2).Py, IsoGammas.at(ig1).Pz + IsoGammas.at(ig2).Pz);
      pi0.Mass = TMath::Sqrt(2 * IsoGammas.at(ig1).E * IsoGammas.at(ig2).E * (1 - TMath::Cos(calculateAngleBetweenPhysicsObjects(IsoGammas.at(ig1), IsoGammas.at(ig2)))));
      Pi0s.push_back(pi0);
    }
  }
}

#endif // _PHYSICSOBJECTS_