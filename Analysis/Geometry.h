#include "TMath.h"

float PFromPxPyPz(float Px, float Py, float Pz)
{
  return TMath::Sqrt(Px * Px + Py * Py + Pz * Pz);
}
float PtFromPxPy(float Px, float Py)
{
  return TMath::Sqrt(Px * Px + Py * Py);
}
float ThetaFromEta(float Eta){ return TMath::Pi()/2+2*TMath::ATan(Eta/fabs(Eta)*TMath::Exp(-fabs(Eta)));}
float PtFromPEta(float P, float Eta)
{
  return fabs(P*TMath::Cos(ThetaFromEta(Eta)));
}
float EtaFromPxPyPz(float Px, float Py, float Pz)
{
  return - 0.5 * TMath::Log((PFromPxPyPz(Px, Py, Pz) + Pz) / (PFromPxPyPz(Px, Py, Pz) - Pz)); // https://en.wikipedia.org/wiki/Pseudorapidity
}
float PhiFromPxPy(float Px, float Py)
{
  float Phi = -1.;
  if (Px >= 0)
  {
    if (Py >= 0)
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

float PxFromPtPhi(float Pt, float Phi){
    return Pt*TMath::Cos(Phi);
}

float PyFromPtPhi(float Pt, float Phi){
    return Pt*TMath::Sin(Phi);
}

float PzFromPtEta(float Pt, float Eta){
    return Pt*TMath::Tan(ThetaFromEta(Eta));
}




void SetAcceptance(TString ClusterAcceptance, float EMCalEtaPhiMinMax[2][2], float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2], float FiducialSpace = 0.){
  if(ClusterAcceptance == "NoCut"){
  // EMCal
  EMCalEtaPhiMinMax[0][0] = -1E3;
  EMCalEtaPhiMinMax[0][1] = 1E3;
  EMCalEtaPhiMinMax[1][0] = -1E3;
  EMCalEtaPhiMinMax[1][1] = 1E3;
  // DCal
  DCalEtaPhiMinMax[0][0] = -1E3;
  DCalEtaPhiMinMax[0][1] = 1E3;
  DCalEtaPhiMinMax[1][0] = -1E3;
  DCalEtaPhiMinMax[1][1] = 1E3;
  // DCal hole
  DCalHoleEtaPhiMinMax[0][0] = 0.;
  DCalHoleEtaPhiMinMax[0][1] = 0.;
  DCalHoleEtaPhiMinMax[1][0] = 0.;
  DCalHoleEtaPhiMinMax[1][1] = 0.;
  } else if(ClusterAcceptance == "EMC"){
  // EMCal
  EMCalEtaPhiMinMax[0][0] = -0.67;
  EMCalEtaPhiMinMax[0][1] = 0.67;
  EMCalEtaPhiMinMax[1][0] = 1.4;
  EMCalEtaPhiMinMax[1][1] = 3.28;
  // DCal
  DCalEtaPhiMinMax[0][0] = -0.67;
  DCalEtaPhiMinMax[0][1] = 0.67;
  DCalEtaPhiMinMax[1][0] = 4.57;
  DCalEtaPhiMinMax[1][1] = 5.70;
  // DCal hole
  DCalHoleEtaPhiMinMax[0][0] = -0.23;
  DCalHoleEtaPhiMinMax[0][1] = 0.23;
  DCalHoleEtaPhiMinMax[1][0] = 4.57;
  DCalHoleEtaPhiMinMax[1][1] = 5.58;
  } else if(ClusterAcceptance == "EMC_Fiducial"){
  // EMCal
  EMCalEtaPhiMinMax[0][0] = -0.67 + FiducialSpace;
  EMCalEtaPhiMinMax[0][1] = 0.67 - FiducialSpace;
  EMCalEtaPhiMinMax[1][0] = 1.4 + FiducialSpace;
  EMCalEtaPhiMinMax[1][1] = 3.28 - FiducialSpace;
  // DCal
  DCalEtaPhiMinMax[0][0] = -0.67 + FiducialSpace;
  DCalEtaPhiMinMax[0][1] = 0.67 - FiducialSpace;
  DCalEtaPhiMinMax[1][0] = 4.57 + FiducialSpace;
  DCalEtaPhiMinMax[1][1] = 5.70 - FiducialSpace;
  // DCal hole
  DCalHoleEtaPhiMinMax[0][0] = -0.23;
  DCalHoleEtaPhiMinMax[0][1] = 0.23;
  DCalHoleEtaPhiMinMax[1][0] = 4.57;
  DCalHoleEtaPhiMinMax[1][1] = 5.58;
  }
}

void SetAcceptance(TString JetAcceptance, float JetEtaPhiMinMax[2][2], float FiducialSpace = 0.){
  if(JetAcceptance == "TPC"){
    JetEtaPhiMinMax[0][0] = -0.9;
    JetEtaPhiMinMax[0][1] = 0.9;
    JetEtaPhiMinMax[1][0] = 0.;
    JetEtaPhiMinMax[1][1] = 2*3.1416;
  } else if(JetAcceptance == "TPC_Fiducial"){
    JetEtaPhiMinMax[0][0] = -0.9 + FiducialSpace;
    JetEtaPhiMinMax[0][1] = 0.9 - FiducialSpace;
    JetEtaPhiMinMax[1][0] = 0.;
    JetEtaPhiMinMax[1][1] = 2*3.1416;
  } else if(JetAcceptance == "NoCut"){
    JetEtaPhiMinMax[0][0] = -1E3;
    JetEtaPhiMinMax[0][1] = 1E3;
    JetEtaPhiMinMax[1][0] = -1E3;
    JetEtaPhiMinMax[1][1] = 1E3;
  // cut for testing purposes that severely restricts the jet acceptance to avoid detector effects
  // should be only used with R =0.2
  } else if(JetAcceptance == "TPC_Fiducial_Restricted"){
    JetEtaPhiMinMax[0][0] = -0.7;
    JetEtaPhiMinMax[0][1] = -0.5;
    JetEtaPhiMinMax[1][0] = 0.;
    JetEtaPhiMinMax[1][1] = 2*3.1416;
  }
}
