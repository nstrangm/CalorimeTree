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