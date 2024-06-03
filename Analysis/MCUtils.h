#include "PlottingClass.h"

enum mcTypes { kMCPhoton,
               kMCPrompt,
               kMCFragmentation,
               kMCISR,
               kMCPi0Decay,
               kMCEtaDecay,
               kMCOtherDecay,
               kMCDecayPairLost,
               kMCDecayPairInCalo,
               kMCDecayDalitz,
               kMCConversion,
               kMCElectron,
               kMCEFromCFromB,
               kMCEFromC,
               kMCEFromB,
               kMCZDecay,
               kMCWDecay,
               kMCMuon,
               kMCPion,
               kMCPi0,
               kMCKaon,
               kMCEta,
               kMCProton,
               kMCAntiProton,
               kMCNeutron,
               kMCAntiNeutron,
               kMCUnknown,
               kMCBadLabel };

bool CheckTagBit(int tag, int test)
{
  if (tag < 0)
    return false;
  if (tag & (1 << test))
    return true;
  else
    return false;
}

void PrintMCTag(Int_t tag)
{
  INFO(Form(
    "Tag %d: photon %d, conv %d, prompt %d, frag %d, isr %d, pi0 decay %d, eta decay %d, other decay %d, lost decay %d, in calo decay %d,  pi0 %d,  eta %d, electron %d, muon %d,pion %d, proton %d, neutron %d, kMCaon %d, a-proton %d, a-neutron %d, unk %d, bad %d",
    tag,
    CheckTagBit(tag, kMCPhoton),
    CheckTagBit(tag, kMCConversion),
    CheckTagBit(tag, kMCPrompt),
    CheckTagBit(tag, kMCFragmentation),
    CheckTagBit(tag, kMCISR),
    CheckTagBit(tag, kMCPi0Decay),
    CheckTagBit(tag, kMCEtaDecay),
    CheckTagBit(tag, kMCOtherDecay),
    CheckTagBit(tag, kMCDecayPairLost),
    CheckTagBit(tag, kMCDecayPairInCalo),
    CheckTagBit(tag, kMCPi0),
    CheckTagBit(tag, kMCEta),
    CheckTagBit(tag, kMCElectron),
    CheckTagBit(tag, kMCMuon),
    CheckTagBit(tag, kMCPion),
    CheckTagBit(tag, kMCProton),
    CheckTagBit(tag, kMCAntiNeutron),
    CheckTagBit(tag, kMCKaon),
    CheckTagBit(tag, kMCAntiProton),
    CheckTagBit(tag, kMCAntiNeutron),
    CheckTagBit(tag, kMCUnknown),
    CheckTagBit(tag, kMCBadLabel)))
}
