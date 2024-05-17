#include "../Plotting/PlottingClass.h"

enum mcTypes { kPhoton,
               kPrompt,
               kFragmentation,
               kISR,
               kPi0Decay,
               kEtaDecay,
               kOtherDecay,
               kDecayPairLost,
               kDecayPairInCalo,
               kDecayDalitz,
               kConversion,
               kElectron,
               kEFromCFromB,
               kEFromC,
               kEFromB,
               kZDecay,
               kWDecay,
               kMuon,
               kPion,
               kPi0,
               kKaon,
               kEta,
               kProton,
               kAntiProton,
               kNeutron,
               kAntiNeutron,
               kUnknown,
               kBadLabel };

bool CheckTagBit(int tag, int test)
{
  if (tag < 0)
    return kFALSE;
  if (tag & (1 << test))
    return kTRUE;
  else
    return kFALSE;
}

void PrintMCTag(Int_t tag)
{
  INFO(Form(
    "Tag %d: photon %d, conv %d, prompt %d, frag %d, isr %d, pi0 decay %d, eta decay %d, other decay %d, lost decay %d, in calo decay %d,  pi0 %d,  eta %d, electron %d, muon %d,pion %d, proton %d, neutron %d, kaon %d, a-proton %d, a-neutron %d, unk %d, bad %d",
    tag,
    CheckTagBit(tag, kPhoton),
    CheckTagBit(tag, kConversion),
    CheckTagBit(tag, kPrompt),
    CheckTagBit(tag, kFragmentation),
    CheckTagBit(tag, kISR),
    CheckTagBit(tag, kPi0Decay),
    CheckTagBit(tag, kEtaDecay),
    CheckTagBit(tag, kOtherDecay),
    CheckTagBit(tag, kDecayPairLost),
    CheckTagBit(tag, kDecayPairInCalo),
    CheckTagBit(tag, kPi0),
    CheckTagBit(tag, kEta),
    CheckTagBit(tag, kElectron),
    CheckTagBit(tag, kMuon),
    CheckTagBit(tag, kPion),
    CheckTagBit(tag, kProton),
    CheckTagBit(tag, kAntiNeutron),
    CheckTagBit(tag, kKaon),
    CheckTagBit(tag, kAntiProton),
    CheckTagBit(tag, kAntiNeutron),
    CheckTagBit(tag, kUnknown),
    CheckTagBit(tag, kBadLabel)))
}
