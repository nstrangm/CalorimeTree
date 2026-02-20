#ifndef TRIGGERALIASES_H_
#define TRIGGERALIASES_H_

enum triggerAliases {
  kINT7 = 0,
  kEMC7,
  kINT7inMUON,
  kMuonSingleLowPt7,
  kMuonSingleHighPt7,
  kMuonUnlikeLowPt7,
  kMuonLikeLowPt7,
  kCUP8,
  kCUP9,
  kMUP10,
  kMUP11,
  kINT1,
  kUnbiased,
  kDMC7,
  kEG1,
  kEJ1,
  kEG2,
  kEJ2,
  kDG1,
  kDJ1,
  kDG2,
  kDJ2,
  kTVXinTRD,
  kTVXinEMC,
  kTVXinPHOS,
  kTVXinHMP,
  kPHOS,
  kALL,
  kNaliases
};

// std::map for a string for each trigger alias
std::map<triggerAliases, TString> triggerAliasStrings = {
  {kINT7, "INT7"},
    {kEMC7, "EMC7"},
    {kINT7inMUON, "INT7inMUON"},
    {kMuonSingleLowPt7, "MuonSingleLowPt7"},
    {kMuonSingleHighPt7, "MuonSingleHighPt7"},
    {kMuonUnlikeLowPt7, "MuonUnlikeLowPt7"},
    {kMuonLikeLowPt7, "MuonLikeLowPt7"},
    {kCUP8, "CUP8"},
    {kCUP9, "CUP9"},
    {kMUP10, "MUP10"},
    {kMUP11, "MUP11"},
    {kINT1, "INT1"},
    {kUnbiased, "Unbiased"},
    {kDMC7, "DMC7"},
    {kEG1, "EG1"},
    {kEJ1, "EJ1"},
    {kEG2, "EG2"},
    {kEJ2, "EJ2"},
    {kDG1, "DG1"},
    {kDJ1, "DJ1"},
    {kDG2, "DG2"},
    {kDJ2, "DJ2"},
    {kTVXinTRD, "TVXinTRD"},
    {kTVXinEMC, "TVXinEMC"},
    {kTVXinPHOS, "TVXinPHOS"},
    {kTVXinHMP, "TVXinHMP"},
    {kPHOS, "PHOS"},
    {kALL, "ALL"}
};

#endif