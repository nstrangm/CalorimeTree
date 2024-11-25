#ifndef _writeanalysistree_h_included
#define _writeanalysistree_h_included

#include <yaml-cpp/yaml.h>

// Histograms for QA purposes
TList *outputhists;

TH1F *hTrackPt;
TH1F *hClusterEnergy;
TH1F *hClusterEta;
TH1F *hClusterPhi;
TH1F *hNEvents;
TH1F *hEvtVtxX;
TH1F *hEvtVtxY;
TH1F *hEvtVtxZ;

TH2F *hClusterM02vsE;

// define global switches
Bool_t fIsPbPb = kFALSE;
Bool_t fIsMC = kFALSE;

//----
// define output tree
//----
TTree *outputTree;

// TODO: MC truth information

// cuts for tree production
YAML::Node treecuts;
YAML::Node eventCuts;
YAML::Node trackCuts;

#endif