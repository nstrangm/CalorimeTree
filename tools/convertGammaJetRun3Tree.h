#ifndef _convertgammajetrun3tree_h_included
#define _convertgammajetrun3tree_h_included

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"




// Tree definitions
Float_t fBuffer_multiplicity;
Float_t fBuffer_centrality;
Float_t fBuffer_rho;
uint8_t fBuffer_eventselection;
uint32_t fBuffer_alias;
Int_t fBuffer_occupancy;
// jet
std::vector<Float_t> *fBuffer_jet_data_pt;
std::vector<Float_t> *fBuffer_jet_data_eta;
std::vector<Float_t> *fBuffer_jet_data_phi;
std::vector<Float_t> *fBuffer_jet_data_radius;
std::vector<Float_t> *fBuffer_jet_data_energy;
std::vector<Float_t> *fBuffer_jet_data_mass;
std::vector<Float_t> *fBuffer_jet_data_area;
std::vector<Float_t> *fBuffer_jet_data_leadingtrackpt;
std::vector<Float_t> *fBuffer_jet_data_perpconerho;
std::vector<UShort_t> *fBuffer_jet_data_nconstituents;


// cluster
std::vector<Float_t> *fBuffer_cluster_data_energy;
std::vector<Int_t> *fBuffer_cluster_data_definition;
std::vector<Float_t> *fBuffer_cluster_data_eta;
std::vector<Float_t> *fBuffer_cluster_data_phi;
std::vector<Float_t> *fBuffer_cluster_data_m02;
std::vector<Float_t> *fBuffer_cluster_data_m20;
std::vector<UShort_t> *fBuffer_cluster_data_ncells;
std::vector<Float_t> *fBuffer_cluster_data_time;
std::vector<Bool_t> *fBuffer_cluster_data_isexotic;
std::vector<Float_t> *fBuffer_cluster_data_distancebadchannel;
std::vector<UShort_t> *fBuffer_cluster_data_nlm;
std::vector<Float_t> *fBuffer_cluster_data_isoraw; 
std::vector<Float_t> *fBuffer_cluster_data_perpconerho;
std::vector<Float_t> *fBuffer_cluster_data_match_deta;
std::vector<Float_t> *fBuffer_cluster_data_match_dphi;
std::vector<Float_t> *fBuffer_cluster_data_match_p;

// define global switches
Bool_t fIsPbPb = kFALSE;
Bool_t fIsMC = kFALSE;

//----
// define output tree
//----
TTree *outputTree;

// TODO: MC truth information


#endif