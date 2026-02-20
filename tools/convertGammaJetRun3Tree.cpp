#include "convertGammaJetRun3Tree.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>

using namespace std;

// flat trees
// event prpoperties tree
// track tree
// cluster tree

// they all have event index and then group by operations to get events and
// iterate accordingly

void createTree() {
  outputTree = new TTree("eventTree", "eventTree");
  outputTree->Branch("event_multiplicity", &fBuffer_multiplicity,
                     "event_multiplicity/F");
    outputTree->Branch("event_centrality", &fBuffer_centrality,
                     "event_centrality/F");
    outputTree->Branch("event_rho", &fBuffer_rho,
                        "event_rho/F");
    outputTree->Branch("event_eventselection", &fBuffer_eventselection,
                        "event_eventselection/s");
    outputTree->Branch("event_occupancy", &fBuffer_occupancy,
                        "event_occupancy/I");
    outputTree->Branch("event_alias", &fBuffer_alias);
    outputTree->Branch("jet_data_pt", &fBuffer_jet_data_pt);
    outputTree->Branch("jet_data_eta", &fBuffer_jet_data_eta);
    outputTree->Branch("jet_data_phi", &fBuffer_jet_data_phi);
    outputTree->Branch("jet_data_radius", &fBuffer_jet_data_radius);
    outputTree->Branch("jet_data_energy", &fBuffer_jet_data_energy);
    outputTree->Branch("jet_data_mass", &fBuffer_jet_data_mass);
    outputTree->Branch("jet_data_area", &fBuffer_jet_data_area);
    outputTree->Branch("jet_data_leadingtrackpt", &fBuffer_jet_data_leadingtrackpt);
    outputTree->Branch("jet_data_perpconerho", &fBuffer_jet_data_perpconerho);
    outputTree->Branch("jet_data_nconstituents", &fBuffer_jet_data_nconstituents);
    outputTree->Branch("cluster_data_energy", &fBuffer_cluster_data_energy);
    outputTree->Branch("cluster_data_definition", &fBuffer_cluster_data_definition);
    outputTree->Branch("cluster_data_eta", &fBuffer_cluster_data_eta);
    outputTree->Branch("cluster_data_phi", &fBuffer_cluster_data_phi);
    outputTree->Branch("cluster_data_m02", &fBuffer_cluster_data_m02);
    outputTree->Branch("cluster_data_m20", &fBuffer_cluster_data_m20);
    outputTree->Branch("cluster_data_ncells", &fBuffer_cluster_data_ncells);
    outputTree->Branch("cluster_data_time", &fBuffer_cluster_data_time);
    outputTree->Branch("cluster_data_isexotic", &fBuffer_cluster_data_isexotic);
    outputTree->Branch("cluster_data_distancebadchannel", &fBuffer_cluster_data_distancebadchannel);
    outputTree->Branch("cluster_data_nlm", &fBuffer_cluster_data_nlm);
    outputTree->Branch("cluster_data_isoraw", &fBuffer_cluster_data_isoraw);
    outputTree->Branch("cluster_data_perpconerho", &fBuffer_cluster_data_perpconerho);
    outputTree->Branch("cluster_data_match_deta", &fBuffer_cluster_data_match_deta);
    outputTree->Branch("cluster_data_match_dphi", &fBuffer_cluster_data_match_dphi);
    outputTree->Branch("cluster_data_match_p", &fBuffer_cluster_data_match_p);
    outputTree->Branch("jetsubstr_data_nentries",    &fBuffer_jetsubstr_data_nentries);
    outputTree->Branch("jetsubstr_data_energymother", &fBuffer_jetsubstr_data_energymother);
    outputTree->Branch("jetsubstr_data_ptleading",    &fBuffer_jetsubstr_data_ptleading);
    outputTree->Branch("jetsubstr_data_ptsubleading", &fBuffer_jetsubstr_data_ptsubleading);
    outputTree->Branch("jetsubstr_data_theta",        &fBuffer_jetsubstr_data_theta);

}

void clearBuffers(){
    fBuffer_jet_data_pt->clear();
    fBuffer_jet_data_eta->clear();
    fBuffer_jet_data_phi->clear();
    fBuffer_jet_data_radius->clear();
    fBuffer_jet_data_energy->clear();
    fBuffer_jet_data_mass->clear();
    fBuffer_jet_data_area->clear();
    fBuffer_jet_data_leadingtrackpt->clear();
    fBuffer_jet_data_perpconerho->clear();
    fBuffer_jet_data_nconstituents->clear();
    fBuffer_cluster_data_energy->clear();
    fBuffer_cluster_data_definition->clear();
    fBuffer_cluster_data_eta->clear();
    fBuffer_cluster_data_phi->clear();
    fBuffer_cluster_data_m02->clear();
    fBuffer_cluster_data_m20->clear();
    fBuffer_cluster_data_ncells->clear();
    fBuffer_cluster_data_time->clear();
    fBuffer_cluster_data_isexotic->clear();
    fBuffer_cluster_data_distancebadchannel->clear();
    fBuffer_cluster_data_nlm->clear();
    fBuffer_cluster_data_isoraw->clear();
    fBuffer_cluster_data_perpconerho->clear();
    fBuffer_cluster_data_match_deta->clear();
    fBuffer_cluster_data_match_dphi->clear();
    fBuffer_cluster_data_match_p->clear();
    fBuffer_jetsubstr_data_nentries->clear();
    fBuffer_jetsubstr_data_energymother->clear();
    fBuffer_jetsubstr_data_ptleading->clear();
    fBuffer_jetsubstr_data_ptsubleading->clear();
    fBuffer_jetsubstr_data_theta->clear();
}

void readAndFillTree(TTree* tEvents, TTree* tClusters, TTree* tJets, TTree* tJetSubstr){
    std::unordered_map<int, std::vector<int>> jetMap;
    std::unordered_map<int, std::vector<int>> clusterMap;
    std::unordered_map<int, std::vector<int>> jetSubstrMap;
    if (!tEvents) {
    cout << "-> Collisions is not found. Skipping this file" << endl;
        return;
    }
    if (!tClusters) {
    cout << "-> Clusters is not found. Skipping this file" << endl;
        return;
    }
    if (!tJets) {
    cout << "-> Jets is not found. Skipping this file" << endl;
        return;
    }
    if (!tJetSubstr) {
    cout << "-> Jet substructure is not found. I will not convert this to the output tree" << endl;
        // return;
    }

    // loop over all jets
     for (int j = 0; j < tJets->GetEntries(); j++) {
      tJets->GetEntry(j);
      int collisionID = tJets->GetBranch("fIndexGjEvents")->GetLeaf("fIndexGjEvents")->GetValue();
      jetMap[collisionID].push_back(j);
    }

    // loop over all clusters
    for (int j = 0; j < tClusters->GetEntries(); j++) {
      tClusters->GetEntry(j);
      int collisionID = tClusters->GetBranch("fIndexGjEvents")->GetLeaf("fIndexGjEvents")->GetValue();
      clusterMap[collisionID].push_back(j);
    }

    // Set up branch addresses for jet substructure (vector<vector<float>> per event)
    const Int_t kMaxSplittings = 100;

    Int_t   substr_energyMother_size = 0;
    Float_t substr_energyMother[kMaxSplittings];
    Int_t   substr_ptLeading_size = 0;
    Float_t substr_ptLeading[kMaxSplittings];
    Int_t   substr_ptSubLeading_size = 0;
    Float_t substr_ptSubLeading[kMaxSplittings];
    Int_t   substr_theta_size = 0;
    Float_t substr_theta[kMaxSplittings];

    if (tJetSubstr) {
      tJetSubstr->SetBranchAddress("fEnergyMother_size", &substr_energyMother_size);
      tJetSubstr->SetBranchAddress("fEnergyMother", substr_energyMother);
      tJetSubstr->SetBranchAddress("fPtLeading_size", &substr_ptLeading_size);
      tJetSubstr->SetBranchAddress("fPtLeading", substr_ptLeading);
      tJetSubstr->SetBranchAddress("fPtSubLeading_size", &substr_ptSubLeading_size);
      tJetSubstr->SetBranchAddress("fPtSubLeading", substr_ptSubLeading);
      tJetSubstr->SetBranchAddress("fTheta_size", &substr_theta_size);
      tJetSubstr->SetBranchAddress("fTheta", substr_theta);

      // Map collision ID to jet substructure entry (one entry per event)
      for (int j = 0; j < tJetSubstr->GetEntries(); j++) {
        tJetSubstr->GetEntry(j);
        int collisionID = tJetSubstr->GetBranch("fIndexGjEvents")->GetLeaf("fIndexGjEvents")->GetValue();
        jetSubstrMap[collisionID].push_back(j);
      }
    }

    // loop over all events
    for (int j = 0; j < tEvents->GetEntries(); j++) {
      tEvents->GetEntry(j);
      int collisionID = j;
      
      fBuffer_multiplicity = tEvents->GetBranch("fMultiplicity")->GetLeaf("fMultiplicity")->GetValue();
      fBuffer_centrality = tEvents->GetBranch("fCentrality")->GetLeaf("fCentrality")->GetValue();
      fBuffer_rho = tEvents->GetBranch("fRho")->GetLeaf("fRho")->GetValue();
      fBuffer_eventselection = tEvents->GetBranch("fEventSel")->GetLeaf("fEventSel")->GetValue();
      fBuffer_alias = tEvents->GetBranch("fAlias")->GetLeaf("fAlias")->GetValue();
      fBuffer_occupancy = tEvents->GetBranch("fOccupancy")->GetLeaf("fOccupancy")->GetValue();
      // 

      // loop over jets for this event
        for (int k = 0; k < jetMap[collisionID].size(); k++) {
            tJets->GetEntry(jetMap[collisionID].at(k));
            fBuffer_jet_data_pt->push_back((Float_t) tJets->GetBranch("fPt")->GetLeaf("fPt")->GetValue());
            fBuffer_jet_data_eta->push_back((Float_t) tJets->GetBranch("fEta")->GetLeaf("fEta")->GetValue());
            fBuffer_jet_data_phi->push_back((Float_t)tJets->GetBranch("fPhi")->GetLeaf("fPhi")->GetValue());
            fBuffer_jet_data_radius->push_back((Float_t)tJets->GetBranch("fRadius")->GetLeaf("fRadius")->GetValue() / 100.);
            fBuffer_jet_data_energy->push_back((Float_t)tJets->GetBranch("fEnergy")->GetLeaf("fEnergy")->GetValue());
            fBuffer_jet_data_mass->push_back((Float_t)tJets->GetBranch("fMass")->GetLeaf("fMass")->GetValue());
            fBuffer_jet_data_area->push_back((Float_t)tJets->GetBranch("fArea")->GetLeaf("fArea")->GetValue());
            fBuffer_jet_data_leadingtrackpt->push_back((Float_t)tJets->GetBranch("fLeadingTrackPt")->GetLeaf("fLeadingTrackPt")->GetValue());
            fBuffer_jet_data_perpconerho->push_back((Float_t)tJets->GetBranch("fPerpConeRho")->GetLeaf("fPerpConeRho")->GetValue());
            fBuffer_jet_data_nconstituents->push_back((UShort_t)tJets->GetBranch("fNConstituents")->GetLeaf("fNConstituents")->GetValue());
        }


        // loop over clusters for this event
        for (int k = 0; k < clusterMap[collisionID].size(); k++) {
            tClusters->GetEntry(clusterMap[collisionID].at(k));
            fBuffer_cluster_data_energy->push_back(tClusters->GetBranch("fEnergy")->GetLeaf("fEnergy")->GetValue());
            fBuffer_cluster_data_definition->push_back(tClusters->GetBranch("fDefinition")->GetLeaf("fDefinition")->GetValue());
            fBuffer_cluster_data_eta->push_back(tClusters->GetBranch("fEta")->GetLeaf("fEta")->GetValue());
            fBuffer_cluster_data_phi->push_back(tClusters->GetBranch("fPhi")->GetLeaf("fPhi")->GetValue());
            fBuffer_cluster_data_m02->push_back(tClusters->GetBranch("fM02")->GetLeaf("fM02")->GetValue());
            fBuffer_cluster_data_m20->push_back(tClusters->GetBranch("fM20")->GetLeaf("fM20")->GetValue());
            fBuffer_cluster_data_ncells->push_back(tClusters->GetBranch("fNCells")->GetLeaf("fNCells")->GetValue());
            fBuffer_cluster_data_time->push_back(tClusters->GetBranch("fTime")->GetLeaf("fTime")->GetValue());
            fBuffer_cluster_data_isexotic->push_back(tClusters->GetBranch("fIsExotic")->GetLeaf("fIsExotic")->GetValue());
            fBuffer_cluster_data_distancebadchannel->push_back(tClusters->GetBranch("fDistanceToBadChannel")->GetLeaf("fDistanceToBadChannel")->GetValue());
            fBuffer_cluster_data_nlm->push_back(tClusters->GetBranch("fNLM")->GetLeaf("fNLM")->GetValue());
            fBuffer_cluster_data_isoraw->push_back(tClusters->GetBranch("fIsoRaw")->GetLeaf("fIsoRaw")->GetValue());
            fBuffer_cluster_data_perpconerho->push_back(tClusters->GetBranch("fPerpConeRho")->GetLeaf("fPerpConeRho")->GetValue());
            fBuffer_cluster_data_match_deta->push_back(tClusters->GetBranch("fTMdeltaEta")->GetLeaf("fTMdeltaEta")->GetValue());
            fBuffer_cluster_data_match_dphi->push_back(tClusters->GetBranch("fTMdeltaPhi")->GetLeaf("fTMdeltaPhi")->GetValue());
            fBuffer_cluster_data_match_p->push_back(tClusters->GetBranch("fTMtrackP")->GetLeaf("fTMtrackP")->GetValue());
        }

        // Read jet substructure for this event (if exists).
        // One substructure entry per jet (same ordering as jetMap).
        // We use the flat array + count pattern:
        //   nentries[k]  = number of splittings for jet k
        //   flat arrays  = all splittings for all jets concatenated in jet order
        // To recover splittings for jet k: offset = sum(nentries[0..k-1])
        if (tJetSubstr) {
            for (int k = 0; k < (int)jetSubstrMap[collisionID].size(); k++) {
                tJetSubstr->GetEntry(jetSubstrMap[collisionID].at(k));

                fBuffer_jetsubstr_data_nentries->push_back(substr_energyMother_size);
                for (int s = 0; s < substr_energyMother_size; s++) {
                    fBuffer_jetsubstr_data_energymother->push_back(substr_energyMother[s]);
                    fBuffer_jetsubstr_data_ptleading->push_back(substr_ptLeading[s]);
                    fBuffer_jetsubstr_data_ptsubleading->push_back(substr_ptSubLeading[s]);
                    fBuffer_jetsubstr_data_theta->push_back(substr_theta[s]);
                }
            }
        }

        outputTree->Fill();
        clearBuffers();

    }
    
}


void processFile(TFile *inputFile) {
  // set all to nullptr
  TTree * O2gjevent = nullptr;
  TTree * O2gjgamma = nullptr;
  TTree * O2gjjet = nullptr;
  TTree * O2gjjetsubstr = nullptr;

  // loop over all directories and print name
  TIter next(inputFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TDirectory"))
      continue;
    TDirectory *dir = (TDirectory *)key->ReadObj();
    TIter next2(dir->GetListOfKeys());
    TKey *key2;
    while ((key2 = (TKey *)next2())) {
      TClass *cl2 = gROOT->GetClass(key2->GetClassName());
      if (!cl2->InheritsFrom("TTree"))
        continue;
      TTree *tree = (TTree *)key2->ReadObj();
      if (strcmp(tree->GetName(), "O2gjevent") == 0) {
        O2gjevent = tree;
      } else if (strcmp(tree->GetName(), "O2gjgamma") == 0) {
        O2gjgamma = tree;
      } else if (strcmp(tree->GetName(), "O2gjchjet") == 0) {
        O2gjjet = tree;
      } else if (strcmp(tree->GetName(), "O2gjjetsubstr") == 0) {
        O2gjjetsubstr = tree;
      }
    }
    cout << "Processing directory " << dir->GetName() << endl;
    // build event
    readAndFillTree(O2gjevent, O2gjgamma, O2gjjet, O2gjjetsubstr);

    // delete all TTrees
    delete O2gjevent;
    delete O2gjgamma;
    delete O2gjjet;
    delete O2gjjetsubstr;
  }
}
void convertGammaJetRun3Tree(TString inputFileList = "",
                       TString outputFileName = "output/test.root",
                      Bool_t isMC = kFALSE) {

  TFile *outFile = new TFile(outputFileName.Data(), "RECREATE");
  fIsMC = isMC;

  createTree();

  // loop over all files in txt file fileList
  std::vector<TString> fileList;
  std::ifstream file(inputFileList.Data());
  std::string str;
  while (std::getline(file, str)) {
    fileList.push_back(str);
  }

  for (Int_t i = 0; i < fileList.size(); i++) {
    TString filePath = fileList.at(i);
    cout << "-> Processing file " << filePath << endl;
    TFile *in = new TFile(filePath.Data());
    processFile(in);
    in->Close();
  }

  outFile->cd();
  // LOOP OVER OUTput hists and write to file
  outputTree->Write();

  return;
}