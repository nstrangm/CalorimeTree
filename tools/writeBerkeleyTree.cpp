#include "writeBerkeleyTree.h"
#include "eventbuilding.h"

// flat trees
// event prpoperties tree
// track tree
// cluster tree

// they all have event index and then group by operations to get events and
// iterate accordingly

void createQAHistos() {
  hTrackPt = new TH1F("hTrackPt", "Track p_{T}", 100, 0, 100);
  hClusterEnergy = new TH1F("hClusterEnergy", "Cluster Energy", 100, 0, 100);
  hClusterEta = new TH1F("hClusterEta", "Cluster #eta", 100, -1, 1);
  hClusterPhi =
      new TH1F("hClusterPhi", "Cluster #phi", 100, 0, 2 * TMath::Pi());
  hNEvents = new TH1F("hNEvents", "Number of events", 2, 0, 2);
  hEvtVtxX = new TH1F("hEvtVtxX", "Event vertex x", 100, -10, 10);
  hEvtVtxY = new TH1F("hEvtVtxY", "Event vertex y", 100, -10, 10);
  hEvtVtxZ = new TH1F("hEvtVtxZ", "Event vertex z", 100, -10, 10);

  hClusterM02vsE =
      new TH2F("hClusterM02vsE", "Cluster M02 vs E", 300, 0, 3, 100, 0, 100);

  outputhists = new TList();
  // add them all to histos list
  outputhists->Add(hTrackPt);
  outputhists->Add(hClusterEnergy);
  outputhists->Add(hClusterEta);
  outputhists->Add(hClusterPhi);
  outputhists->Add(hNEvents);
  outputhists->Add(hEvtVtxX);
  outputhists->Add(hEvtVtxY);
  outputhists->Add(hEvtVtxZ);
  outputhists->Add(hClusterM02vsE);
}
void createTree() {
  outputTree = new TTree("eventTree", "eventTree");
  outputTree->Branch("run_number", &fBuffer_RunNumber, "RunNumber/I");
  outputTree->Branch("event_selection", &fBuffer_eventselection,
                     "eventselection/s"); // todo check format
  outputTree->Branch("triggersel", &fBuffer_triggersel, "triggersel/l");
  outputTree->Branch("centrality", &fBuffer_centrality, "centrality/F");
  outputTree->Branch("multiplicity", &fBuffer_multiplicity, "multiplicity/F");
  // track
  outputTree->Branch("track_data_eta", &fBuffer_track_data_eta);
  outputTree->Branch("track_data_phi", &fBuffer_track_data_phi);
  outputTree->Branch("track_data_pt", &fBuffer_track_data_pt);
  outputTree->Branch("track_data_label", &fBuffer_track_data_label);
  outputTree->Branch("track_data_tracksel", &fBuffer_track_data_tracksel);
  // cluster
  outputTree->Branch("cluster_data_energy", &fBuffer_cluster_data_energy);
  outputTree->Branch("cluster_data_eta", &fBuffer_cluster_data_eta);
  outputTree->Branch("cluster_data_phi", &fBuffer_cluster_data_phi);
  outputTree->Branch("cluster_data_m02", &fBuffer_cluster_data_m02);
  outputTree->Branch("cluster_data_m20", &fBuffer_cluster_data_m20);
  outputTree->Branch("cluster_data_ncells", &fBuffer_cluster_data_ncells);
  outputTree->Branch("cluster_data_time", &fBuffer_cluster_data_time);
  outputTree->Branch("cluster_data_isexotic", &fBuffer_cluster_data_isexotic);
  outputTree->Branch("cluster_data_distancebadchannel",
                     &fBuffer_cluster_data_distancebadchannel);
  outputTree->Branch("cluster_data_nlm", &fBuffer_cluster_data_nlm);
  outputTree->Branch("cluster_data_clusterdef",
                     &fBuffer_cluster_data_clusterdef);
  outputTree->Branch("cluster_data_matchedTrackIndex",
                     &fBuffer_cluster_data_matchedTrackIndex);
}

void clearBuffers() {
  fBuffer_track_data_eta->clear();
  fBuffer_track_data_phi->clear();
  fBuffer_track_data_pt->clear();
  fBuffer_track_data_label->clear();
  fBuffer_track_data_tracksel->clear();

  fBuffer_cluster_data_energy->clear();
  fBuffer_cluster_data_eta->clear();
  fBuffer_cluster_data_phi->clear();
  fBuffer_cluster_data_m02->clear();
  fBuffer_cluster_data_m20->clear();
  fBuffer_cluster_data_ncells->clear();
  fBuffer_cluster_data_time->clear();
  fBuffer_cluster_data_isexotic->clear();
  fBuffer_cluster_data_distancebadchannel->clear();
  fBuffer_cluster_data_nlm->clear();
  fBuffer_cluster_data_clusterdef->clear();
  fBuffer_cluster_data_matchedTrackIndex->clear();
}

// write events to TTree
void writeEvents(TTree *tree, std::vector<event> &events) {
  for (auto &ev : events) {
    // clear all buffers
    clearBuffers();

    if (treecuts["event_cuts"]["min_clus_E"].as<float>() > 0.) {
      bool acc = false;
      for (auto &cl : ev.clusters) {
        if (cl.energy > treecuts["event_cuts"]["min_clus_E"].as<float>()) {
          acc = true;
          break;
        }
      }
      if (!acc)
        continue;
    }

    // fill event level properties
    fBuffer_RunNumber = (Int_t)ev.col.runNumber;
    fBuffer_eventselection = (uint16_t)ev.col.eventsel;
    fBuffer_centrality = (Float_t)ev.col.centrality;
    fBuffer_multiplicity = (Float_t)ev.col.multiplicity;
    fBuffer_triggersel = (uint64_t)ev.col.triggersel;

    trackCuts = treecuts["track_cuts"];
    // fill track properties
    for (auto &tr : ev.tracks) {
      if (tr.pt < trackCuts["min_pt"].as<float>())
        continue;
      if (tr.eta < trackCuts["min_eta"].as<float>())
        continue;
      if (tr.eta > trackCuts["max_eta"].as<float>())
        continue;

      fBuffer_track_data_eta->push_back((Float_t)tr.eta);
      fBuffer_track_data_phi->push_back((Float_t)tr.phi);
      fBuffer_track_data_pt->push_back((Float_t)tr.pt);
      fBuffer_track_data_label->push_back(0.); // TODO placeholder
      fBuffer_track_data_tracksel->push_back((uint8_t)tr.trackSel);
    }

    // fill cluster properties
    for (auto &cl : ev.clusters) {
      fBuffer_cluster_data_energy->push_back((Float_t)cl.energy);
      fBuffer_cluster_data_eta->push_back((Float_t)cl.eta);
      fBuffer_cluster_data_phi->push_back((Float_t)cl.phi);
      fBuffer_cluster_data_m02->push_back((Float_t)cl.m02);
      fBuffer_cluster_data_m20->push_back((Float_t)cl.m20);
      fBuffer_cluster_data_ncells->push_back((UShort_t)cl.ncells);
      fBuffer_cluster_data_time->push_back((Float_t)cl.time);
      fBuffer_cluster_data_isexotic->push_back((Bool_t)cl.isexotic);
      fBuffer_cluster_data_distancebadchannel->push_back(
          (UShort_t)cl.distancebadchannel);
      fBuffer_cluster_data_nlm->push_back((UShort_t)cl.nlm);
      fBuffer_cluster_data_clusterdef->push_back((UShort_t)cl.clusterdef);
      fBuffer_cluster_data_matchedTrackIndex->push_back(
          (UShort_t)cl.matchedTrackIndex);
    }

    // fill tree
    tree->Fill();
  }
}

void doEventSelection(std::vector<event> &events) {
  // loop over events and remove them from vector if they don't fulfull a
  // certain cut possibility to do event selection if wanted
  eventCuts = treecuts["event_cuts"];

  // loop over event vector and erase if not fulfilling cuts
  std::vector<event>::iterator it = events.begin();
  while (it != events.end()) {
    if (TMath::Abs(it->col.posz) > eventCuts["z_vtx_cut"].as<float>()) {
      it = events.erase(it);
    } else if ((it->col.multiplicity < eventCuts["mult_min"].as<float>()) ||
               (it->col.multiplicity > eventCuts["mult_max"].as<float>())) {
      it = events.erase(it);

    } else {
      ++it;
    }
  }
  return;
}

// can be used to do analysis (if needed)
void doAnalysis(std::vector<event> &events) {
  for (auto &ev : events) {
    hNEvents->Fill(1);
    hEvtVtxX->Fill(ev.col.posx);
    hEvtVtxY->Fill(ev.col.posy);
    hEvtVtxZ->Fill(ev.col.posz);

    // plot pt of all tracks
    for (auto &tr : ev.tracks) {
      hTrackPt->Fill(tr.pt);
    }

    // plot energy, eta and phi of all clusters
    for (auto &cl : ev.clusters) {
      hClusterEnergy->Fill(cl.energy);
      hClusterEta->Fill(cl.eta);
      hClusterPhi->Fill(cl.phi);

      hClusterM02vsE->Fill(cl.m02, cl.energy);
    }
  }
}
void processFile(TFile *inputFile) {
  // define all TTrees that are in the file
  TTree *O2jclustertrack, *O2jcollision, *O2jtrack, *O2jcluster, *O2jbc;
  // set all to nullptr
  O2jclustertrack = nullptr;
  O2jcollision = nullptr;
  O2jtrack = nullptr;
  O2jcluster = nullptr;
  O2jbc = nullptr;

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
      if (strcmp(tree->GetName(), "O2jclustertrack") == 0) {
        O2jclustertrack = tree;
      } else if (strcmp(tree->GetName(), "O2jcollision") == 0) {
        O2jcollision = tree;
      } else if (strcmp(tree->GetName(), "O2jtrack") == 0) {
        O2jtrack = tree;
      } else if (strcmp(tree->GetName(), "O2jcluster") == 0) {
        O2jcluster = tree;
      } else if (strcmp(tree->GetName(), "O2jbc") == 0) {
        O2jbc = tree;
      }
    }

    // build event
    std::vector<event> events =
        buildEvents(O2jcollision, O2jbc, O2jtrack, O2jcluster, O2jclustertrack);

    // do event selection
    doEventSelection(events);

    // loop over events
    doAnalysis(events);

    // write events to TTree
    writeEvents(outputTree, events);

    // delete all TTrees
    delete O2jclustertrack;
    delete O2jcollision;
    delete O2jtrack;
    delete O2jcluster;

    // delete all events
    events.clear();
  }
}
void writeBerkeleyTree(TString inputFileList = "",
                       TString outputFileName = "output/test.root",
                       TString configFile = "treeCutsBerkeley.yaml",
                       Bool_t isPbPb = kFALSE, Bool_t isMC = kFALSE) {

  TFile *outFile = new TFile(outputFileName.Data(), "RECREATE");
  fIsPbPb = isPbPb;
  fIsMC = isMC;

  treecuts = YAML::LoadFile(configFile.Data());

  createQAHistos();
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
  outputhists->Write();
  outputTree->Write();

  return;
}