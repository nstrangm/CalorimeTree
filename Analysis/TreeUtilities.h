#ifndef _TREEUTILITIES_
#define _TREEUTILITIES_

#include <fstream>
#include <TTree.h>
#include <TKey.h>
#include <TChain.h>

enum TreeTypes
{
  kRun2Tree,
  kRun3Tree,
  kBerkeleyTree
};

// Small helper function for readTree
TString GetTreeName(TString InputFileNamesFile)
{
  std::ifstream fileRootFileList;
  fileRootFileList.open(InputFileNamesFile, std::ios_base::in);
  if (!fileRootFileList)
    FATAL(Form("File listing input files %s not found!", InputFileNamesFile.Data()));

  for (TString tempLine; tempLine.ReadLine(fileRootFileList, kTRUE);)
  {
    if (tempLine.Contains("histos")) // Skip the file that contains the histograms
      continue;
    // before using check if file is health
    TFile *f = TFile::Open(tempLine.Data());
    if ((!f) || f->IsZombie())
    {
      ERROR(Form("Input file %s is not healthy! Skipping this one...", tempLine.Data()))
      continue;
    }
    TKey *key;
    TIter next(f->GetListOfKeys());
    key = (TKey *)next();
    return key->GetName();
  }

  FATAL("Not a single tree was found to find the name of the tree")
}

// Function that creates a chain, adds all trees listed in "InputFileNamesFile" to the chain and returns it
TChain *readTree(TString InputFileNamesFile, GlobalOptions optns)
{
  TString treeName = GetTreeName(InputFileNamesFile);
  TChain *tree = new TChain(treeName);

  std::ifstream fileRootFileList;
  fileRootFileList.open(InputFileNamesFile, std::ios_base::in);
  if (!fileRootFileList)
    FATAL(Form("File listing input files %s not found!", InputFileNamesFile.Data()));

  for (TString tempLine; tempLine.ReadLine(fileRootFileList, kTRUE);)
  {
    if (tempLine.Contains("histos")) // Skip the file that contains the histograms
      continue;
    // before adding check if file is health
    TFile *f = TFile::Open(tempLine.Data());
    if ((!f) || f->IsZombie())
    {
      ERROR(Form("Input file %s is not healthy! Skipping this one...", tempLine.Data()))
      continue;
    }
    f->Close();
    INFO(Form("Adding tree %s from %s to chain", treeName.Data(), tempLine.Data()))
    tree->Add(Form("%s", tempLine.Data()));
    if (optns.isDebugRun)
      break;
  }
  return tree;
}

int listTreeBranches(TChain *tree)
{
  TObjArray *branchList = tree->GetListOfBranches();
  TString sBranchList = "";
  bool ClustersInTree = false, JetsInTree = false, Pi0InTree = false;
  for (Int_t b = 0; b < branchList->GetEntries(); b++)
  {
    TString branchName = branchList->At(b)->GetName();
    if (branchName.Contains("Cluster") || branchName.Contains("Gamma") || branchName.Contains("Photon"))
      ClustersInTree = true;
    if (branchName.Contains("Jet") || branchName.Contains("Track"))
      JetsInTree = true;
    sBranchList += branchName;
    if (b < branchList->GetEntries() - 1)
      sBranchList += ", ";
  }
  INFO(Form("Your tree %s contains the following %d branches:%s\n", tree->GetName(), branchList->GetEntries(), sBranchList.Data()))
  if (sBranchList.Contains(", Event_NPrimaryTracks, Event_IsTriggered, Event_ZVertex, Event_Quality, Event_NotAccepted, "))
    return kRun2Tree;
  else if (sBranchList.Contains("event_multiplicity, event_centrality, event_rho, event_eventselection, event_alias")) // Add Run3 tree specific branches here
    return kRun3Tree;
  else if (sBranchList.Contains("Berkely")) // Add Berkely tree specific branches here
    return kBerkeleyTree;
  else
    FATAL("Tree not recognized based on its branches. Please add your tree format in the lines above!")
}

class TreeBuffer
{
private:
  void ReadRun2TreeIntoBuffer(TChain *tree, GlobalOptions optns);
  void ReadRun3TreeIntoBuffer(TChain *tree, GlobalOptions optns);
  void ReadBerkeleyTreeIntoBuffer(TChain *tree, GlobalOptions optns);

public:
  int NEvents = 0;

  // #############################################
  // Members that will be matched to DATA branches
  // #############################################

  // -------------- Event branches ---------------
  float Event_Rho = 0;                     // Run 2 + 3
  unsigned short Event_NPrimaryTracks = 0; // Run 2
  bool Event_IsTriggered = 0;              // Run 2
  float Event_ZVertex = 0;                 // Run 2
  unsigned short Event_Quality = 0;        // Run 2
  unsigned short Event_NotAccepted = 0;    // Run 2
  int Event_Multiplicity = 0;              // Run     3
  float Event_Centrality = 0;              // Run     3
  uint16_t Event_Selection = 0;            // Run     3
  uint32_t Event_Alias = 0;                // Run     3

  // -------------- Cluster branches ---------------
  std::vector<float> *Cluster_E = 0;                    // Run 2 + 3
  std::vector<float> *Cluster_Px = 0;                   // Run 2
  std::vector<float> *Cluster_Py = 0;                   // Run 2
  std::vector<float> *Cluster_Pz = 0;                   // Run 2
  std::vector<float> *Cluster_Eta = 0;                  // Run     3
  std::vector<float> *Cluster_Phi = 0;                  // Run     3
  std::vector<float> *Cluster_M02 = 0;                  // Run 2 + 3
  std::vector<float> *Cluster_M02Recalc = 0;            // Run 2
  std::vector<float> *Cluster_M20 = 0;                  // Run 2 + 3
  std::vector<unsigned short> *Cluster_NCells = 0;      // Run 2 + 3
  std::vector<float> *Cluster_V1SplitMass = 0;          // Run 2
  std::vector<float> *Cluster_Time = 0;                 // Run     3
  std::vector<float> *Cluster_MinMassDiffToPi0 = 0;     // Run 2
  std::vector<float> *Cluster_MinMassDiffToEta = 0;     // Run 2
  std::vector<unsigned short> *Cluster_NLM = 0;         // Run 2 + 3
  std::vector<unsigned short> *Cluster_SM = 0;          // Run 2
  std::vector<float> *Cluster_EFrac = 0;                // Run 2
  std::vector<float> *Cluster_IsoCharged1 = 0;          // Run 2     R = 0.2
  std::vector<float> *Cluster_IsoCharged2 = 0;          // Run 2     R = 0.3
  std::vector<float> *Cluster_IsoCharged3 = 0;          // Run 2 + 3 R = 0.4 (standard)
  std::vector<float> *Cluster_IsoBckPerp = 0;           // Run 2 + 3 In run2 raw, in run3 corrected
  std::vector<float> *Cluster_MatchTrackdEta = 0;       // Run 2
  std::vector<float> *Cluster_MatchTrackdPhi = 0;       // Run 2
  std::vector<float> *Cluster_MatchTrackP = 0;          // Run 2
  std::vector<float> *Cluster_MatchTrackPt = 0;         // Run 2
  std::vector<bool> *Cluster_MatchTrackIsConv = 0;      // Run 2
  std::vector<bool> *Cluster_IsExotic = 0;              // Run     3
  std::vector<float> *Cluster_DistanceToBadChannel = 0; // Run 2 Also available in run3 but as ushort and not used currently -> Omitted for now

  // -------------- Jet branches ---------------
  std::vector<float> *Jet_Pt = 0;                 // Run     3
  std::vector<float> *Jet_E = 0;                  // Run     3
  std::vector<float> *Jet_M = 0;                  // Run     3
  std::vector<float> *Jet_Px = 0;                 // Run 2
  std::vector<float> *Jet_Py = 0;                 // Run 2
  std::vector<float> *Jet_Pz = 0;                 // Run 2
  std::vector<float> *Jet_Radius = 0;                // Run     3
  std::vector<float> *Jet_Eta = 0;                // Run     3
  std::vector<float> *Jet_Phi = 0;                // Run     3
  std::vector<float> *Jet_Area = 0;               // Run 2 + 3
  std::vector<unsigned short> *Jet_Nch = 0;       // Run 2
  std::vector<unsigned short> *Jet_Nclus = 0;     // Run 2
  std::vector<float> *Jet_LeadingTrackPt = 0;     // Run     3
  std::vector<float> *Jet_PerpConeRho = 0;        // Run     3
  std::vector<unsigned short> *Jet_Nconstits = 0; // Run     3

  // ###########################################
  // Members that will be matched to MC branches
  // ###########################################

  // ------------ Event branches ---------------
  float Event_RhoMC;
  double Event_Weight;
  float Event_Xsection;
  unsigned short Event_Ntrials;

  // ------------ Cluster branches ---------------
  std::vector<float> *TrueCluster_E = 0;
  std::vector<float> *TrueCluster_Px = 0;
  std::vector<float> *TrueCluster_Py = 0;
  std::vector<float> *TrueCluster_Pz = 0;
  std::vector<float> *TrueCluster_LeadingEFrac = 0;
  std::vector<float> *TrueCluster_IsoCharged1 = 0;
  std::vector<float> *TrueCluster_IsoCharged2 = 0;
  std::vector<float> *TrueCluster_IsoCharged3 = 0;
  std::vector<float> *TrueCluster_IsoBckPerp = 0;
  std::vector<int> *TrueCluster_MCTag = 0;
  std::vector<bool> *TrueCluster_IsConv = 0;

  std::vector<float> *GenPhoton_E = 0;
  std::vector<float> *GenPhoton_Px = 0;
  std::vector<float> *GenPhoton_Py = 0;
  std::vector<float> *GenPhoton_Pz = 0;
  std::vector<float> *GenPhoton_IsoCharged1 = 0;
  std::vector<float> *GenPhoton_IsoCharged2 = 0;
  std::vector<float> *GenPhoton_IsoCharged3 = 0;
  std::vector<float> *GenPhoton_IsoBckPerp = 0;
  std::vector<int> *GenPhoton_MCTag = 0;
  std::vector<bool> *GenPhoton_IsConv = 0;

  // ------------ Jet branches ---------------
  std::vector<float> *PLJet_Px = 0; // PL = Particle Level
  std::vector<float> *PLJet_Py = 0;
  std::vector<float> *PLJet_Pz = 0;
  std::vector<float> *PLJet_Area = 0;
  std::vector<unsigned short> *PLJet_NPart = 0;

  TreeBuffer(TChain *tree, GlobalOptions optns);
  ~TreeBuffer() {};
};

TreeBuffer::TreeBuffer(TChain *tree, GlobalOptions optns)
{
  NEvents = tree->GetEntries();
  switch (optns.TreeFormat)
  {
  case kRun2Tree:
    ReadRun2TreeIntoBuffer(tree, optns);
    break;
  case kRun3Tree:
    ReadRun3TreeIntoBuffer(tree, optns);
    break;
  case kBerkeleyTree:
  default:
    FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
  }
}

void TreeBuffer::ReadRun2TreeIntoBuffer(TChain *tree, GlobalOptions optns)
{
  tree->SetBranchAddress("Event_Rho", &Event_Rho);
  tree->SetBranchAddress("Event_NPrimaryTracks", &Event_NPrimaryTracks);
  tree->SetBranchAddress("Event_IsTriggered", &Event_IsTriggered);
  tree->SetBranchAddress("Event_ZVertex", &Event_ZVertex);
  tree->SetBranchAddress("Event_Quality", &Event_Quality);
  tree->SetBranchAddress("Event_NotAccepted", &Event_NotAccepted);
  if (optns.isMC)
  {
    tree->SetBranchAddress("Event_RhoMC", &Event_RhoMC);
    tree->SetBranchAddress("Event_Weight", &Event_Weight);
    tree->SetBranchAddress("Event_Xsection", &Event_Xsection);
    tree->SetBranchAddress("Event_Ntrials", &Event_Ntrials);
  }
  if (optns.doIsoGamma)
  {
    tree->SetBranchAddress("Cluster_E", &Cluster_E);
    tree->SetBranchAddress("Cluster_Px", &Cluster_Px);
    tree->SetBranchAddress("Cluster_Py", &Cluster_Py);
    tree->SetBranchAddress("Cluster_Pz", &Cluster_Pz);
    tree->SetBranchAddress("Cluster_M02", &Cluster_M02);
    tree->SetBranchAddress("Cluster_M02Recalc", &Cluster_M02Recalc);
    tree->SetBranchAddress("Cluster_M20", &Cluster_M20);
    tree->SetBranchAddress("Cluster_NCells", &Cluster_NCells);
    tree->SetBranchAddress("Cluster_V1SplitMass", &Cluster_V1SplitMass);
    tree->SetBranchAddress("Cluster_MinMassDiffToPi0", &Cluster_MinMassDiffToPi0);
    tree->SetBranchAddress("Cluster_MinMassDiffToEta", &Cluster_MinMassDiffToEta);
    tree->SetBranchAddress("Cluster_NLM", &Cluster_NLM);
    tree->SetBranchAddress("Cluster_SM", &Cluster_SM);
    tree->SetBranchAddress("Cluster_EFrac", &Cluster_EFrac);
    tree->SetBranchAddress("Cluster_IsoCharged1", &Cluster_IsoCharged1);
    tree->SetBranchAddress("Cluster_IsoCharged2", &Cluster_IsoCharged2);
    tree->SetBranchAddress("Cluster_IsoCharged3", &Cluster_IsoCharged3);
    tree->SetBranchAddress("Cluster_IsoBckPerp", &Cluster_IsoBckPerp);
    tree->SetBranchAddress("Cluster_MatchTrackdEta", &Cluster_MatchTrackdEta);
    tree->SetBranchAddress("Cluster_MatchTrackdPhi", &Cluster_MatchTrackdPhi);
    tree->SetBranchAddress("Cluster_MatchTrackP", &Cluster_MatchTrackP);
    tree->SetBranchAddress("Cluster_MatchTrackPt", &Cluster_MatchTrackPt);
    tree->SetBranchAddress("Cluster_MatchTrackIsConv", &Cluster_MatchTrackIsConv);
    tree->SetBranchAddress("Cluster_DistanceToBadChannel", &Cluster_DistanceToBadChannel);
    if (optns.isMC)
    {
      tree->SetBranchAddress("TrueCluster_E", &TrueCluster_E);
      tree->SetBranchAddress("TrueCluster_Px", &TrueCluster_Px);
      tree->SetBranchAddress("TrueCluster_Py", &TrueCluster_Py);
      tree->SetBranchAddress("TrueCluster_Pz", &TrueCluster_Pz);
      tree->SetBranchAddress("TrueCluster_LeadingEFrac", &TrueCluster_LeadingEFrac);
      tree->SetBranchAddress("TrueCluster_MCIsoCharged1", &TrueCluster_IsoCharged1);
      tree->SetBranchAddress("TrueCluster_MCIsoCharged2", &TrueCluster_IsoCharged2);
      tree->SetBranchAddress("TrueCluster_MCIsoCharged3", &TrueCluster_IsoCharged3);
      tree->SetBranchAddress("TrueCluster_MCIsoBckPerp", &TrueCluster_IsoBckPerp);
      tree->SetBranchAddress("TrueCluster_MCTagMore", &TrueCluster_MCTag);
      tree->SetBranchAddress("TrueCluster_IsConv", &TrueCluster_IsConv);
      tree->SetBranchAddress("GenPhoton_E", &GenPhoton_E);
      tree->SetBranchAddress("GenPhoton_Px", &GenPhoton_Px);
      tree->SetBranchAddress("GenPhoton_Py", &GenPhoton_Py);
      tree->SetBranchAddress("GenPhoton_Pz", &GenPhoton_Pz);
      tree->SetBranchAddress("GenPhoton_MCIsoCharged1", &GenPhoton_IsoCharged1);
      tree->SetBranchAddress("GenPhoton_MCIsoCharged2", &GenPhoton_IsoCharged2);
      tree->SetBranchAddress("GenPhoton_MCIsoCharged3", &GenPhoton_IsoCharged3);
      tree->SetBranchAddress("GenPhoton_MCIsoBckPerp", &GenPhoton_IsoBckPerp);
      tree->SetBranchAddress("GenPhoton_MCTag", &GenPhoton_MCTag);
      tree->SetBranchAddress("GenPhoton_IsConv", &GenPhoton_IsConv);
    }
  }
  if (optns.doJets)
  {
    tree->SetBranchAddress("Jet_Px", &Jet_Px);
    tree->SetBranchAddress("Jet_Py", &Jet_Py);
    tree->SetBranchAddress("Jet_Pz", &Jet_Pz);
    tree->SetBranchAddress("Jet_Area", &Jet_Area);
    tree->SetBranchAddress("Jet_Nch", &Jet_Nch);
    tree->SetBranchAddress("Jet_Nclus", &Jet_Nclus);
    if (optns.isMC)
    {
      tree->SetBranchAddress("TrueJet_Px", &PLJet_Px);
      tree->SetBranchAddress("TrueJet_Py", &PLJet_Py);
      tree->SetBranchAddress("TrueJet_Pz", &PLJet_Pz);
      tree->SetBranchAddress("TrueJet_Area", &PLJet_Area);
      tree->SetBranchAddress("TrueJet_NPart", &PLJet_NPart);
    }
  }
}

void TreeBuffer::ReadRun3TreeIntoBuffer(TChain *tree, GlobalOptions optns)
{
  if (optns.isMC)
  {
    FATAL("The branches of run3 MC trees are not yet implemented. Please do so here!") // TODO: MC run3 tree branches
  }
  tree->SetBranchAddress("event_multiplicity", &Event_Multiplicity);
  tree->SetBranchAddress("event_centrality", &Event_Centrality);
  tree->SetBranchAddress("event_rho", &Event_Rho);
  tree->SetBranchAddress("event_eventselection", &Event_Selection);
  tree->SetBranchAddress("event_alias", &Event_Alias);
  if (optns.doIsoGamma)
  {
    tree->SetBranchAddress("cluster_data_energy", &Cluster_E);
    tree->SetBranchAddress("cluster_data_eta", &Cluster_Eta);
    tree->SetBranchAddress("cluster_data_phi", &Cluster_Phi);
    tree->SetBranchAddress("cluster_data_m02", &Cluster_M02);
    tree->SetBranchAddress("cluster_data_m20", &Cluster_M20);
    tree->SetBranchAddress("cluster_data_ncells", &Cluster_NCells);
    tree->SetBranchAddress("cluster_data_time", &Cluster_Time);
    tree->SetBranchAddress("cluster_data_nlm", &Cluster_NLM);
    tree->SetBranchAddress("cluster_data_isoraw", &Cluster_IsoCharged3);
    tree->SetBranchAddress("cluster_data_isexotic", &Cluster_IsExotic);
    tree->SetBranchAddress("cluster_data_perpconerho", &Cluster_IsoBckPerp);
    tree->SetBranchAddress("cluster_data_match_deta", &Cluster_MatchTrackdEta);
    tree->SetBranchAddress("cluster_data_match_dphi", &Cluster_MatchTrackdPhi);
    tree->SetBranchAddress("cluster_data_match_p", &Cluster_MatchTrackP);
  }
  if (optns.doJets)
  {
    tree->SetBranchAddress("jet_data_pt", &Jet_Pt);
    tree->SetBranchAddress("jet_data_eta", &Jet_Eta);
    tree->SetBranchAddress("jet_data_phi", &Jet_Phi);
    tree->SetBranchAddress("jet_data_radius", &Jet_Radius);
    tree->SetBranchAddress("jet_data_energy", &Jet_E);
    tree->SetBranchAddress("jet_data_mass", &Jet_M);
    tree->SetBranchAddress("jet_data_area", &Jet_Area);
    tree->SetBranchAddress("jet_data_leadingtrackpt", &Jet_LeadingTrackPt);
    tree->SetBranchAddress("jet_data_perpconerho", &Jet_PerpConeRho);
    tree->SetBranchAddress("jet_data_nconstituents", &Jet_Nconstits);
  }
}

#endif