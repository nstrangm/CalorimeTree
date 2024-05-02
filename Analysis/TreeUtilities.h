#ifndef _TREEUTILITIES_
#define _TREEUTILITIES_

#include "Utilities.h"
#include <fstream>
#include "TTree.h"
#include "TChain.h"

enum TreeTypes { kRun2Tree,
                 kRun3Tree,
                 kBerkeleyTree };

// Function that creates a chain, adds all trees listed in "InputFileNamesFile" to the chain and returns it
TChain* readTree(TString InputFileNamesFile)
{
  TString treeName = "CaloTree";
  TChain* tree = new TChain(treeName);

  std::ifstream fileRootFileList;
  fileRootFileList.open(InputFileNamesFile, std::ios_base::in);
  if (!fileRootFileList)
    FATAL(Form("File listing input files %s not found!", InputFileNamesFile.Data()));

  for (TString tempLine; tempLine.ReadLine(fileRootFileList, kTRUE);) {
    if (tempLine.Contains("histos")) // Skip the file that contains the histograms
      continue;
    // before adding check if file is health
    TFile* f = TFile::Open(tempLine.Data());
    if ((!f) || f->IsZombie()) {
      ERROR(Form("Input file %s is not healthy! Skipping this one...", tempLine.Data()))
      continue;
    }
    f->Close();
    INFO(Form("Adding tree %s from %s to chain", treeName.Data(), tempLine.Data()))
    tree->Add(Form("%s", tempLine.Data()));
  }
  return tree;
}

int listTreeBranches(TChain* tree)
{
  INFO("------------------------------------------------------")
  INFO(Form("Your tree %s contains the following branches:", tree->GetName()))
  TObjArray* branchList = tree->GetListOfBranches();
  for (Int_t b = 0; b < branchList->GetEntries(); b++) {
    TString branchName = branchList->At(b)->GetName();
    INFO(Form("Branch %02d name: %s", b, branchName.Data())) // TODO: Check that branch names of requested analyses are there
  }
  INFO("-----------------------------------------------------")
  return kRun2Tree; // TODO: Check what tree type this is based on specific branch names
}

class TreeBuffer
{
 private:
  void ReadRun2TreeIntoBuffer(TChain* tree, GlobalOptions optns);
  void ReadRun3TreeIntoBuffer(TChain* tree, GlobalOptions optns);
  void ReadBerkeleyTreeIntoBuffer(TChain* tree, GlobalOptions optns);

 public:
  int NEvents;
  std::vector<Float_t>* Cluster_E = 0;
  std::vector<Float_t>* Cluster_Px = 0;
  std::vector<Float_t>* Cluster_Py = 0;
  std::vector<Float_t>* Cluster_Pz = 0; // TODO: Add M02
  // TODO: Add Jet stuff
  TreeBuffer(TChain* tree, GlobalOptions optns);
  ~TreeBuffer(){};
};

TreeBuffer::TreeBuffer(TChain* tree, GlobalOptions optns)
{
  NEvents = tree->GetEntries();
  switch (optns.TreeFormat) {
    case kRun2Tree:
      ReadRun2TreeIntoBuffer(tree, optns);
      break;
    case kRun3Tree:
    case kBerkeleyTree:
    default:
      FATAL(Form("Unknown treeFormat %d", optns.TreeFormat))
  }
}

void TreeBuffer::ReadRun2TreeIntoBuffer(TChain* tree, GlobalOptions optns)
{
  if (optns.doIsoGamma) {
    tree->SetBranchAddress("Cluster_E", &Cluster_E);
    tree->SetBranchAddress("Cluster_Px", &Cluster_Px);
    tree->SetBranchAddress("Cluster_Py", &Cluster_Py);
    tree->SetBranchAddress("Cluster_Pz", &Cluster_Pz);
  }
}

#endif