#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "Utilities.h"
// TODO: Turn doIsoGamma, ... arguments into "Cutstrings" like Daikis O2Physics CutLibrary
// TODO: Add HistogramRegistry (as found in O2Physics)
void makeHistosFromTree(bool doQA, bool doIsoGamma = false, bool doJets = false, bool doPi0s = false)
{
  ENTER

  GlobalOptions optns(doQA, doIsoGamma, doJets, doPi0s);

  THashList hList = DefineHistograms();

  // These vectors store all information about all selected (by cuts) physics objects within a given event
  std::vector<IsoGamma> IsoGammas;
  std::vector<Jet> Jets;
  std::vector<Pi0> Pi0s;

  TChain* chain = readTree("Input/InputFiles.txt");
  // TODO: Read histograms from input file

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);

  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++) {
    chain->GetEntry(iEvent);

    // TODO: Read event as its own class

    fillEventHistograms();
    // TODO: Event cut here
    fillEventHistograms();

    if (optns.doIsoGamma) {
      saveClustersFromEventInVector(tree, IsoGammas);
      fillClusterHistograms(IsoGammas, hList, 0);
      doClusterCuts(IsoGammas);
      doIsolationCuts(IsoGammas);
      fillClusterHistograms(IsoGammas, hList, 1);
    }
    if (optns.doJets) {
        // TODO
    }
    if (optns.doIsoGamma && optns.doJets) {
        // TODO: Correlations
    }
  }

  EXIT
}