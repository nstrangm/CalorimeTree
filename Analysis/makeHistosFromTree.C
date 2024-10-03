#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"

void makeHistosFromTree(TString AnalysisDirectory, int jobId = 0)
{
  ENTER

  GlobalOptions optns(AnalysisDirectory, jobId);
  // Create output file
  TFile *fOut = new TFile(Form("%s/HistosFromTree_%d.root", AnalysisDirectory.Data(), jobId), "RECREATE");

  // Set up directories for output
  TDirectory *hDirEvents = DefineEventHistograms(fOut, optns);
  TDirectory *hQADirEvents = DefineEventQAHistograms(fOut, optns);
  TDirectory *hDirClusters = DefineIsoGammaHistograms(fOut, "Clusters", optns);
  TDirectory *hQADirClusters = DefineIsoGammaQAHistograms(fOut, "ClusterQA", optns);
  TDirectory *hDirGammas = DefineIsoGammaHistograms(fOut, "Gammas", optns);
  TDirectory *hQADirGammas = DefineIsoGammaQAHistograms(fOut, "GammaQA", optns);
  TDirectory *hDirIsoGammas = DefineIsoGammaHistograms(fOut, "IsoGammas", optns);
  TDirectory *hQADirIsoGammas = DefineIsoGammaQAHistograms(fOut, "IsoGammaQA", optns);
  TDirectory *hDirggPi0s = DefineGGPi0Histograms(fOut, "ggPi0", optns);
  TDirectory *hQADirggPi0s = DefineGGPi0QAHistograms(fOut, "ggPi0QA", optns);
  TDirectory *hDirJets = DefineJetHistograms(fOut, optns);
  TDirectory *hQADirJets = DefineJetQAHistograms(fOut, optns);
  TDirectory *hDirmPi0s = DefineIsoGammaHistograms(fOut, "mPi0s", optns);
  TDirectory *hQADirmPi0s = DefineIsoGammaQAHistograms(fOut, "mPi0QA", optns);
  TDirectory *hDirGammaJetCorrelations = DefineGammaJetHistograms(fOut, "GammaJetCorrelations", optns);
  TDirectory *hQADirGammaJetCorrelations = DefineGammaJetQAHistograms(fOut, "GammaJetCorrelationQA", optns);
  TDirectory *hDirmPi0JetCorrelations = DefineGammaJetHistograms(fOut, "mPi0JetCorrelations", optns);
  TDirectory *hQADirmPi0JetCorrelations = DefineGammaJetQAHistograms(fOut, "mPi0JetCorrelationQA", optns);
  TDirectory *hDirGGPi0JetCorrelations = DefineGammaJetHistograms(fOut, "ggPi0JetCorrelations", optns);
  TDirectory *hQADirGGPi0JetCorrelations = DefineGammaJetQAHistograms(fOut, "ggPi0JetCorrelationQA", optns);

  // Load the cuts:
  EventCuts eventCuts(optns);
  IsoGammaCuts isoGammaCuts(optns, hQADirIsoGammas); // Pass the QA dir so that the "cut passed" function can fill a "loss histogram"
  GammaGenCuts GammaGenCuts(optns);
  DLJetCuts dljetCuts(optns);
  Pi0Cuts ggpi0Cuts(optns);
  IsoGammaCuts mPi0Cuts(optns, hQADirmPi0s);

  // These vectors store all information about all selected (by cuts) physics objects within a given event
  std::vector<Cluster> IsoGammas;
  std::vector<GammaGen> GammaGens;
  std::vector<DLJet> DLJets;
  std::vector<PLJet> PLJets; // Particle Level Jets -> Will only be filled if this is a MC
  std::vector<GammaJetPair> GammaJetPairs;
  std::vector<GGPi0JetPair> GGPi0JetPairs;
  std::vector<GammaJetPair> mPi0JetPairs;
  std::vector<Pi0> ggPi0s;
  std::vector<Cluster> mPi0s;
  std::vector<Cluster> gPi0s;

  TChain *chain = readTree(Form("%s/../InputFiles/InputFiles_group_%d.txt", AnalysisDirectory.Data(), !jobId ? 1 : jobId), optns);

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);

  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++) // Start Event loop
  {
    chain->GetEntry(iEvent); // Get event attributes from TChain:

    PrintProgress(iEvent, tree.NEvents, 1000, jobId == 0);

    Event event(tree, optns);
    if (!eventCuts.PassedCuts(event))
      continue;
    fillHistograms(event, hDirEvents, hQADirEvents, event.weight, optns);

    // ###################### Isolated Gammas ######################
    if (optns.doIsoGamma)
    {
      saveClustersFromEventInVector(tree, IsoGammas, optns); // Load all clusters
      calculateIsolation(IsoGammas, event, isoGammaCuts.useRhoInsteadOfPerpCone);
      if (optns.isMC)
      {
        saveGenPhotonsFromEventInVector(tree, GammaGens, optns); // Load generated photons.
        GammaGencalculateIsolation(GammaGens, event);

        fillHistograms(GammaGens, hDirGammas, hQADirGammas, event.weight, optns, GammaGenCuts); // All generated photons.
        doGammaGenIsolationCuts(GammaGens, GammaGenCuts);
        fillHistograms(GammaGens, hDirIsoGammas, hQADirIsoGammas, event.weight, optns, GammaGenCuts); // Isolated generated photons.
      }
      if (isoGammaCuts.applyNonLin)
        applyNonLinAndFineTuningCorrection(IsoGammas, isoGammaCuts, optns);
      fillHistograms(IsoGammas, hDirClusters, hQADirClusters, event.weight, optns, isoGammaCuts); // Fill hists (raw clusters)
      doIsoGammaClusterCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirGammas, hQADirGammas, event.weight, optns, isoGammaCuts); // Fill hists (cluster cuts, not isolated)
      doIsoGammaCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirIsoGammas, hQADirIsoGammas, event.weight, optns, isoGammaCuts); // Fill hists (cluster cuts and isolated)
    }

    // ###################### Jets ######################
    if (optns.doJets)
    {
      saveJetsFromEventInVector(tree, DLJets, optns);
      doJetCuts(DLJets, dljetCuts);
      fillHistograms(DLJets, hDirJets, hQADirJets, event.weight, optns);
      if (optns.isMC)
      {
        savePLJetsFromEventInVector(tree, PLJets, optns);
        mapPLtoDLjets(DLJets, PLJets, dljetCuts.R);
        // TODO: Do I need PL jet cuts?
        fillHistograms(PLJets, hDirJets, hQADirJets, event.weight, optns);
      }
    }

    // ###################### Merged Pi0s ######################
    if (optns.domPi0)
    {
      saveClustersFromEventInVector(tree, mPi0s, optns);
      doIsoGammaClusterCuts(mPi0s, mPi0Cuts);
      doIsoGammaCuts(mPi0s, mPi0Cuts);
      fillHistograms(mPi0s, hDirmPi0s, hQADirmPi0s, event.weight, optns, isoGammaCuts);
    }

    // ###################### GammaGamma Pi0s ######################
    if (optns.doGGPi0)
    {
      saveClustersFromEventInVector(tree, gPi0s, optns);
      pairGammasFromEventInVector(gPi0s, ggPi0s);
      doPi0Cuts(ggPi0s, ggpi0Cuts);
      fillHistograms(ggPi0s, hDirggPi0s, nullptr, event.weight, optns);
    }

    // ###################### Isolated Gamma - Jet Correlations ######################
    if (optns.doIsoGamma && optns.doJets)
    {
      pairXWithYIntoZ(IsoGammas, DLJets, GammaJetPairs);
      // TODO: CorrelationCuts
      fillHistograms(GammaJetPairs, hDirGammaJetCorrelations, hQADirGammaJetCorrelations, event.weight, optns);
    }

    // ###################### GammaGamma Pi0 - Jet Correlations ######################
    if (optns.doGGPi0 && optns.doJets)
    {
      pairXWithYIntoZ(ggPi0s, DLJets, GGPi0JetPairs);
      // TODO: CorrelationCuts
      fillHistograms(GGPi0JetPairs, hDirGGPi0JetCorrelations, hQADirGGPi0JetCorrelations, event.weight, optns);
    }

    // ###################### Merged Pi0 - Jet Correlations ######################
    if (optns.domPi0 && optns.doJets)
    {
      pairXWithYIntoZ(mPi0s, DLJets, mPi0JetPairs);
      // TODO: CorrelationCuts
      fillHistograms(mPi0JetPairs, hDirmPi0JetCorrelations, hQADirmPi0JetCorrelations, event.weight, optns);
    }

    IsoGammas.clear();
    GammaGens.clear();
    DLJets.clear();
    PLJets.clear();
    mPi0s.clear();
    gPi0s.clear();
    ggPi0s.clear();
    GammaJetPairs.clear();
    GGPi0JetPairs.clear();
    mPi0JetPairs.clear();
  } // End of the event loop

  fOut->Write(); // Write output

  EXIT
}