#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"
#include "JetCorrections.h"

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
  TDirectory *hDirJets = DefineJetHistograms(fOut, "Jets" ,optns);
  TDirectory *hDirJetsRaw = DefineJetHistograms(fOut, "JetsRaw" ,optns);
  TDirectory *hQADirJets = DefineJetQAHistograms(fOut, "JetQA" ,optns);
  TDirectory *hQADirJetsRaw = DefineJetQAHistograms(fOut, "JetQARaw" ,optns);
  TDirectory *hDirmPi0s = DefineIsoGammaHistograms(fOut, "mPi0s", optns);
  TDirectory *hQADirmPi0s = DefineIsoGammaQAHistograms(fOut, "mPi0QA", optns);
  TDirectory *hQADirgPi0s = DefineIsoGammaQAHistograms(fOut, "gammaForPi0QA", optns);
  TDirectory *hDirGammaJetCorrelations = DefineGammaJetHistograms(fOut, "GammaJetCorrelations", optns);
  TDirectory *hQADirGammaJetCorrelations = DefineGammaJetQAHistograms(fOut, "GammaJetCorrelationQA", optns);
  TDirectory *hDirmPi0JetCorrelations = DefineGammaJetHistograms(fOut, "mPi0JetCorrelations", optns);
  TDirectory *hQADirmPi0JetCorrelations = DefineGammaJetQAHistograms(fOut, "mPi0JetCorrelationQA", optns);
  TDirectory *hDirGGPi0JetCorrelations = DefineGammaJetHistograms(fOut, "ggPi0JetCorrelations", optns);
  TDirectory *hQADirGGPi0JetCorrelations = DefineGammaJetQAHistograms(fOut, "ggPi0JetCorrelationQA", optns);
  
  // histograms for exclusive trigger photons (all correlations done for signal and reference trigger)
  TDirectory *hDirTriggerPhotonsSignal = DefineIsoGammaHistograms(fOut, "TriggerPhotonsSignal", optns);
  TDirectory *hQADirTriggerPhotonsSignal = DefineIsoGammaQAHistograms(fOut, "TriggerPhotonSignalQA", optns);
  TDirectory *hDirTriggerPhotonsReference = DefineIsoGammaHistograms(fOut, "TriggerPhotonsReference", optns);
  TDirectory *hQADirTriggerPhotonsReference = DefineIsoGammaQAHistograms(fOut, "TriggerPhotonReferenceQA", optns);
  TDirectory *hDirTriggerGammaJetCorrelationsSignal = DefineGammaJetHistograms(fOut, "TriggerGammaJetCorrelationsSignal", optns);
  TDirectory *hQADirTriggerGammaJetCorrelationsSignal = DefineGammaJetQAHistograms(fOut, "TriggerGammaJetCorrelationSignalQA", optns);
  TDirectory *hDirTriggerGammaJetCorrelationsReference = DefineGammaJetHistograms(fOut, "TriggerGammaJetCorrelationsReference", optns);
  TDirectory *hQADirTriggerGammaJetCorrelationsReference = DefineGammaJetQAHistograms(fOut, "TriggerGammaJetCorrelationReferenceQA", optns);
  // histograms for jets with exclusive trigger photons
  TDirectory *hDirJetsWithPhotonSignalTrigger = DefineJetHistograms(fOut, "JetsWithPhotonSignalTrigger" ,optns);
  TDirectory *hDirJetsWithPhotonReferenceTrigger = DefineJetHistograms(fOut, "JetsWithPhotonReferenceTrigger" ,optns);
  TDirectory *hQADirJetsWithPhotonSignalTrigger = DefineJetQAHistograms(fOut, "JetQAWithPhotonSignalTrigger" ,optns);
  TDirectory *hQADirJetsWithPhotonReferenceTrigger = DefineJetQAHistograms(fOut, "JetQAWithPhotonReferenceTrigger" ,optns);

  // histograms for exclusive trigger merged pi0s
  TDirectory *hDirTriggerMergedPi0sSignal = DefineIsoGammaHistograms(fOut, "TriggerMergedPi0sSignal", optns);
  TDirectory *hQADirTriggerMergedPi0sSignal = DefineIsoGammaQAHistograms(fOut, "TriggerMergedPi0sSignalQA", optns);
  TDirectory *hDirTriggerMergedPi0sReference = DefineIsoGammaHistograms(fOut, "TriggerMergedPi0sReference", optns);
  TDirectory *hQADirTriggerMergedPi0sReference = DefineIsoGammaQAHistograms(fOut, "TriggerMergedPi0sReferenceQA", optns);
  TDirectory *hDirTriggerMergedPi0JetCorrelationsSignal = DefineGammaJetHistograms(fOut, "TriggerMergedPi0JetCorrelationsSignal", optns);
  TDirectory *hQADirTriggerMergedPi0JetCorrelationsSignal = DefineGammaJetQAHistograms(fOut, "TriggerMergedPi0JetCorrelationSignalQA", optns);
  TDirectory *hDirTriggerMergedPi0JetCorrelationsReference = DefineGammaJetHistograms(fOut, "TriggerMergedPi0JetCorrelationsReference", optns);
  TDirectory *hQADirTriggerMergedPi0JetCorrelationsReference = DefineGammaJetQAHistograms(fOut, "TriggerMergedPi0JetCorrelationReferenceQA", optns);

  TDirectory *hDirJetsWithMergedPi0SignalTrigger = DefineJetHistograms(fOut, "JetsWithMergedPi0SignalTrigger" ,optns);
  TDirectory *hDirJetsWithMergedPi0ReferenceTrigger = DefineJetHistograms(fOut, "JetsWithMergedPi0ReferenceTrigger" ,optns);
  TDirectory *hQADirJetsWithMergedPi0SignalTrigger = DefineJetQAHistograms(fOut, "JetQAWithMergedPi0SignalTrigger" ,optns);
  TDirectory *hQADirJetsWithMergedPi0ReferenceTrigger = DefineJetQAHistograms(fOut, "JetQAWithMergedPi0ReferenceTrigger" ,optns);

  // Load the cuts:
  EventCuts eventCuts(optns);
  IsoGammaCuts isoGammaCuts(optns, hQADirIsoGammas); // Pass the QA dir so that the "cut passed" function can fill a "loss histogram"
  GammaGenCuts GammaGenCuts(optns);
  DLJetCuts dljetCuts(optns);
  Pi0Cuts ggpi0Cuts(optns);
  IsoGammaCuts mPi0Cuts(optns, hQADirmPi0s);
  IsoGammaCuts gPi0Cuts(optns, hQADirgPi0s);
  // exclusive single trigger photon selection like in h-jet analysis
  ExclusiveTriggerParticleSelection exclusiveTriggerPhotonSelection(optns);
  // exclusive single trigger merged pi0 selection like in h-jet analysis
  // this needs to be a separate class to ensure also independence between mPi0 and photon trigger
  ExclusiveTriggerParticleSelection exclusiveTriggermPi0Selection(optns);

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
  // vectors for exclusive trigger particle selection
  std::vector<Cluster> triggerIsoPhotons;
  std::vector<Cluster> triggermPi0s;
  std::vector<GammaJetPair> triggerGammaJetPairs;
  std::vector<GammaJetPair> triggerMergedPi0JetPairs;
  Cluster trigIsoPhoton;
  Cluster trigmPi0;
  bool isSigEventGamma = false;
  bool isSigEventmPi0 = false;

  TChain *chain = readTree(Form("%s/../InputFiles/InputFiles_group_%d.txt", AnalysisDirectory.Data(), !jobId ? 1 : jobId), optns);

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);

  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++) // Start Event loop
  {
    chain->GetEntry(iEvent); // Get event attributes from TChain:

    PrintProgress(iEvent, tree.NEvents, 1000, jobId == 0);
    Event event(tree, optns);
    fillHistograms(event, hDirEvents, hQADirEvents, event.weight, eventCuts, optns);
    if (!eventCuts.PassedCuts(event))
      continue;

    // ###################### Isolated Gammas ######################
    if (optns.doIsoGamma)
    {
      saveClustersFromEventInVector(tree, IsoGammas, optns); // Load all clusters
      calculateIsolation(IsoGammas, event, isoGammaCuts.useRhoInsteadOfPerpCone, isoGammaCuts.IsolationConeRadius);
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
      fillHistograms(IsoGammas,event, hDirClusters, hQADirClusters, event.weight, optns, isoGammaCuts); // Fill hists (raw clusters)
      doIsoGammaClusterCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas,event, hDirGammas, hQADirGammas, event.weight, optns, isoGammaCuts); // Fill hists (cluster cuts, not isolated)
      doIsoGammaCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas,event, hDirIsoGammas, hQADirIsoGammas, event.weight, optns, isoGammaCuts); // Fill hists (cluster cuts and isolated)
    }


    // ###################### Jets ######################
    if (optns.doJets)
    {
      saveJetsFromEventInVector(tree, DLJets, optns, dljetCuts.R);
      if (dljetCuts.doUEsubtraction){
        fillHistograms(DLJets, event, hDirJetsRaw, hQADirJetsRaw, event.weight, optns);
        // UE event subtraction
        applyJetPtUESubtraction(DLJets, event, dljetCuts);
      }
      doJetCuts(DLJets, dljetCuts);
      fillHistograms(DLJets, event, hDirJets, hQADirJets, event.weight, optns);
      if (optns.isMC)
      {
        savePLJetsFromEventInVector(tree, PLJets, optns);
        mapPLtoDLjets(DLJets, PLJets, dljetCuts.R);
        // TODO: Do I need PL jet cuts?
        fillHistograms(PLJets, hDirJets, hQADirJets, event.weight, optns);
      }
    }

    // ###################### Exclusive Trigger Photon Particle Selection ######################
    if (exclusiveTriggerPhotonSelection.doExclusiveSelections)
    {
      isSigEventGamma = exclusiveTriggerPhotonSelection.isSignalEvent();
      bool found = exclusiveTriggerPhotonSelection.getTriggerParticle(IsoGammas, trigIsoPhoton);
      if(found)
        triggerIsoPhotons.push_back(trigIsoPhoton);
      if(isSigEventGamma)
        fillHistograms(triggerIsoPhotons,event, hDirTriggerPhotonsSignal, hQADirTriggerPhotonsSignal, event.weight, optns, isoGammaCuts);
      else
        fillHistograms(triggerIsoPhotons,event, hDirTriggerPhotonsReference, hQADirTriggerPhotonsReference, event.weight, optns, isoGammaCuts);
      
      // Fill new JET histograms that are only for events where a signal trigger particle was found
      // cuts were already applied in main track loop
      if(isSigEventGamma && found && optns.doJets){
        fillHistograms(DLJets, event, hDirJetsWithPhotonSignalTrigger, hQADirJetsWithPhotonSignalTrigger , event.weight, optns);
      }
      else if(!isSigEventGamma && found && optns.doJets){
        fillHistograms(DLJets, event, hDirJetsWithPhotonReferenceTrigger, hQADirJetsWithPhotonReferenceTrigger , event.weight, optns);
      }
    }

    // ###################### Merged Pi0s ######################
    if (optns.domPi0)
    {
      saveClustersFromEventInVector(tree, mPi0s, optns);
      doIsoGammaClusterCuts(mPi0s, mPi0Cuts);
      doIsoGammaCuts(mPi0s, mPi0Cuts);
      fillHistograms(mPi0s,event, hDirmPi0s, hQADirmPi0s, event.weight, optns, isoGammaCuts);
    }

    // ###################### Exclusive Trigger Merged Pi0 Particle Selection ######################
    if (exclusiveTriggermPi0Selection.doExclusiveSelections){
      isSigEventmPi0 = exclusiveTriggermPi0Selection.isSignalEvent();
      bool found = exclusiveTriggermPi0Selection.getTriggerParticle(mPi0s, trigmPi0);
      if(found)
        triggermPi0s.push_back(trigmPi0);
      if(isSigEventmPi0)
        fillHistograms(triggermPi0s,event, hDirTriggerMergedPi0sSignal, hQADirTriggerMergedPi0sSignal, event.weight, optns, isoGammaCuts);
      else
        fillHistograms(triggermPi0s,event, hDirTriggerMergedPi0sReference, hQADirTriggerMergedPi0sReference, event.weight, optns, isoGammaCuts);
      
      // Fill new JET histograms that are only for events where a signal trigger particle was found
      // cuts were already applied in main track loop
      if(isSigEventmPi0 && found && optns.doJets){
          fillHistograms(DLJets, event, hDirJetsWithMergedPi0SignalTrigger, hQADirJetsWithMergedPi0SignalTrigger , event.weight, optns);
      }
      else if(!isSigEventmPi0 && found && optns.doJets){
        fillHistograms(DLJets, event, hDirJetsWithMergedPi0ReferenceTrigger, hQADirJetsWithMergedPi0ReferenceTrigger , event.weight, optns);
      }
    }

    // ###################### GammaGamma Pi0s ######################
    if (optns.doGGPi0)
    {
      saveClustersFromEventInVector(tree, gPi0s, optns);
      doIsoGammaClusterCuts(gPi0s, gPi0Cuts);
      doIsoGammaShowerShapeCuts(gPi0s, gPi0Cuts);
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

    // ##################### Trigger Photon - Jet Correlations #####################
    // This is done for the exclusive selection with only one trigger particle per event
    if(exclusiveTriggerPhotonSelection.doExclusiveSelections){
      pairXWithYIntoZ(triggerIsoPhotons, DLJets, triggerGammaJetPairs);
      if(isSigEventGamma){
        fillHistograms(triggerGammaJetPairs, hDirTriggerGammaJetCorrelationsSignal, hQADirTriggerGammaJetCorrelationsSignal, event.weight, optns);
      }
      else{
        fillHistograms(triggerGammaJetPairs, hDirTriggerGammaJetCorrelationsReference, hQADirTriggerGammaJetCorrelationsReference, event.weight, optns);
      }

    }

    // ###################### Trigger Merged Pi0 - Jet Correlations ######################
    if(exclusiveTriggermPi0Selection.doExclusiveSelections){
      pairXWithYIntoZ(triggermPi0s, DLJets, triggerMergedPi0JetPairs);
      if(isSigEventmPi0){
        fillHistograms(triggerMergedPi0JetPairs, hDirTriggerMergedPi0JetCorrelationsSignal, hQADirTriggerMergedPi0JetCorrelationsSignal, event.weight, optns);
      }
      else{
        fillHistograms(triggerMergedPi0JetPairs, hDirTriggerMergedPi0JetCorrelationsReference, hQADirTriggerMergedPi0JetCorrelationsReference, event.weight, optns);
      }
    }

    // ###################### Clear vectors ######################

    IsoGammas.clear();
    GammaGens.clear();
    triggerIsoPhotons.clear();
    triggermPi0s.clear();
    DLJets.clear();
    PLJets.clear();
    mPi0s.clear();
    gPi0s.clear();
    ggPi0s.clear();
    GammaJetPairs.clear();
    GGPi0JetPairs.clear();
    mPi0JetPairs.clear();
    triggerGammaJetPairs.clear();
    triggerMergedPi0JetPairs.clear();
  } // End of the event loop

  fOut->Write(); // Write output

  EXIT
}