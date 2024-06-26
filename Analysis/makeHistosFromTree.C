#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"

void makeHistosFromTree(TString AnalysisDirectory, int jobId = 0)
{
  ENTER

  if (jobId < 0)
    FATAL("Negative jobId")
  if(!jobId)
    LOG("This is a debug run")

  GlobalOptions optns(AnalysisDirectory, jobId);

  TFile *fOut = new TFile(Form("%s/HistosFromTree_%d.root", AnalysisDirectory.Data(), jobId), "RECREATE");

  TDirectory *hDirEvents = DefineEventHistograms(fOut, optns);
  TDirectory *hQADirEvents = DefineEventQAHistograms(fOut, optns);
  TDirectory *hDirIsoGammasRaw = DefineIsoGammaHistograms(fOut, "Clusters" ,optns);
  TDirectory *hDirIsoGammasClusterCuts = DefineIsoGammaHistograms(fOut, "Gammas", optns);
  TDirectory *hDirIsoGammas = DefineIsoGammaHistograms(fOut,"IsoGammas" , optns);
  TDirectory *hQADirIsoGammas = DefineIsoGammaQAHistograms(fOut, optns);
  TDirectory *hDirJets = DefineJetHistograms(fOut, optns);
  TDirectory *hQADirJets = DefineJetQAHistograms(fOut, optns);
  TDirectory *hDirGammaJetCorrelations = DefineGammaJetHistograms(fOut, optns);
  TDirectory *hQADirGammaJetCorrelations = DefineGammaJetQAHistograms(fOut, optns);

  EventCuts eventCuts(optns);
  IsoGammaCuts isoGammaCuts(optns, hQADirIsoGammas); // Pass the QA dir so that the "cut passed" function can fill a "loss histogram"
  JetCuts jetCuts(optns);
  Pi0Cuts pi0Cuts(optns);

  // These vectors store all information about all selected (by cuts) physics objects within a given event
  std::vector<IsoGamma> IsoGammas;
  std::vector<GammaGen> GammaGens;
  std::vector<Pi0> Pi0sForIsoGammaQA;
  std::vector<Jet> Jets;
  std::vector<PLJet> PLJets; // Particle Level Jets -> Will only be filled if this is a MC
  std::vector<GammaJetPair> GammaJetPairs;
  std::vector<Pi0> Pi0s;

  TChain *chain = readTree(Form("%s/../InputFiles/InputFiles_group_%d.txt", AnalysisDirectory.Data(), jobId));

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);

  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++)
  {
    chain->GetEntry(iEvent);

    if (jobId == 0)
      PrintProgress(iEvent, tree.NEvents);
    else
      PrintProgressNumber(iEvent, tree.NEvents, 1000);

    Event event(tree, optns);

    if (!eventCuts.PassedCuts(event))
      continue;
    fillHistograms(event, hDirEvents, event.weight, optns);
    if (optns.doQA)
      fillQAHistograms(event, hQADirEvents, event.weight, optns);

    if (optns.doIsoGamma)
    {
      saveClustersFromEventInVector(tree, IsoGammas, optns);
      calculateIsolation(IsoGammas, event, isoGammaCuts.useRhoInsteadOfPerpCone);
      //fillnoCutsHistograms(IsoGammas, hDirIsoGammas, event.weight, optns);
      if(optns.isMC){
        saveGenPhotonsFromEventInVector(tree, GammaGens, optns);
        fillGammaGenHistograms(GammaGens, hDirIsoGammas, event.weight, optns);
      }
      if (isoGammaCuts.applyNonLin)
        applyNonLinAndFineTuningCorrection(IsoGammas, isoGammaCuts, optns);
      pairIsoGammasFromEventInVector(IsoGammas, Pi0sForIsoGammaQA);
      //Fill hists (cluster cuts, not isolated)
      doIsoGammaClusterCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirIsoGammasClusterCuts, event.weight, optns);
      //Fill hists (cluster cuts and isolated)
      doIsoGammaCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirIsoGammas, event.weight, optns);
      if (optns.doQA)
      {
        fillQAHistograms(IsoGammas, hQADirIsoGammas, event.weight, optns);
        fillQAHistograms(Pi0sForIsoGammaQA, hQADirIsoGammas, event.weight, optns);
      }
    }
    if (optns.doJets)
    {
      saveJetsFromEventInVector(tree, Jets);
      // TODO: JetCuts
      fillHistograms(Jets, hDirJets, event.weight, optns);
      if (optns.doQA)
        fillQAHistograms(Jets, hQADirJets, event.weight, optns);
      if (optns.isMC)
      {
        savePLJetsFromEventInVector(tree, PLJets);
        mapPLtoDLjets(Jets, PLJets, jetCuts.R);
        // TODO: Do I need PL jet cuts?
        fillHistograms(PLJets, hDirJets, event.weight, optns);
        if (optns.doQA)
          fillQAHistograms(PLJets, hQADirJets, event.weight, optns);
      }
    }
    if (optns.doIsoGamma && optns.doJets)
    {
      pairGammasWithJets(IsoGammas, Jets, GammaJetPairs);
      // TODO: CorrelationCuts
      fillHistograms(GammaJetPairs, hDirGammaJetCorrelations, event.weight, optns);
      if (optns.doQA)
        fillQAHistograms(GammaJetPairs, hQADirGammaJetCorrelations, event.weight, optns);
    }
    IsoGammas.clear();
    GammaGens.clear();
    Jets.clear();
    PLJets.clear();
    Pi0s.clear();
    Pi0sForIsoGammaQA.clear();
    GammaJetPairs.clear();
  }

  fOut->Write();

  EXIT
}