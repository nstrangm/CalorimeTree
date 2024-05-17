#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"


void makeHistosFromTree(bool isMC, bool doQA, TString EventCutString = "", TString IsoGammaCutString = "", TString JetCutString = "", TString Pi0CutString = "")
{
  ENTER

  GlobalOptions optns(isMC, doQA, EventCutString, IsoGammaCutString, JetCutString, Pi0CutString);

  EventCuts eventCuts(EventCutString);
  IsoGammaCuts isoGammaCuts(IsoGammaCutString);
  JetCuts jetCuts(JetCutString);
  Pi0Cuts pi0Cuts(Pi0CutString);

  TFile* fOut = new TFile("Output/CalorimeTree.root", "RECREATE");

  TDirectory* hDirEvents = DefineEventHistograms(fOut);
  TDirectory* hQADirEvents = doQA ? DefineEventQAHistograms(fOut) : nullptr;
  TDirectory* hDirIsoGammas = DefineIsoGammaHistograms(fOut);
  TDirectory* hQADirIsoGammas = doQA ? DefineIsoGammaQAHistograms(fOut, optns) : nullptr;
  TDirectory* hDirJets = DefineJetHistograms(fOut);
  TDirectory* hQADirJets = doQA ? DefineJetQAHistograms(fOut) : nullptr;


  // These vectors store all information about all selected (by cuts) physics objects within a given event
  std::vector<IsoGamma> IsoGammas;
  std::vector<Pi0> Pi0sForIsoGammaQA;
  std::vector<Jet> Jets;
  std::vector<PLJet> PLJets; // Particle Level Jets -> Will only be filled if this is a MC
  std::vector<Pi0> Pi0s;

  TChain* chain = readTree("Input/InputFiles.txt");
  // TODO: Read histograms from input file

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);

  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++) {
    chain->GetEntry(iEvent);

    PrintProgress(iEvent, tree.NEvents);

    Event event(tree, optns);

    if (!eventCuts.PassedCuts(event))
      continue;
    fillHistograms(event, hDirEvents, event.weight);
    if (doQA)
      fillQAHistograms(event, hQADirEvents, event.weight, optns);

    if (optns.doIsoGamma) {
      saveClustersFromEventInVector(tree, IsoGammas, optns);
      applyNonLinAndFineTuningCorrection(IsoGammas, isoGammaCuts, optns);
      pairIsoGammasFromEventInVector(IsoGammas, Pi0sForIsoGammaQA);
      calculateIsolation(IsoGammas, event, isoGammaCuts.useRhoInsteadOfPerpCone);
      doIsoGammaCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirIsoGammas, event.weight);
      if (doQA){
        fillQAHistograms(IsoGammas, hQADirIsoGammas, event.weight, optns);
        fillQAHistograms(Pi0sForIsoGammaQA, hQADirIsoGammas, event.weight, optns);
      }
    }
    if (optns.doJets) {
      saveJetsFromEventInVector(tree, Jets);
      // TODO: JetCuts
      fillHistograms(Jets, hDirJets, event.weight);
      if (doQA)
        fillQAHistograms(Jets, hQADirJets, event.weight, optns);
      if (optns.isMC) {
        savePLJetsFromEventInVector(tree, PLJets);
        mapPLtoDLjets(Jets, PLJets, jetCuts.R);
        // TODO: Do I need PL jet cuts?
        fillHistograms(PLJets, hDirJets, event.weight);
        if (doQA)
          fillQAHistograms(PLJets, hQADirJets, event.weight, optns);
      }
    }
    if (optns.doIsoGamma && optns.doJets) {
      // TODO: Correlations
    }
    IsoGammas.clear();
    Jets.clear();
    PLJets.clear();
    Pi0s.clear();
    Pi0sForIsoGammaQA.clear();
  }

  fOut->Write();

  EXIT
}