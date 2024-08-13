#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"

void makeHistosFromTree(TString AnalysisDirectory, int jobId = 0)
{
//General:
  ENTER

  if (jobId < 0)
    FATAL("Negative jobId")
  if(!jobId)
    LOG("This is a debug run")
  //Set global options
  GlobalOptions optns(AnalysisDirectory, jobId);
  //Create output file
  TFile *fOut = new TFile(Form("%s/HistosFromTree_%d.root", AnalysisDirectory.Data(), jobId), "RECREATE");
  //Set up directories for output
  TDirectory *hDirEvents = DefineEventHistograms(fOut, optns);
  TDirectory *hQADirEvents = DefineEventQAHistograms(fOut, optns);
  TDirectory *hDirClusters = DefineIsoGammaHistograms(fOut, "Clusters" ,optns);
  TDirectory *hDirGammas = DefineIsoGammaHistograms(fOut, "Gammas", optns);
  TDirectory *hDirIsoGammas = DefineIsoGammaHistograms(fOut,"IsoGammas" , optns);
  TDirectory *hQADirClusters = DefineIsoGammaQAHistograms(fOut, "ClusterQA", optns);
  TDirectory *hQADirGammas = DefineIsoGammaQAHistograms(fOut,"GammaQA", optns);
  TDirectory *hQADirIsoGammas = DefineIsoGammaQAHistograms(fOut,"IsoGammaQA", optns);
  TDirectory *hDirJets = DefineJetHistograms(fOut, optns);
  TDirectory *hQADirJets = DefineJetQAHistograms(fOut, optns);
  TDirectory *hDirGammaJetCorrelations = DefineGammaJetHistograms(fOut, optns);
  TDirectory *hQADirGammaJetCorrelations = DefineGammaJetQAHistograms(fOut, optns);
  //Supply the cuts:
  EventCuts eventCuts(optns);
  IsoGammaCuts isoGammaCuts(optns, hQADirIsoGammas); // Pass the QA dir so that the "cut passed" function can fill a "loss histogram"
  GammaGenCuts GammaGenCuts(optns);
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
  //Load input data into TChain.
  TChain *chain = readTree(Form("%s/../InputFiles/InputFiles_group_%d.txt", AnalysisDirectory.Data(), jobId));

  optns.TreeFormat = listTreeBranches(chain);

  TreeBuffer tree(chain, optns);
//Start Event loop
  for (int iEvent = 0; iEvent < tree.NEvents; iEvent++)
  {
    //Get event attributes from TChain:
    chain->GetEntry(iEvent);

    if (jobId == 0)
      PrintProgress(iEvent, tree.NEvents);
    else
      PrintProgressNumber(iEvent, tree.NEvents, 1000);
    
    Event event(tree, optns);
    //Apply event cuts and fill output hists:
    if (!eventCuts.PassedCuts(event))
      continue;
    fillHistograms(event, hDirEvents, event.weight, optns, isoGammaCuts);
    if (optns.doQA)
      fillQAHistograms(event, hQADirEvents, event.weight, optns, isoGammaCuts);
//If do Isogammas, calculate the relevant output.
    if (optns.doIsoGamma)
    {
      saveClustersFromEventInVector(tree, IsoGammas, optns);
      calculateIsolation(IsoGammas, event, isoGammaCuts.useRhoInsteadOfPerpCone);
      if(optns.isMC){
        //Load generated photons.
        saveGenPhotonsFromEventInVector(tree, GammaGens, optns);
        GammaGencalculateIsolation(GammaGens, event);
        
        //All generated photons.
        fillGammaGenHistograms(GammaGens, hDirGammas, event.weight, optns, GammaGenCuts);
        if(optns.doQA){
          fillGammaGenQAHistograms(GammaGens, hQADirGammas, event.weight, optns, GammaGenCuts);
        }
        //Isolated generated photons.
        doGammaGenIsolationCuts(GammaGens, GammaGenCuts);
        fillGammaGenHistograms(GammaGens, hDirIsoGammas, event.weight, optns, GammaGenCuts);
        if(optns.doQA){
          fillGammaGenQAHistograms(GammaGens, hQADirIsoGammas, event.weight, optns, GammaGenCuts);
        }
        
      }
      if (isoGammaCuts.applyNonLin)
        applyNonLinAndFineTuningCorrection(IsoGammas, isoGammaCuts, optns);
      pairIsoGammasFromEventInVector(IsoGammas, Pi0sForIsoGammaQA);
      //Fill hists (raw clusters)
      fillHistograms(IsoGammas, hDirClusters, event.weight, optns, isoGammaCuts);
      if (optns.doQA)
      {
        fillQAHistograms(IsoGammas, hQADirClusters, event.weight, optns, isoGammaCuts);
      }
      //Fill hists (cluster cuts, not isolated)
      doIsoGammaClusterCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirGammas, event.weight, optns, isoGammaCuts);
      if (optns.doQA)
      {
        fillQAHistograms(IsoGammas, hQADirGammas, event.weight, optns, isoGammaCuts);
      }
      //Fill hists (cluster cuts and isolated)
      doIsoGammaCuts(IsoGammas, isoGammaCuts);
      fillHistograms(IsoGammas, hDirIsoGammas, event.weight, optns, isoGammaCuts);
      if (optns.doQA)
      {
        fillQAHistograms(IsoGammas, hQADirIsoGammas, event.weight, optns, isoGammaCuts);
        fillQAHistograms(Pi0sForIsoGammaQA, hQADirIsoGammas, event.weight, optns, isoGammaCuts);
      }
    }
//Do jets:
    if (optns.doJets)
    {
      saveJetsFromEventInVector(tree, Jets, optns);
      doJetCuts(Jets, jetCuts);
      fillHistograms(Jets, hDirJets, event.weight, optns, isoGammaCuts);
      if (optns.doQA)
        fillQAHistograms(Jets, hQADirJets, event.weight, optns, isoGammaCuts);
      if (optns.isMC)
      {
        savePLJetsFromEventInVector(tree, PLJets, optns);
        mapPLtoDLjets(Jets, PLJets, jetCuts.R);
        // TODO: Do I need PL jet cuts?
        fillHistograms(PLJets, hDirJets, event.weight, optns, isoGammaCuts);
        if (optns.doQA)
          fillQAHistograms(PLJets, hQADirJets, event.weight, optns, isoGammaCuts);
      }
    }
    if (optns.doIsoGamma && optns.doJets)
    {
      pairGammasWithJets(IsoGammas, Jets, GammaJetPairs);
      // TODO: CorrelationCuts
      fillHistograms(GammaJetPairs, hDirGammaJetCorrelations, event.weight, optns, isoGammaCuts);
      if (optns.doQA)
        fillQAHistograms(GammaJetPairs, hQADirGammaJetCorrelations, event.weight, optns, isoGammaCuts);
    }
    IsoGammas.clear();
    GammaGens.clear();
    Jets.clear();
    PLJets.clear();
    Pi0s.clear();
    Pi0sForIsoGammaQA.clear();
    GammaJetPairs.clear();
  }
//End event loop:
//Write output
  fOut->Write();

  EXIT
}