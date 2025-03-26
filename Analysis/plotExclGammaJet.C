#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "plotExclGammaJet.h"
#include "analyseExclGammaJet.h"

Int_t fancyColors[8];

void extractConfigurationInformation(TString AnalysisDirectory, GlobalOptions &optns)
{
  // Extract centrality information from directory name
  TString trainConfig = AnalysisDirectory;
  if (trainConfig.Contains("Run3-")) {
    trainConfig = trainConfig(trainConfig.Index("Run3-"), trainConfig.Length());
    trainConfig = trainConfig(0, trainConfig.First("/"));
    
    // Run3-0-10
    INFO(Form("Train config: %s", trainConfig.Data()));
    // tokenize the Runstring with dashes
    TObjArray *runTokens = trainConfig.Tokenize("-");
    centralityRange.first = ((TObjString*)runTokens->At(1))->GetString().Atof();
    centralityRange.second = ((TObjString*)runTokens->At(2))->GetString().Atof();
    delete runTokens;
    
    DLJetCuts dljetCuts(optns);
    ExclusiveTriggerParticleSelection exclusiveTriggerPhotonSelection(optns);

    jetRadius = dljetCuts.R;
    signalTriggerPt = std::make_pair(exclusiveTriggerPhotonSelection.getPtMinSignal(), exclusiveTriggerPhotonSelection.getPtMaxSignal());
    referenceTriggerPt = std::make_pair(exclusiveTriggerPhotonSelection.getPtMinReference(), exclusiveTriggerPhotonSelection.getPtMaxReference());


    INFO("Extracted configuration information:");
    INFO(Form("Centrality range: %f - %f", centralityRange.first, centralityRange.second));
    INFO(Form("Jet radius: %.2f", jetRadius));
    INFO(Form("Signal trigger pt: %.1f - %.1f", signalTriggerPt.first, signalTriggerPt.second));
    INFO(Form("Reference trigger pt: %.1f - %.1f", referenceTriggerPt.first, referenceTriggerPt.second));
  }

}

void plotRecoilJets(TDirectory *dRecoilJets, TString triggerParticle, GlobalOptions optns)
{
   DLJetCuts jetCuts(optns);
   // create output subdirectory
   TString savePath = Form("%s/%sTrigger", outputDir.Data(), triggerParticle.Data());
   createDirectory(savePath.Data());

   TString triggerLabel = "";
   if (triggerParticle == "IsoPhoton") triggerLabel = "#gamma_{iso}";
   else if (triggerParticle == "MergedPi0") triggerLabel = "#pi^{0}_{merged}";
   
   TH2F* hIsoGammaJetDeltaPhiJetPtSignal = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtSignal");
   TH2F* hIsoGammaJetDeltaPhiJetPtReference = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtReference");

   // signal trigger
   TH2F* hIsoGammaJetDeltaPhi2piJetPtSignal = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhi2piJetPtSignal");
   Plotting2D PDeltaPhiJetPtSignal;
   PDeltaPhiJetPtSignal.New(hIsoGammaJetDeltaPhi2piJetPtSignal);
   PDeltaPhiJetPtSignal.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtSignal.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; signal %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtSignal.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtSignal.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignal.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignal.Plot(Form("%s/DeltaPhiJetPtSignal.%s", savePath.Data(), suffix.Data()), false, false, true);

   // reference trigger
   TH2F* hIsoGammaJetDeltaPhi2piJetPtReference = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhi2piJetPtReference");
   Plotting2D PDeltaPhiJetPtReference;
   PDeltaPhiJetPtReference.New(hIsoGammaJetDeltaPhi2piJetPtReference);
   PDeltaPhiJetPtReference.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtReference.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; reference %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtReference.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtReference.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReference.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReference.Plot(Form("%s/DeltaPhiJetPtReference.%s", savePath.Data(), suffix.Data()), false, false, true);

   // mirrored histograms - signal trigger
   TH2F* hIsoGammaJetDeltaPhiJetPtSignalMirrored = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtSignalMirrored");
   Plotting2D PDeltaPhiJetPtSignalMirrored;
   PDeltaPhiJetPtSignalMirrored.New(hIsoGammaJetDeltaPhiJetPtSignalMirrored);
   PDeltaPhiJetPtSignalMirrored.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtSignalMirrored.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; signal %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtSignalMirrored.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtSignalMirrored.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignalMirrored.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignalMirrored.Plot(Form("%s/DeltaPhiJetPtSignalMirrored.%s", savePath.Data(), suffix.Data()), false, false, true);

   // mirrored histograms - reference trigger
   TH2F* hIsoGammaJetDeltaPhiJetPtReferenceMirrored = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtReferenceMirrored");
   Plotting2D PDeltaPhiJetPtReferenceMirrored;
   PDeltaPhiJetPtReferenceMirrored.New(hIsoGammaJetDeltaPhiJetPtReferenceMirrored);
   PDeltaPhiJetPtReferenceMirrored.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtReferenceMirrored.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; reference %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtReferenceMirrored.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtReferenceMirrored.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReferenceMirrored.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReferenceMirrored.Plot(Form("%s/DeltaPhiJetPtReferenceMirrored.%s", savePath.Data(), suffix.Data()), false, false, true);


   // do the same plots but projected onto deltaPhi in bins of jet pt
   // Create PlottingGrid for signal trigger deltaPhi projections
   PlottingGrid PDeltaPhiJetPtGridSignal;
   PDeltaPhiJetPtGridSignal.SetMargins(0.12, 0.12, 0.05, 0.125);
   PDeltaPhiJetPtGridSignal.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   TH1F* hRecoilJetdPhiSignal[nJetPtBins];
   TH1F* hRecoilJetdPhiReference[nJetPtBins];
   for (int i = 0; i < nJetPtBins; i++)
   {
      int minBin = hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->FindBin(jetPtBins[i]);
      int maxBin = hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->FindBin(jetPtBins[i+1]);
      hRecoilJetdPhiSignal[i] = (TH1F*) hIsoGammaJetDeltaPhiJetPtSignal->ProjectionX(Form("hRecoilJetdPhiSignal_%d", i), minBin, maxBin);
      PDeltaPhiJetPtGridSignal.New(hRecoilJetdPhiSignal[i], "", -1, 1, kBlack, "hist e");
      PDeltaPhiJetPtGridSignal.NextPad(Form("%.1f < #it{p}_{T}^{jet} < %.1f GeV/#it{c}", jetPtBins[i], jetPtBins[i+1]));
   }
   PDeltaPhiJetPtGridSignal.Plot(Form("%s/DeltaPhiJetPtGridSignal.%s", savePath.Data(), suffix.Data()), false, false);

   // Create PlottingGrid for reference trigger deltaPhi projections
   PlottingGrid PDeltaPhiJetPtGridReference;
   PDeltaPhiJetPtGridReference.SetMargins(0.12, 0.12, 0.05, 0.125);
   PDeltaPhiJetPtGridReference.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   for (int i = 0; i < nJetPtBins; i++)
   {
      int minBin = hIsoGammaJetDeltaPhiJetPtReference->GetYaxis()->FindBin(jetPtBins[i]);
      int maxBin = hIsoGammaJetDeltaPhiJetPtReference->GetYaxis()->FindBin(jetPtBins[i+1]);
      hRecoilJetdPhiReference[i] = (TH1F*) hIsoGammaJetDeltaPhiJetPtReference->ProjectionX(Form("hRecoilJetdPhiReference_%d", i), minBin, maxBin);
      PDeltaPhiJetPtGridReference.New(hRecoilJetdPhiReference[i], "", -1, 1, kBlack, "hist e");
      PDeltaPhiJetPtGridReference.NextPad(Form("%.1f < #it{p}_{T}^{jet} < %.1f GeV/#it{c}", jetPtBins[i], jetPtBins[i+1]));
   }
   PDeltaPhiJetPtGridReference.Plot(Form("%s/DeltaPhiJetPtGridReference.%s", savePath.Data(), suffix.Data()), false, false);

   // loop over all phi bins in phiBins and get the plots
   for(auto phiBin : phiBins)
   {
     TString phiBinLabel = Form("#Delta#phi_{%.1f#pi,%.1f#pi}", phiBin.first/(2*TMath::Pi()), phiBin.second/(2*TMath::Pi()));
     TH1F* hRecoilJetPtSignal = (TH1F*) dRecoilJets->Get(Form("hRecoilJetPtSignal_%.2f_%.2f", phiBin.first, phiBin.second));
     TH1F* hRecoilJetPtReference = (TH1F*) dRecoilJets->Get(Form("hRecoilJetPtReference_%.2f_%.2f", phiBin.first, phiBin.second));
     TH1F* hRecoilJetPtSignalReferenceRatio = (TH1F*) dRecoilJets->Get(Form("hRecoilJetPtSignalReferenceRatio_%.2f_%.2f", phiBin.first, phiBin.second));
     TF1*  fCRefFit = (TF1*) dRecoilJets->Get(Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second));
      double phiMin = phiBin.first;
      double phiMax = phiBin.second;

      // if phiMin is close to pi-0.6 and phiMax is close to pi+0.6, set phiMin to pi-0.6 and phiMax to pi+0.6
      TString phiLabel = Form("%.2f < #varphi < %.2f", phiMin,phiMax);
      // check if phiMin is close to pi-0.6
      if(TMath::Abs(phiMin - (TMath::Pi() - 0.6)) < 0.05)
      {
        phiMin = TMath::Pi() - 0.6;
        phiLabel = "|#varphi - #pi| < 0.6";
      }

      double maxY = hRecoilJetPtReference->GetMaximum();
      if (hRecoilJetPtSignal->GetMaximum() > maxY) maxY = hRecoilJetPtSignal->GetMaximum();
      PlottingRatio PRecoilJets;
      PRecoilJets.SetMargins(0.1, 0.2, 0.025, 0.025, 2000, 1500, 0.33);
      PRecoilJets.NewRatioLine(-50, 1., 200., 1., 2, 2, kGray + 2);

      PRecoilJets.New(hRecoilJetPtSignal, Form("sig %s TT{%.0f,%.0f}", triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second), 20, 2, kRed, "p");                                               
      PRecoilJets.New(hRecoilJetPtReference, Form("ref %s TT{%.0f,%.0f}", triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second), 4, 2, kBlack, "p");
      PRecoilJets.NewRatio(hRecoilJetPtSignalReferenceRatio, "", 20, 2, kBlack, "p");                                          // -> Statistical Uncertainty
      PRecoilJets.NewRatio(fCRefFit, Form("#it{c}_{ref} = %.3f #pm %.3f", fCRefFit->GetParameter(0), fCRefFit->GetParError(0)), 20, 2, fancyColors[0], "l");
      PRecoilJets.SetLegend(0.5, 0.8, 0.6, 0.7, true);
      PRecoilJets.SetLegendR(0.65, 0.9, 0.3, 0.45, true);
      PRecoilJets.SetAxisRange(-50,201,8E-4,maxY*1.8,0.,2.6);
      PRecoilJets.NewLatex(0.9, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger; %s", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), phiLabel.Data()));
      PRecoilJets.SetAxisLabel("#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi}", "signal / reference", 1, 1.4);
      PRecoilJets.Plot(Form("%s/RecoilJetPt_%.2f_%.2f.%s", savePath.Data(), phiBin.first, phiBin.second, suffix.Data()), false, true);
   }

   // plot the 2D histogram hIsoGammaJetDeltaEtaJetPt
   TH2F* hIsoGammaJetDeltaEtaJetPt = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaEtaJetPt");
   Plotting2D PDeltaEtaJetPt;
   PDeltaEtaJetPt.New(hIsoGammaJetDeltaEtaJetPt);
   PDeltaEtaJetPt.SetAxisLabel("#Delta#eta", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaEtaJetPt.Plot(Form("%s/DeltaEtaJetPt.%s", savePath.Data(), suffix.Data()), false, false, true);

  //  Create a PlottingGrid that shows the projection onto deltaEta for different jet pt
  PlottingGrid PDeltaEtaJetPtGrid;
   const int nJetPtBins = 7;
   double jetPtBins[nJetPtBins+1] = {-50, 0, 5,10, 20, 30,40,50};
   TH1F* hDeltaEtaProjections[nJetPtBins];

   PDeltaEtaJetPtGrid.SetAxisLabel("#bf{#Delta#eta}", "#bf{N_{jets}}");
   PDeltaEtaJetPtGrid.SetMargins(0.12, 0.12, 0.05, 0.125);

   for (int iPtBin = 0; iPtBin < nJetPtBins; iPtBin++) {
     int minBin = hIsoGammaJetDeltaEtaJetPt->GetYaxis()->FindBin(jetPtBins[iPtBin]);
     int maxBin = hIsoGammaJetDeltaEtaJetPt->GetYaxis()->FindBin(jetPtBins[iPtBin + 1]);
     
     hDeltaEtaProjections[iPtBin] = (TH1F*)hIsoGammaJetDeltaEtaJetPt->ProjectionX(Form("hDeltaEtaProjection_%d", iPtBin), minBin, maxBin);
     hDeltaEtaProjections[iPtBin]->SetName(Form("hDeltaEtaProjection_%d", iPtBin));

     PDeltaEtaJetPtGrid.New(hDeltaEtaProjections[iPtBin], "", -1, 1, kBlack, "hist e");
     PDeltaEtaJetPtGrid.NextPad(Form("%.1f < #it{p}_{T}^{jet} < %.1f GeV/#it{c}", jetPtBins[iPtBin], jetPtBins[iPtBin + 1]));
   }

   PDeltaEtaJetPtGrid.SetLegend(0.1,0.9,0.1,0.4);
   PDeltaEtaJetPtGrid.NewLatex(0.9, 0.9, Form("#bf{ALICE work in progress};Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data()), StdTextSize*0.5, 0.06);

   PDeltaEtaJetPtGrid.Plot(Form("%s/DeltaEtaProjections.%s", savePath.Data(), suffix.Data()), false, false);




   TH1F* hRecoilJetPtDEta[10];
   TH1F* hRecoilJetPtDEtaRebinned[10];
   // get recoil jet pt for different dEta bins
   PlottingRatio PRecoilJetPtDEta;
   PRecoilJetPtDEta.SetMargins(0.1, 0.2, 0.025, 0.025, 2000, 1500, 0.33);
   PRecoilJetPtDEta.NewRatioLine(-50, 1., 200., 1., 2, 2, kGray + 2);

   double maxY = 0;

   int iEtaBin = 0;

   // get the correct dEta bins ( we need to adapt top boundary according )
   std::vector<std::pair<double, double>> dEtaBinsCorr;
   for (auto etaBin : dEtaBins)
   {
     float dEtaMax = jetCuts.JetEtaPhiMinMax[0][1] + 0.67;
      if(etaBin.second > dEtaMax)
      {
        etaBin.second = dEtaMax;
      }
      dEtaBinsCorr.emplace_back(etaBin.first,etaBin.second);
   }


   for (auto etaBin : dEtaBinsCorr)
   {
     hRecoilJetPtDEta[iEtaBin] = (TH1F*) dRecoilJets->Get(Form("hRecoilJetPtDEta_%.2f_%.2f", etaBin.first, etaBin.second));
     hRecoilJetPtDEta[iEtaBin]->SetName(Form("hRecoilJetPtDEta_%.2f_%.2f", etaBin.first, etaBin.second));
     hRecoilJetPtDEta[iEtaBin]->GetXaxis()->SetTitle("#it{p}_{T,jet}^{reco} (GeV/#it{c})");
     hRecoilJetPtDEta[iEtaBin]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{2}N}{d#it{p}_{T}^{jet} d#Delta#phi}");
     if (hRecoilJetPtDEta[iEtaBin]->GetMaximum() > maxY) maxY = hRecoilJetPtDEta[iEtaBin]->GetMaximum();
    
     hRecoilJetPtDEtaRebinned[iEtaBin] = (TH1F*) hRecoilJetPtDEta[iEtaBin]->Clone(Form("hRecoilJetPtDEtaRebinned_%.2f_%.2f", etaBin.first, etaBin.second));
     hRecoilJetPtDEtaRebinned[iEtaBin]->Rebin(4); // even more rebin
     iEtaBin++;
   }

   TH1F* hRatiosToLargestdEta[dEtaBinsCorr.size()-1];
   TH1F* hRatiosToLargestdEtaRebinned[dEtaBinsCorr.size()-1];
   
    // Divide all dEta bins by largest dEta bin
    for (int i = 0; i < dEtaBinsCorr.size()-1; i++) {
        hRatiosToLargestdEta[i] = (TH1F*) hRecoilJetPtDEta[dEtaBinsCorr.size()-1]->Clone(Form("hRatiosToLargestdEta_%d", i));
        hRatiosToLargestdEta[i]->Divide(hRecoilJetPtDEta[i], hRecoilJetPtDEta[dEtaBinsCorr.size()-1], 1.0, 1.0);
        // hRatiosToLargestdEta[i]->GetYaxis()->SetRangeUser(0.5, 4.);
        hRatiosToLargestdEtaRebinned[i] = (TH1F*) hRecoilJetPtDEtaRebinned[dEtaBinsCorr.size()-1]->Clone(Form("hRatiosToLargestdEtaRebinned_%d", i));
        hRatiosToLargestdEtaRebinned[i]->Divide(hRecoilJetPtDEtaRebinned[i], hRecoilJetPtDEtaRebinned[dEtaBinsCorr.size()-1], 1.0, 1.0);
        hRatiosToLargestdEtaRebinned[i]->GetYaxis()->SetRangeUser(0.5, 3.);
    }

    // Plot dEta bins
    PRecoilJetPtDEta.SetLegend(0.5, 0.8, 0.6, 0.7, true);
    PRecoilJetPtDEta.SetLegendR(0.65, 0.9, 0.3, 0.45, true);
    PRecoilJetPtDEta.SetAxisRange(-50, 201, 8E-4, maxY*1.8, 0., 10.1);
    PRecoilJetPtDEta.NewLatex(0.9, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data()));
    PRecoilJetPtDEta.SetAxisLabel("#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi d#Delta#eta}", "ratio to largest #Delta#eta", 1, 1.4);

    // Add histograms to plot
    for (int i = 0; i < dEtaBinsCorr.size(); i++) {
        PRecoilJetPtDEta.New(hRecoilJetPtDEta[i], Form("%.2f < #Delta#eta < %.2f", dEtaBinsCorr[i].first, dEtaBinsCorr[i].second), 20+i, 2, fancyColors[i], "p");
    }

    // Plot ratios
    for (int i = 0; i < dEtaBinsCorr.size()-1; i++) {
        PRecoilJetPtDEta.NewRatio(hRatiosToLargestdEta[i], Form("%.2f < #Delta#eta < %.2f", dEtaBinsCorr[i].first, dEtaBinsCorr[i].second), 20+i, 2, fancyColors[i], "p");
    }

    PRecoilJetPtDEta.Plot(Form("%s/RecoilJetPtDEta.%s", savePath.Data(), suffix.Data()), false, true);
    

    //Do a 1D plot of only the ratios
    Plotting1D PRecoilJetPtDEtaRatios;
    for (int i = 0; i < dEtaBinsCorr.size()-1; i++) {
        PRecoilJetPtDEtaRatios.New(hRatiosToLargestdEtaRebinned[i], Form("%.2f < #Delta#eta < %.2f", dEtaBinsCorr[i].first, dEtaBinsCorr[i].second), 20+i, 2, fancyColors[i], "p");
    }
    PRecoilJetPtDEtaRatios.Plot(Form("%s/RecoilJetPtDEtaRatiosRebinned.%s", savePath.Data(), suffix.Data()), false, false);
   
    
    // Plot recoil jet yield as a function of phi
    TGraphAsymmErrors* gRecoilJetPtIntegralPhi[recoilJetPtBinsIntegration.size()];
    for (int i = 0; i < recoilJetPtBinsIntegration.size(); i++) {
        gRecoilJetPtIntegralPhi[i] = (TGraphAsymmErrors*) dRecoilJets->Get(Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBinsIntegration[i].first, recoilJetPtBinsIntegration[i].second));
    }
    Plotting1D PRecoilJetPtPhi[recoilJetPtBinsIntegration.size()];
    for (int i = 0; i < recoilJetPtBinsIntegration.size(); i++) {
        TString jetPtLabel = Form("%.0f < #it{p}_{T}^{jet} < %.0f GeV/#it{c}", recoilJetPtBinsIntegration[i].first, recoilJetPtBinsIntegration[i].second);
        PRecoilJetPtPhi[i].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        PRecoilJetPtPhi[i].New(gRecoilJetPtIntegralPhi[i]);
        PRecoilJetPtPhi[i].SetAxisLabel("#Delta#phi", "#frac{1}{N_{trig}} #frac{N}{d#Delta#phi}",1.,1.4); // TODO: check axis labels
        PRecoilJetPtPhi[i].NewLatex(0.2, 0.9, Form("ALICE work in progress;Pb-Pb LHC23; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger; %s", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), jetPtLabel.Data()));
        PRecoilJetPtPhi[i].Plot(Form("%s/RecoilJetPtPhi_%d.%s", savePath.Data(), i, suffix.Data()), false, false);
    }


}

template <typename T>
void drawABCDRegions(T& plot)
{
  plot.NewLine(regionA.minM02, -20, regionA.minM02, regionA.maxIso, 2, 3, kBlack);
  plot.NewLine(regionA.maxM02, -20, regionA.maxM02, regionA.maxIso, 2, 3, kBlack);
  plot.NewLine(regionA.minM02, regionA.maxIso, regionA.maxM02, regionA.maxIso, 2, 3, kBlack);
  // region B
  plot.NewLine(regionB.minM02, -20, regionB.minM02, regionB.maxIso, 2, 3, kBlack);
  plot.NewLine(regionB.maxM02, -20, regionB.maxM02, regionB.maxIso, 2, 3, kBlack);
  plot.NewLine(regionB.minM02, regionB.maxIso, regionB.maxM02, regionB.maxIso, 2, 3, kBlack);
  // region C
  plot.NewLine(regionC.minM02, regionC.minIso, regionC.minM02, 80., 2, 3, kBlack);
  plot.NewLine(regionC.maxM02, regionC.minIso, regionC.maxM02, 80., 2, 3, kBlack);
  plot.NewLine(regionC.minM02, regionC.minIso, regionC.maxM02, regionC.minIso, 2, 3, kBlack);
  // region D
  plot.NewLine(regionD.minM02, regionD.minIso, regionD.minM02, regionD.maxIso, 2, 3, kBlack);
  plot.NewLine(regionD.maxM02, regionD.minIso, regionD.maxM02, regionD.maxIso, 2, 3, kBlack);
  plot.NewLine(regionD.minM02, regionD.minIso, regionD.maxM02, regionD.minIso, 2, 3, kBlack);


}

void plotRho(TDirectory *dRho, TString triggerParticle = "IsoPhoton")
{
    TString savePath = Form("%s/%sTrigger", outputDir.Data(), triggerParticle.Data());
   createDirectory(savePath.Data());

   TString triggerLabel = "";
   if (triggerParticle == "IsoPhoton") triggerLabel = "#gamma_{iso}";
   else if (triggerParticle == "MergedPi0") triggerLabel = "#pi^{0}_{merged}";

    // get hRhoSignal
    TH1F* hRhoSignal = (TH1F*) dRho->Get("hRhoSignal");
    TH1F* hRhoReference = (TH1F*) dRho->Get("hRhoReference");
    TH1F* hRhoSignalReferenceRatio = (TH1F*) dRho->Get("hRhoRatio_-0.0");
    TF1* fRhoFit = (TF1*) dRho->Get("fRhoRatiosFit_-0.0");



    double maxY = hRhoReference->GetMaximum();
    if (hRhoSignal->GetMaximum() > maxY) maxY = hRhoSignal->GetMaximum();

    PlottingRatio PRho;
    PRho.SetMargins(0.1, 0.2, 0.025, 0.025, 2000, 1500, 0.33);
    PRho.NewRatioLine(0, 1., 200., 1., 2, 2, kGray + 2);

    PRho.New(hRhoSignal, Form("sig %s TT{%.0f,%.0f}", triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second), 20, 2, kRed, "p");
    PRho.New(hRhoReference, Form("ref %s TT{%.0f,%.0f}", triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second), 4, 2, kBlack, "p");
    PRho.NewRatio(hRhoSignalReferenceRatio, "", 20, 2, kBlack, "p");
    PRho.SetLegendR(0.65, 0.9, 0.3, 0.45, true);
    PRho.NewRatio(fRhoFit, Form("pol1: %.3f #rho + %.2f", fRhoFit->GetParameter(1), fRhoFit->GetParameter(0)), 20, 2, fancyColors[0], "l");

    PRho.SetLegend(0.5, 0.8, 0.6, 0.7, true);
    PRho.SetAxisRange(0, 200, 1E-4, maxY*1.8, 0., 2.6);
    PRho.NewLatex(0.9, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger", centralityRange.first, centralityRange.second, jetRadius, triggerParticle.Data()));
    PRho.SetAxisLabel("#rho (GeV/#it{c})", "#frac{1}{N_{trig}} #frac{dN}{d#rho}", "signal / reference", 1, 1.4);
    PRho.Plot(Form("%s/Rho.%s", savePath.Data(), suffix.Data()), false, true);
}

void plotPurity(TDirectory *dPhotons, TString triggerParticle = "IsoPhoton")
{
    TString savePath = Form("%s/%sTrigger", outputDir.Data(), triggerParticle.Data());
   createDirectory(savePath.Data());

      TString triggerLabel = "";
   if (triggerParticle == "IsoPhoton") triggerLabel = "#gamma_{iso}";
   else if (triggerParticle == "MergedPi0") triggerLabel = "#pi^{0}_{merged}";

    TH2F* hM02vsIsoSignal = (TH2F*) dPhotons->Get("hM02vsIsoSignal");
    TH2F* hM02vsIsoReference = (TH2F*) dPhotons->Get("hM02vsIsoReference");

    Plotting2D PM02vsIsoSignal;
    PM02vsIsoSignal.New(hM02vsIsoSignal);
    PM02vsIsoSignal.SetAxisLabel("#it{#sigma}^{2}_{long}", "p_{T}^{iso} (GeV/#it{c})");
    PM02vsIsoSignal.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second));
    PM02vsIsoSignal.SetMargins(0.12, 0.12, 0.025, 0.17);
    drawABCDRegions(PM02vsIsoSignal);
    PM02vsIsoSignal.Plot(Form("%s/M02vsIsoSignal.%s", savePath.Data(), suffix.Data()), false, false, true);
    
    Plotting2D PM02vsIsoReference;
    PM02vsIsoReference.New(hM02vsIsoReference);
    PM02vsIsoReference.SetAxisLabel("#it{#sigma}^{2}_{long}", "p_{T}^{iso} (GeV/#it{c})");
    PM02vsIsoReference.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; ch. jets, anti-k_{T}, R=%.1f; %s trigger TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second));
    PM02vsIsoReference.SetMargins(0.12, 0.12, 0.025, 0.17);
    drawABCDRegions(PM02vsIsoReference);
    PM02vsIsoReference.Plot(Form("%s/M02vsIsoReference.%s", savePath.Data(), suffix.Data()), false, false, true);


}

void plotExclGammaJet(TString AnalysisDirectory, bool isDebugRun = false)
{
  TString fancyColorHex[8] = {"#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff"};
  for (int i = 0; i < 8; i++){
      fancyColors[i] = TColor::GetColor(fancyColorHex[i]);
  }

  INFO(Form("Plotting excl. gamma-jet analysis for %s", AnalysisDirectory.Data()));
  GlobalOptions optns(AnalysisDirectory, isDebugRun);
  extractConfigurationInformation(AnalysisDirectory, optns);

  outputDir = Form("%s/PlottingExclGammaJet", optns.analysisDirPath.Data());
  createDirectory(outputDir);

  TString inputFilePath = Form("%s/ExclGammaJet.root", AnalysisDirectory.Data());

  TFile *fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  // get all needed directories
  TDirectory *dRho_Photons = fIn->GetDirectory("Rho_Photons");
  TDirectory *dRho_MergedPi0s = fIn->GetDirectory("Rho_MergedPi0s");
  TDirectory *dRecoilJet_Photons = fIn->GetDirectory("RecoilJet_Photons");
  TDirectory *dRecoilJet_MergedPi0s = fIn->GetDirectory("RecoilJet_MergedPi0s");
  TDirectory *dPurity_Photons = fIn->GetDirectory("Purity_Photons");

  plotRecoilJets(dRecoilJet_Photons, "IsoPhoton",optns);
  plotRecoilJets(dRecoilJet_MergedPi0s, "MergedPi0",optns);

  plotRho(dRho_Photons, "IsoPhoton");
  plotRho(dRho_MergedPi0s, "MergedPi0");
  plotPurity(dPurity_Photons, "IsoPhoton");
}
