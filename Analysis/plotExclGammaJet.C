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

void plotRecoilJets(TDirectory *dRecoilJets, TString triggerParticle = "IsoPhoton")
{
   // create output subdirectory
   TString savePath = Form("%s/%sTrigger", outputDir.Data(), triggerParticle.Data());
   createDirectory(savePath.Data());

   TString triggerLabel = "";
   if (triggerParticle == "IsoPhoton") triggerLabel = "#gamma_{iso}";
   else if (triggerParticle == "MergedPi0") triggerLabel = "#pi^{0}_{merged}";


   // signal trigger
   TH2F* hIsoGammaJetDeltaPhiJetPtSignal = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtSignal");
   Plotting2D PDeltaPhiJetPtSignal;
   PDeltaPhiJetPtSignal.New(hIsoGammaJetDeltaPhiJetPtSignal);
   PDeltaPhiJetPtSignal.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtSignal.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; signal %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), signalTriggerPt.first, signalTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtSignal.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtSignal.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignal.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtSignal.Plot(Form("%s/DeltaPhiJetPtSignal.%s", savePath.Data(), suffix.Data()), false, false, true);

   // reference trigger
   TH2F* hIsoGammaJetDeltaPhiJetPtReference = (TH2F*) dRecoilJets->Get("hIsoGammaJetDeltaPhiJetPtReference");
   Plotting2D PDeltaPhiJetPtReference;
   PDeltaPhiJetPtReference.New(hIsoGammaJetDeltaPhiJetPtReference);
   PDeltaPhiJetPtReference.SetAxisLabel("#Delta#phi", "#it{p}_{T,jet}^{reco} (GeV/#it{c})");
   PDeltaPhiJetPtReference.NewLatex(0.8, 0.9, Form("ALICE work in progress;Pb-Pb LHC23zzm; %.0f - %.0f%% centrality; Ch-particle jets, anti-k_{T}, R=%.1f; reference %s trigger, TT{%.0f,%.0f}", centralityRange.first, centralityRange.second, jetRadius, triggerLabel.Data(), referenceTriggerPt.first, referenceTriggerPt.second), StdTextSize, 0.05);
   PDeltaPhiJetPtReference.SetMargins(0.12, 0.12, 0.025, 0.17);
   PDeltaPhiJetPtReference.NewLine(TMath::Pi()-0.6, -100, TMath::Pi()-0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReference.NewLine(TMath::Pi()+0.6, -100, TMath::Pi()+0.6, 200, 2, 3, kBlack);
   PDeltaPhiJetPtReference.Plot(Form("%s/DeltaPhiJetPtReference.%s", savePath.Data(), suffix.Data()), false, false, true);

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

  plotRecoilJets(dRecoilJet_Photons, "IsoPhoton");
  plotRecoilJets(dRecoilJet_MergedPi0s, "MergedPi0");

  plotRho(dRho_Photons, "IsoPhoton");
  plotRho(dRho_MergedPi0s, "MergedPi0");
  plotPurity(dPurity_Photons, "IsoPhoton");
}
