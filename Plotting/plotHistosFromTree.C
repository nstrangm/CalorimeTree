#include "PlottingClass.h"

const char* suffix = "png";

void plotJets(TDirectory* dJets, TString JetCutString)
{
  ENTER

  const char* outputDir = "Output/Jets";

  // Plot response matrix
  TH2F* hDLJetPtVsPLJetPt = (TH2F*)dJets->Get("hDLJetPtVsPLJetPt");
  Plotting2D PResponseMatrix;
  PResponseMatrix.New(hDLJetPtVsPLJetPt);
  PResponseMatrix.SetAxisLabel("#bf{#it{p}_{T}^{DL} (GeV/#it{c})}", "#bf{#it{p}_{T}^{PL} (GeV/#it{c})}");
  PResponseMatrix.Plot(Form("%s/ResponseMatrix.%s", outputDir, suffix));

  EXIT
}

void plotJetQA(TDirectory* dJetQA, TString JetCutString)
{
  ENTER

  const char* outputDir = "Output/Jets";

  // Plot eta phi map of reconstructed jets (DL)
  TH2F* hDLJetEtaPhi = (TH2F*)dJetQA->Get("hJetEtaPhi");
  Plotting2D PDLEtaPhiMap;
  PDLEtaPhiMap.SetMargins(0.12, 0.1, 0.05, 0.125);
  PDLEtaPhiMap.New(hDLJetEtaPhi);
  PDLEtaPhiMap.AddEMCalOutline();
  PDLEtaPhiMap.AddPHOSOutline();
  PDLEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PDLEtaPhiMap.Plot(Form("%s/DLEtaPhiMap.%s", outputDir, suffix));

  EXIT
}

void plotIsoGammaQA(TDirectory* dIsoGammaQA, TString IsoGammaCutString, bool isMC)
{
  ENTER

  const char* outputDir = "Output/IsoGammas";

  // Plot response matrix
  Plotting1D PMinMassDiffToPi0;
  TH1F* hMinMassDiffToPi0Signal;
  TH1F* hMinMassDiffToPi0Background;
  TH1F* hMinMassDiffToPi0;

  if (isMC) {
    hMinMassDiffToPi0Signal = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0Signal")->Clone("hMinMassDiffToPi0Signal");
    hMinMassDiffToPi0Background = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0Background")->Clone("hMinMassDiffToPi0Background");
    float integralSignal = hMinMassDiffToPi0Signal->Integral(0, hMinMassDiffToPi0Signal->GetNbinsX());
    float integralBG = hMinMassDiffToPi0Background->Integral(0, hMinMassDiffToPi0Background->GetNbinsX());
    hMinMassDiffToPi0Signal->Scale(1. / integralSignal);
    hMinMassDiffToPi0Background->Scale(1. / integralBG);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Signal, "Signal");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Background, "Background");
  } else {
    hMinMassDiffToPi0 = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0);
  }
  PMinMassDiffToPi0.SetAxisLabel("#bf{#Delta#it{m}_{inv} (GeV/#it{c}^{2})}");
  PMinMassDiffToPi0.Plot(Form("%s/MinMassDiffToPi0.%s", outputDir, suffix));

  EXIT
}

void plotHistosFromTree(bool isMC, bool doQA, TString EventCutString = "", TString IsoGammaCutString = "", TString JetCutString = "", TString Pi0CutString = "")
{
  TString inputFilePath = "Output/CalorimeTree.root";
  TFile* fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL("File not found")
  TDirectory* dJets = (TDirectory*)fIn->Get("Jets");
  if (!dJets)
    FATAL("Dir Jets not found")

  plotJets(dJets, JetCutString);

  TDirectory* dJetQA = (TDirectory*)fIn->Get("JetQA");
  if (!dJetQA)
    FATAL("Dir JetQA not found")

  plotJetQA(dJetQA, JetCutString);

  TDirectory* dIsoGammaQA = (TDirectory*)fIn->Get("IsoGammaQA");
  if (!dIsoGammaQA)
    FATAL("Dir IsoGammaQA not found")

  plotIsoGammaQA(dIsoGammaQA, IsoGammaCutString, isMC);

  fIn->Close();
}