#include "PlottingClass.h"

const char *suffix = "png";

void plotJets(TDirectory *dJets, TString JetCutString)
{
  ENTER

  const char *outputDir = "Output/Jets";

  // Plot response matrix
  TH2F *hDLJetPtVsPLJetPt = (TH2F *)dJets->Get("hDLJetPtVsPLJetPt");
  hDLJetPtVsPLJetPt->Sumw2();
  Plotting2D PResponseMatrix;
  PResponseMatrix.New(hDLJetPtVsPLJetPt);
  PResponseMatrix.SetAxisLabel("#bf{#it{p}_{T}^{DL} (GeV/#it{c})}", "#bf{#it{p}_{T}^{PL} (GeV/#it{c})}");
  PResponseMatrix.SetAxisRange(5, 200, 5, 200, 1E-10, 1E1);
  PResponseMatrix.Plot(Form("%s/ResponseMatrix.%s", outputDir, suffix), 1, 1, 1);

  // Plot response matrix projections
  const int NProjections = 5;
  double projectionRanges[NProjections][2] = {{9, 11}, {19, 21}, {49, 51}, {99, 101}, {149, 151}};
  int rebinFactors[NProjections] = {2, 2, 2, 2, 5};
  // const int NProjections = 10;
  // double projectionRanges[NProjections][2] = {{5, 7}, {9, 11}, {14, 16}, {19, 21}, {29, 31}, {49, 51}, {74, 76}, {99, 101}, {149, 151}, {198, 200}};
  TH1F *hDLJetPtVsPLJetPtProjection[NProjections];
  TGraphAsymmErrors *gRatioDLJetPtToPLJetPtProjection[NProjections];
  Plotting1D PResponseMatrixProjections;
  for (int iProj = 0; iProj < NProjections; iProj++)
  {
    int binProjFrom = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][0] + 0.001);
    int binProjTo = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][1] + 0.001);
    hDLJetPtVsPLJetPtProjection[iProj] = (TH1F *)hDLJetPtVsPLJetPt->ProjectionX(Form("hDLJetPtVsPLJetPtProjection_%d", iProj), binProjFrom, binProjTo);
    hDLJetPtVsPLJetPtProjection[iProj]->Rebin(rebinFactors[iProj]);
    gRatioDLJetPtToPLJetPtProjection[iProj] = new TGraphAsymmErrors();
    float PLJetpT = (projectionRanges[iProj][0] + projectionRanges[iProj][1]) / 2.;
    for (int iDLJetPtBin = 0; iDLJetPtBin < hDLJetPtVsPLJetPtProjection[iProj]->GetNbinsX(); iDLJetPtBin++)
    {
      float DLJetpT = hDLJetPtVsPLJetPtProjection[iProj]->GetBinCenter(iDLJetPtBin);
      float AddBinContent = hDLJetPtVsPLJetPtProjection[iProj]->GetBinContent(iDLJetPtBin);
      float AddBinError = hDLJetPtVsPLJetPtProjection[iProj]->GetBinError(iDLJetPtBin);
      gRatioDLJetPtToPLJetPtProjection[iProj]->AddPoint(DLJetpT / PLJetpT, AddBinContent);
      gRatioDLJetPtToPLJetPtProjection[iProj]->SetPointError(iDLJetPtBin, 0, 0, AddBinError, AddBinError);
    }
    gRatioDLJetPtToPLJetPtProjection[iProj]->Scale(1. / gRatioDLJetPtToPLJetPtProjection[iProj]->Integral(0, hDLJetPtVsPLJetPtProjection[iProj]->GetNbinsX()));
    PResponseMatrixProjections.New(gRatioDLJetPtToPLJetPtProjection[iProj], Form("#it{p}_{T}^{PL} = %.0f #pm 1 GeV/#it{c}", PLJetpT));
  }
  PResponseMatrixProjections.SetAxisRange(0, 2.5);
  PResponseMatrixProjections.SetLegend(0.5,0.9,0.5,0.9);
  PResponseMatrixProjections.SetAxisLabel("#bf{#it{p}_{T}^{DL}/#it{p}_{T}^{PL}} #hat{=} #bf{#it{p}_{T}^{rec}/#it{p}_{T}^{gen}}");
  PResponseMatrixProjections.Plot(Form("%s/ResponseMatrixProjections.%s", outputDir, suffix));

  EXIT
}

void plotJetQA(TDirectory *dJetQA, TString JetCutString)
{
  ENTER

  const char *outputDir = "Output/Jets";

  // Plot eta phi map of reconstructed jets (DL)
  TH2F *hDLJetEtaPhi = (TH2F *)dJetQA->Get("hJetEtaPhi");
  Plotting2D PDLEtaPhiMap;
  PDLEtaPhiMap.SetMargins(0.12, 0.1, 0.05, 0.125);
  PDLEtaPhiMap.New(hDLJetEtaPhi);
  PDLEtaPhiMap.AddEMCalOutline();
  PDLEtaPhiMap.AddPHOSOutline();
  PDLEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PDLEtaPhiMap.Plot(Form("%s/DLEtaPhiMap.%s", outputDir, suffix));

  EXIT
}

void plotIsoGammaQA(TDirectory *dIsoGammaQA, TString IsoGammaCutString, bool isMC)
{
  ENTER

  const char *outputDir = "Output/IsoGammas";

  // Plot response matrix
  Plotting1D PMinMassDiffToPi0;
  TH1F *hMinMassDiffToPi0Signal;
  TH1F *hMinMassDiffToPi0Background;
  TH1F *hMinMassDiffToPi0;

  if (isMC)
  {
    hMinMassDiffToPi0Signal = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0Signal")->Clone("hMinMassDiffToPi0Signal");
    hMinMassDiffToPi0Background = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0Background")->Clone("hMinMassDiffToPi0Background");
    float integralSignal = hMinMassDiffToPi0Signal->Integral(0, hMinMassDiffToPi0Signal->GetNbinsX());
    float integralBG = hMinMassDiffToPi0Background->Integral(0, hMinMassDiffToPi0Background->GetNbinsX());
    hMinMassDiffToPi0Signal->Scale(1. / integralSignal);
    hMinMassDiffToPi0Background->Scale(1. / integralBG);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Signal, "Signal");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Background, "Background");
  }
  else
  {
    hMinMassDiffToPi0 = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0);
  }
  PMinMassDiffToPi0.SetAxisLabel("#bf{#Delta#it{m}_{inv} (GeV/#it{c}^{2})}");
  PMinMassDiffToPi0.Plot(Form("%s/MinMassDiffToPi0.%s", outputDir, suffix));

  EXIT
}

void plotHistosFromTree(bool isMC, bool doQA, TString EventCutString = "", TString IsoGammaCutString = "", TString JetCutString = "", TString Pi0CutString = "")
{
  TString inputFilePath = "Output/CalorimeTree.root";
  TFile *fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL("File not found")
  TDirectory *dJets = (TDirectory *)fIn->Get("Jets");
  if (!dJets)
    FATAL("Dir Jets not found")

  plotJets(dJets, JetCutString);

  TDirectory *dJetQA = (TDirectory *)fIn->Get("JetQA");
  if (!dJetQA)
    FATAL("Dir JetQA not found")

  plotJetQA(dJetQA, JetCutString);

  TDirectory *dIsoGammaQA = (TDirectory *)fIn->Get("IsoGammaQA");
  if (!dIsoGammaQA)
    FATAL("Dir IsoGammaQA not found")

  plotIsoGammaQA(dIsoGammaQA, IsoGammaCutString, isMC);

  fIn->Close();
}