#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TKey.h"
#include "Logging.h"
#include "TFile.h"
#include "Utilities.h"
#include "PhysicsObjects.h"


TDirectory *DefineEventHistograms(TFile *f, GlobalOptions optns)
{
  TDirectory *dir = f->mkdir("Events");
  dir->cd();
  TH1F *hNEvents = new TH1F("hNEvents", "hNEvents", 5, -0.5, 4.5);
  hNEvents->GetXaxis()->SetBinLabel(1, "All");
  hNEvents->GetXaxis()->SetBinLabel(2, "Triggered");
  hNEvents->GetXaxis()->SetBinLabel(3, "Z-Vertex ok");
  hNEvents->GetXaxis()->SetBinLabel(4, "Quality ok");
  hNEvents->GetXaxis()->SetBinLabel(5, "Selected");

  return dir;
}

TDirectory *DefineEventQAHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doQA)
    return nullptr;
  TDirectory *dir = f->mkdir("EventQA");
  dir->cd();
  TH1F *hZVtx = new TH1F("hZVtx", "hZVtx; Vtx_{z} (cm)", 100, -20, 20);

  const int nWeightBins = 100;
  Double_t weightBinning[nWeightBins + 1];
  setExpArray(weightBinning, nWeightBins + 1, 1E-9, 1E-2);
  TH1F *hWeights = new TH1F("hWeights", "hWeights;Weight", nWeightBins, weightBinning);

  TH1F *hNEventswoWeights = new TH1F("hNEventswoWeights", "hNEventswoWeights", 5, -0.5, 4.5);
  hNEventswoWeights->GetXaxis()->SetBinLabel(1, "All");
  hNEventswoWeights->GetXaxis()->SetBinLabel(2, "Triggered");
  hNEventswoWeights->GetXaxis()->SetBinLabel(3, "Z-Vertex ok");
  hNEventswoWeights->GetXaxis()->SetBinLabel(4, "Quality ok");
  hNEventswoWeights->GetXaxis()->SetBinLabel(5, "Selected");

  return dir;
}

TDirectory *DefineIsoGammaHistograms(TFile *f, string dirname, GlobalOptions optns)
{
  if (!optns.doIsoGamma)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  TH1F *hIsoGammaPt = new TH1F("hIsoGammaPt", "hPt", 1000, 0., 100.);
  //dir->Add(hIsoGammaPt);
  TH1F *hIsoGammaE = new TH1F("hIsoGammaE", "hE", 1000, 0., 100.);
  //dir->Add(hIsoGammaE);

  if(optns.isMC){
    TH1F *hGammaGenPt = new TH1F("hGammaGenPt", "hPt", 1000, 0., 100.);
    TH1F *hGammaGenEAcceptanceCut = new TH1F("hGammaGenEAcceptanceCut", "hEAcceptanceCut", 1000, 0., 100.);
    TH1F *hGammaGenE = new TH1F("hGammaGenE", "hE", 1000, 0., 100.);
    TH1F *hGammaGenPtAcceptanceCut = new TH1F("hGammaGenPtAcceptanceCut", "hPtAcceptanceCut", 1000, 0., 100.);
    TH1F *hIsoGammaPtSignal = new TH1F("hIsoGammaPtSignal", "hPtSignal", 1000, 0., 100.);
    TH1F *hIsoGammaESignal = new TH1F("hIsoGammaESignal", "hESignal", 1000, 0., 100.);
    TH1F *hIsoGammaPtBackground = new TH1F("hIsoGammaPtBackground", "hPtBackground", 1000, 0., 100.);
    TH1F *hIsoGammaEBackground = new TH1F("hIsoGammaEBackground", "hEBackground", 1000, 0., 100.);
  }

  int nBins[3] = {1000, 200, 2000};
  double xmin[3] = {-20, 0, 0};
  double xmax[3] = {80, 2, 200};
  THnSparseF *hIsoGammaIsovsM02vsPt = new THnSparseF("hIsoGammaIsovsM02vsPt", "hIsoGammaIsovsM02vsPt", 3, nBins, xmin, xmax);
  hIsoGammaIsovsM02vsPt->Sumw2();
  dir->Add(hIsoGammaIsovsM02vsPt);

  if(optns.isMC){
    //no-cuts signal
    THnSparseF *hIsoGammaIsovsM02vsPtSignal = new THnSparseF("hIsoGammaIsovsM02vsPtSignal", "hIsoGammaIsovsM02vsPtSignal", 3, nBins, xmin, xmax);
    hIsoGammaIsovsM02vsPtSignal->Sumw2();
    dir->Add(hIsoGammaIsovsM02vsPtSignal);
    //no-cuts background
    THnSparseF *hIsoGammaIsovsM02vsPtBackground = new THnSparseF("hIsoGammaIsovsM02vsPtBackground", "hIsoGammaIsovsM02vsPtBackground", 3, nBins, xmin, xmax);
    hIsoGammaIsovsM02vsPtBackground->Sumw2();
    dir->Add(hIsoGammaIsovsM02vsPtBackground);
  }

  return dir;
}

TDirectory *DefineIsoGammaQAHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doQA)
    return nullptr;

  const float M02Range[2] = {0, 5};
  const float M20Range[2] = {0, 5};
  const float pxyzRange[2] = {-50, 50};
  const float pTRange[2] = {0, 50};

  // QA plots for all IsoGammas
  TDirectory *dir = f->mkdir("IsoGammaQA");
  dir->cd();

  TH1F *hNIsoGamma = new TH1F("hNIsoGamma", "hNIsoGamma", 21, -0.5, 20.5);

  TH2F *hpTSpectraAfterSubsequentCuts = new TH2F("hpTSpectraAfterSubsequentCuts", "hpTSpectraAfterSubsequentCuts", 10, -0.5, 9.5, 100, pTRange[0], pTRange[1]);
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(1, "All");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(2, "Accepted");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(3, "E_{min}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(4, "N_{cells}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(5, "NLM");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(6, "Dist BC");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(7, "Track matching");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(8, "Exotic");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(9, "M_{02}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(10, "#it{p}_{T}^{iso}");

  TH2F *hpTSpectrumLossFromIndividualCuts = new TH2F("hpTSpectrumLossFromIndividualCuts", "hpTSpectrumLossFromIndividualCuts", 10, -0.5, 8.5, 100, pTRange[0], pTRange[1]);
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(1, "All");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(2, "Acceptance");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(3, "E_{min}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(4, "N_{cells}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(5, "NLM");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(6, "Dist BC");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(7, "Track matching");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(8, "Exotic");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(9, "M_{02}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(10, "#it{p}_{T}^{iso}");

  TH1F *hIsoGammaIsoCharged = new TH1F("hIsoGammaIsoCharged", "hIsoGammaIsoCharged", 100, 0, 200);
  TH1F *hIsoGammaIsoChargedCorrected = new TH1F("hIsoGammaIsoChargedCorrected", "hIsoGammaIsoChargedCorrected", 100, -10, 30);
  TH1F *hIsoGammaE = new TH1F("hIsoGammaE", "hIsoGammaE", 100, 0, 156);
  TH1F *hIsoGammaM02 = new TH1F("hIsoGammaM02", "hIsoGammaM02", 100, M02Range[0], M02Range[1]);
  TH1F *hIsoGammaM20 = new TH1F("hIsoGammaM20", "hIsoGammaM20", 100, M20Range[0], M20Range[1]);
  TH1F *hIsoGammaPx = new TH1F("hIsoGammaPx", "hPx", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hIsoGammaPy = new TH1F("hIsoGammaPy", "hPy", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hIsoGammaPz = new TH1F("hIsoGammaPz", "hPz", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hIsoGammaNLM = new TH1F("hIsoGammaNLM", "hIsoGammaNLM", 6, 0, 6);
  TH2F *hIsoGammaEtaPhi = new TH2F("hIsoGammaEtaPhi", "hEtaPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());
  TH2F *hIsoGammaM02pT = new TH2F("hIsoGammaM02pT", "hIsoGammaM02pT", 100, M02Range[0], M02Range[1], 100, 0, 50);
  TH2F *hIsoGammaEBeforeAfterNL = new TH2F("hIsoGammaEBeforeAfterNL", "hIsoGammaEBeforeAfterNL;#bf{E_{cls}^{before} (GeV)};#bf{E_{cls}^{after}/E_{cls}^{before}}", 2000, 0, 200, 200, 0.8, 1.2);
  TH2F *hGGMinvDist = new TH2F("hGGMinvDist", "hGGMinvDist", 100, 0, 1, 500, 0, 50);
  TH2F *hIsoGammadEtadphi = new TH2F("hIsoGammadEtadphi","hIsoGammadEtadphi",100,-0.1,0.1,100,-0.1,0.1);
  TH2F *hIsoGammadEtadphiCut = new TH2F("hIsoGammadEtadphiCut","hIsoGammadEtadphiCut",100,-0.1,0.1,100,-0.1,0.1);

  if (optns.isMC)
  { // ------------------> MC IsoGammaQA histograms
    TH1F *hMinMassDiffToPi0Signal = new TH1F("hMinMassDiffToPi0Signal", "hMinMassDiffToPi0Signal", 50, 0., 0.1);
    TH1F *hMinMassDiffToPi0Background = new TH1F("hMinMassDiffToPi0Background", "hMinMassDiffToPi0Background", 50, 0., 0.1);

    TH1F *hTrueIsoGammaMCTag = new TH1F("hTrueIsoGammaMCTag", "hMCTag", 20, 0.5, 20.5);
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(1, "Prompt #gamma");
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(2, "Frag #gamma");
    //Generated photons QA:
    TH1F *hGammaGenIsoCharged = new TH1F("hGammaGenIsoCharged","hGammaGenIsoCharged",100,0,200);
    TH1F *hGammaGenBckPerb = new TH1F("hGammaGenBckPerb","hGammaGenBckPerb",100,0,200);

    // Signal IsoGammas
    TH1F *hNIsoGammaSignal = new TH1F("hNIsoGammaSignal", "hNIsoGammaSignal", 21, -0.5, 20.5);
    TH1F *hIsoGammaIsoChargedSignal = new TH1F("hIsoGammaIsoChargedSignal", "hIsoGammaIsoChargedSignal", 100, 0, 200);
    TH1F *hIsoGammaIsoChargedCorrectedSignal = new TH1F("hIsoGammaIsoChargedCorrectedSignal", "hIsoGammaIsoChargedCorrectedSignal", 100, -10, 30);
    TH1F *hIsoGammaESignal = new TH1F("hIsoGammaESignal", "hIsoGammaESignal", 100, 0, 156);
    TH1F *hIsoGammaM02Signal = new TH1F("hIsoGammaM02Signal", "hIsoGammaM02Signal", 100, M02Range[0], M02Range[1]);
    TH1F *hIsoGammaM20Signal = new TH1F("hIsoGammaM20Signal", "hIsoGammaM20Signal", 100, M20Range[0], M20Range[1]);
    TH1F *hIsoGammaPxSignal = new TH1F("hIsoGammaPxSignal", "hPxSignal", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaPySignal = new TH1F("hIsoGammaPySignal", "hPySignal", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaPzSignal = new TH1F("hIsoGammaPzSignal", "hPzSignal", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaNLMSignal = new TH1F("hIsoGammaNLMSignal", "hIsoGammaNLMSignal", 6, 0, 6);
    TH2F *hIsoGammaEtaPhiSignal = new TH2F("hIsoGammaEtaPhiSignal", "hEtaPhiSignal", 100, -1., 1., 100, 0., 2 * TMath::Pi());
    TH2F *hIsoGammaM02pTSignal = new TH2F("hIsoGammaM02pTSignal", "hIsoGammaM02pTSignal", 100, M02Range[0], M02Range[1], 100, 0, 50);
    TH2F *hIsoGammaEBeforeAfterNLSignal = new TH2F("hIsoGammaEBeforeAfterNLSignal", "hIsoGammaEBeforeAfterNLSignal;#bf{E_{cls}^{before} (GeV)};#bf{E_{cls}^{after}/E_{cls}^{before}}", 2000, 0, 200, 200, 0.8, 1.2);
    TH2F *hGGMinvDistSignal = new TH2F("hGGMinvDistSignal", "hGGMinvDistSignal", 100, 0, 1, 500, 0, 50);

    // NonSignal IsGammas
    TH1F *hNIsoGammaBackground = new TH1F("hNIsoGammaBackground", "hNIsoGammaBackground", 21, -0.5, 20.5);
    TH1F *hIsoGammaIsoChargedBackground = new TH1F("hIsoGammaIsoChargedBackground", "hIsoGammaIsoChargedBackground", 100, 0, 200);
    TH1F *hIsoGammaIsoChargedCorrectedBackground = new TH1F("hIsoGammaIsoChargedCorrectedBackground", "hIsoGammaIsoChargedCorrectedBackground", 100, -10, 30);
    TH1F *hIsoGammaEBackground = new TH1F("hIsoGammaEBackground", "hIsoGammaEBackground", 100, 0, 156);
    TH1F *hIsoGammaM02Background = new TH1F("hIsoGammaM02Background", "hIsoGammaM02Background", 100, M02Range[0], M02Range[1]);
    TH1F *hIsoGammaM20Background = new TH1F("hIsoGammaM20Background", "hIsoGammaM20Background", 100, M20Range[0], M20Range[1]);
    TH1F *hIsoGammaPxBackground = new TH1F("hIsoGammaPxBackground", "hPxBackground", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaPyBackground = new TH1F("hIsoGammaPyBackground", "hPyBackground", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaPzBackground = new TH1F("hIsoGammaPzBackground", "hPzBackground", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hIsoGammaNLMBackground = new TH1F("hIsoGammaNLMBackground", "hIsoGammaNLMBackground", 6, 0, 6);
    TH2F *hIsoGammaEtaPhiBackground = new TH2F("hIsoGammaEtaPhiBackground", "hEtaPhiBackground", 100, -1., 1., 100, 0., 2 * TMath::Pi());
    TH2F *hIsoGammaM02pTBackground = new TH2F("hIsoGammaM02pTBackground", "hIsoGammaM02pTBackground", 100, M02Range[0], M02Range[1], 100, 0, 50);
    TH2F *hIsoGammaEBeforeAfterNLBackground = new TH2F("hIsoGammaEBeforeAfterNLBackground", "hIsoGammaEBeforeAfterNLBackground;#bf{E_{cls}^{before} (GeV)};#bf{E_{cls}^{after}/E_{cls}^{before}}", 2000, 0, 200, 200, 0.8, 1.2);
    TH2F *hGGMinvDistBackground = new TH2F("hGGMinvDistBackground", "hGGMinvDistBackground", 100, 0, 1, 500, 0, 50);

    // All Gammas
    TH2F *hIsoGammaM02pTAllGamma = new TH2F("hIsoGammaM02pTAllGamma", "hM02pTAllGamma", 100, M02Range[0], M02Range[1], 100, 0, 50);
    // Pion-decay Gammas
    TH2F *hIsoGammaM02pTPionDecayGamma = new TH2F("hIsoGammaM02pTPionDecayGamma", "hM02pTPionDecayGamma", 100, M02Range[0], M02Range[1], 100, 0, 50);
    // Eta-decay Gammas
    TH2F *hIsoGammaM02pTEtaDecayGamma = new TH2F("hIsoGammaM02pTEtaDecayGamma", "hM02pTEtaDecayGamma", 100, M02Range[0], M02Range[1], 100, 0, 50);
    // Merged-pion-Gammas
    TH2F *hIsoGammaM02pTMergedPionGamma = new TH2F("hIsoGammaM02pTMergedPionGamma", "hM02pTMergedPionGamma", 100, M02Range[0], M02Range[1], 100, 0, 50);
  }
  else
  { // ------------------> Only in data isogamma histograms
    TH1F *hMinMassDiffToPi0 = new TH1F("hMinMassDiffToPi0", "hMinMassDiffToPi0", 50, 0., 0.1);
  }

  return dir;
}

TDirectory *DefineJetHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doJets)
    return nullptr;
  TDirectory *dir = f->mkdir("Jets");
  dir->cd();

  TH1F *hJetPt = new TH1F("hJetPt", "hPt", 1000, 0., 100.);

  // ------------------> MC Jet histograms
  if (optns.isMC)
  {
    TH1F *hPLJetPt = new TH1F("hPLJetPt", "hPt", 1000, 0., 100.);

    TH2F *hDLJetPtVsPLJetPt = new TH2F("hDLJetPtVsPLJetPt", "hPt;#it{p}_{T}^{DL Jet} (GeV/#it{c});#it{p}_{T}^{PL Jet} (GeV/#it{c})", 400, 0., 200., 400, 0., 200.);
    TH2F *hDLSubPLDivPLJetPtVsPLJetPt = new TH2F("hDLSubPLDivPLJetPtVsPLJetPt", "hPt;(#it{p}_{T}^{DL Jet}-#it{p}_{T}^{PL Jet})/#it{p}_{T}^{PL Jet} ;#it{p}_{T}^{PL Jet} (GeV/#it{c})", 400, -1, 3., 400, 0., 200.);
  }

  return dir;
}

TDirectory *DefineJetQAHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doJets || !optns.doQA)
    return nullptr;
  const float pxyzRange[2] = {-50, 50};
  TDirectory *dir = f->mkdir("JetQA");
  dir->cd();
  TH1F *hJetPx = new TH1F("hJetPx", "hPx", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetPy = new TH1F("hJetPy", "hPy", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetPz = new TH1F("hJetPz", "hPz", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetArea = new TH1F("hJetArea", "hArea", 100, 0., 1.);
  TH1F *hJetNch = new TH1F("hJetNch", "hNch", 51, -0.5, 50.5);
  TH1F *hJetNclus = new TH1F("hJetNclus", "hNclus", 51, -0.5, 50.5);
  TH2F *hJetEtaPhi = new TH2F("hJetEtaPhi", "hJetPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());

  // ------------------> MC Jet QA histograms
  if (optns.isMC)
  {
    TH1F *hPLJetPx = new TH1F("hPLJetPx", "hPx", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hPLJetPy = new TH1F("hPLJetPy", "hPy", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hPLJetPz = new TH1F("hPLJetPz", "hPz", 1000, pxyzRange[0], pxyzRange[1]);
    TH1F *hPLJetArea = new TH1F("hPLJetArea", "hArea", 100, 0., 1.);
    TH1F *hPLJetNPart = new TH1F("hPLJetNPart", "hNPart", 51, -0.5, 50.5);
    TH2F *hPLJetEtaPhi = new TH2F("hPLJetEtaPhi", "hPLJetPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());
    TH2F *hDLJetPVsPLJetP = new TH2F("hDLJetPVsPLJetP", "hPt;#it{p}_{T}^{DL Jet} (GeV/#it{c});#it{p}_{T}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
    TH2F *hDLJetPxVsPLJetPx = new TH2F("hDLJetPxVsPLJetPx", "hPx;#it{p}_{x}^{DL Jet} (GeV/#it{c});#it{p}_{x}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
    TH2F *hDLJetPyVsPLJetPy = new TH2F("hDLJetPyVsPLJetPy", "hPx;#it{p}_{y}^{DL Jet} (GeV/#it{c});#it{p}_{y}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
  }

  return dir;
}

TDirectory *DefineGammaJetHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doJets)
    return nullptr;
  TDirectory *dir = f->mkdir("GammaJetCorrelations");
  dir->cd();

  TH2F *hpTImbalancevsDeltaPhi = new TH2F("hpTImbalancevsDeltaPhi", "hpTImbalancevsDeltaPhi;#bf{#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}};#bf{#Delta#phi = |#phi_{#gamma}-#phi_{jet}|}", 100, 0., 2., 314, 0, TMath::Pi());

  return dir;
}

TDirectory *DefineGammaJetQAHistograms(TFile *f, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doJets || !optns.doQA)
    return nullptr;
  TDirectory *dir = f->mkdir("GammaJetCorrelationQA");
  dir->cd();

  TH2F *hNIsoGammaJets = new TH2F("hNIsoGammaJets", "hNIsoGammaJets;#bf{N_{#gamma}};#bf{N_{jet}}", 11, -0.5, 10.5, 11, -0.5, 10.5);

  TH1F *hpTImbalance = new TH1F("hpTImbalance", "hpTImbalance;#bf{#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}}", 100, 0., 2.);
  TH1F *hDeltaPhi = new TH1F("hDeltaPhi", "hDeltaPhi;#bf{#Delta#phi = |#phi_{#gamma}-#phi_{jet}|}", 314, 0, TMath::Pi());
  
  return dir;
}

template <typename T>
void fillHistograms(T obj, TDirectory *dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<IsoGamma>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hIsoGammaPt"))->Fill(obj.at(i).Pt(), eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaE"))->Fill(obj.at(i).E, eventWeight);
      ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPt"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
      if(optns.isMC){
        if(obj.at(i).isSignal()){
          ((TH1F *)dir->FindObject("hIsoGammaPtSignal"))->Fill(obj.at(i).Pt(), eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaESignal"))->Fill(obj.at(i).E, eventWeight);
          ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPtSignal"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }else{
          ((TH1F *)dir->FindObject("hIsoGammaPtBackground"))->Fill(obj.at(i).Pt(), eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaEBackground"))->Fill(obj.at(i).E, eventWeight);
          ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPtBackground"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<Jet>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, std::vector<PLJet>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hPLJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
      if (obj.at(i).ClosestDLJet != nullptr){
        ((TH2F *)dir->FindObject("hDLJetPtVsPLJetPt"))->Fill(obj.at(i).ClosestDLJet->Pt(), obj.at(i).Pt(), eventWeight);
        ((TH2F *)dir->FindObject("hDLSubPLDivPLJetPtVsPLJetPt"))->Fill((obj.at(i).ClosestDLJet->Pt()-obj.at(i).Pt())/obj.at(i).Pt(), obj.at(i).Pt(), eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<GammaJetPair>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH2F *)dir->FindObject("hpTImbalancevsDeltaPhi"))->Fill(obj.at(i).pTImbalance, obj.at(i).DPhi, eventWeight);
    }
  }
  if constexpr (std::is_same<T, Event>::value)
  {
    ((TH1F *)dir->FindObject("hNEvents"))->Fill((int)0, eventWeight);
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(1, ((int)(obj.IsTriggered)) * eventWeight);
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(2, ((int)(obj.ZVtxOK)) * eventWeight);
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(3, ((int)(obj.QualityOK)) * eventWeight);
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(4, ((int)(obj.Selected)) * eventWeight);
  }
}

template <typename T>
void fillGammaGenHistograms(T obj, TDirectory *dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<GammaGen>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hGammaGenE"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenPt"))->Fill(obj.at(i).Pt(), eventWeight);
    }
  }
}

template <typename T>
void fillGammaGenAcceptanceCutHistograms(T obj, TDirectory *dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<GammaGen>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hGammaGenEAcceptanceCut"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenPtAcceptanceCut"))->Fill(obj.at(i).Pt(), eventWeight);
    }
  }
}

template <typename T>
void fillGammaGenQAHistograms(T obj, TDirectory *dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<GammaGen>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hGammaGenIsoCharged"))->Fill(obj.at(i).IsoCharged,eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenBckPerb"))->Fill(obj.at(i).IsoCharged,eventWeight);
    }
  }
}

template <typename T>
void fillQAHistograms(T obj, TDirectory *dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<IsoGamma>>::value)
  {
    ((TH1F *)dir->FindObject("hNIsoGamma"))->Fill((int)obj.size(), eventWeight);
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      //
      ((TH1F *)dir->FindObject("hIsoGammaIsoCharged"))->Fill(obj.at(i).IsoCharged, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaIsoChargedCorrected"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaPy"))->Fill(obj.at(i).Py, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaPz"))->Fill(obj.at(i).Pz, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaE"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaM02"))->Fill(obj.at(i).M02, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaM20"))->Fill(obj.at(i).M20, eventWeight);
      ((TH1F *)dir->FindObject("hIsoGammaNLM"))->Fill(obj.at(i).NLM, eventWeight);
      ((TH2F *)dir->FindObject("hIsoGammaEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      ((TH2F *)dir->FindObject("hIsoGammaM02pT"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
      ((TH2F *)dir->FindObject("hIsoGammaEBeforeAfterNL"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
      ((TH2F *)dir->FindObject("hIsoGammadEtadphi"))->Fill(obj.at(i).MatchedTrack.dPhi, obj.at(i).MatchedTrack.dEta, eventWeight);
      //TH2F *hIsoGammadEtadphi = new TH2F("hIsoGammadEtadphiCut","hIsoGammadEtadphiCut",100,-1,1,100,-1,1);
      if (optns.isMC)
      {
        ((TH2F *)dir->FindObject("hTrueIsoGammaMCTag"))->Fill(obj.at(i).MCTag, eventWeight);
        if (obj.at(i).isSignal())
        {
          ((TH1F *)dir->FindObject("hIsoGammaIsoChargedSignal"))->Fill(obj.at(i).IsoCharged, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaIsoChargedCorrectedSignal"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
          ((TH1F *)dir->FindObject("hNIsoGammaSignal"))->Fill((int)obj.size(), eventWeight);
          ((TH1F *)dir->FindObject("hMinMassDiffToPi0Signal"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPxSignal"))->Fill(obj.at(i).Px, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPySignal"))->Fill(obj.at(i).Py, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPzSignal"))->Fill(obj.at(i).Pz, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaESignal"))->Fill(obj.at(i).E, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaM02Signal"))->Fill(obj.at(i).M02, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaM20Signal"))->Fill(obj.at(i).M20, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaNLMSignal"))->Fill(obj.at(i).NLM, eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaEtaPhiSignal"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaM02pTSignal"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaEBeforeAfterNLSignal"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
        }
        else
        {
          ((TH1F *)dir->FindObject("hIsoGammaIsoChargedBackground"))->Fill(obj.at(i).IsoCharged, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaIsoChargedCorrectedBackground"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
          ((TH1F *)dir->FindObject("hNIsoGammaBackground"))->Fill((int)obj.size(), eventWeight);
          ((TH1F *)dir->FindObject("hMinMassDiffToPi0Background"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPxBackground"))->Fill(obj.at(i).Px, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPyBackground"))->Fill(obj.at(i).Py, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaPzBackground"))->Fill(obj.at(i).Pz, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaEBackground"))->Fill(obj.at(i).E, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaM02Background"))->Fill(obj.at(i).M02, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaM20Background"))->Fill(obj.at(i).M20, eventWeight);
          ((TH1F *)dir->FindObject("hIsoGammaNLMBackground"))->Fill(obj.at(i).NLM, eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaEtaPhiBackground"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaM02pTBackground"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
          ((TH2F *)dir->FindObject("hIsoGammaEBeforeAfterNLBackground"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
        }
        if (obj.at(i).isPhoton())
        {
          ((TH2F *)dir->FindObject("hIsoGammaM02pTAllGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isGammaFromPion())
        {
          ((TH2F *)dir->FindObject("hIsoGammaM02pTPionDecayGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isGammaFromEta())
        {
          ((TH2F *)dir->FindObject("hIsoGammaM02pTEtaDecayGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isPion())
        {
          ((TH2F *)dir->FindObject("hIsoGammaM02pTMergedPionGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
      }
      else
      {
        ((TH1F *)dir->FindObject("hMinMassDiffToPi0"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<Pi0>>::value)
  {
    if (((TString)dir->GetName()) == "IsoGammaQA")
    {
      for (unsigned long i = 0; i < obj.size(); i++)
      {
        ((TH2F *)dir->FindObject("hGGMinvDist"))->Fill(obj.at(i).Mass, obj.at(i).Pt(), eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<Jet>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hJetPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH1F *)dir->FindObject("hJetPy"))->Fill(obj.at(i).Py, eventWeight);
      ((TH1F *)dir->FindObject("hJetPz"))->Fill(obj.at(i).Pz, eventWeight);
      ((TH1F *)dir->FindObject("hJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F *)dir->FindObject("hJetNch"))->Fill(obj.at(i).Nch, eventWeight);
      ((TH1F *)dir->FindObject("hJetNclus"))->Fill(obj.at(i).Nclus, eventWeight);
      ((TH2F *)dir->FindObject("hJetEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, std::vector<PLJet>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hPLJetPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH1F *)dir->FindObject("hPLJetPy"))->Fill(obj.at(i).Py, eventWeight);
      ((TH1F *)dir->FindObject("hPLJetPz"))->Fill(obj.at(i).Pz, eventWeight);
      ((TH1F *)dir->FindObject("hPLJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F *)dir->FindObject("hPLJetNPart"))->Fill(obj.at(i).NPart, eventWeight);
      ((TH2F *)dir->FindObject("hPLJetEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      if (obj.at(i).ClosestDLJet != nullptr)
      {
        ((TH2F *)dir->FindObject("hDLJetPVsPLJetP"))->Fill(obj.at(i).ClosestDLJet->P(), obj.at(i).P(), eventWeight);
        ((TH2F *)dir->FindObject("hDLJetPxVsPLJetPx"))->Fill(fabs(obj.at(i).ClosestDLJet->Px), fabs(obj.at(i).Px), eventWeight);
        ((TH2F *)dir->FindObject("hDLJetPyVsPLJetPy"))->Fill(fabs(obj.at(i).ClosestDLJet->Py), fabs(obj.at(i).Py), eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<GammaJetPair>>::value)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)dir->FindObject("hpTImbalance"))->Fill(obj.at(i).pTImbalance, eventWeight);
      ((TH1F *)dir->FindObject("hDeltaPhi"))->Fill(obj.at(i).DPhi, eventWeight);
    }
  }
  if constexpr (std::is_same<T, Event>::value)
  {
    ((TH1F *)dir->FindObject("hZVtx"))->Fill(obj.ZVtx, eventWeight);
    if (dir->FindObject("hWeights"))
    {
      ((TH1F *)dir->FindObject("hWeights"))->Fill(obj.weight);
    }
    ((TH1F *)dir->FindObject("hNEventswoWeights"))->Fill((int)0, 1);
    ((TH1F *)dir->FindObject("hNEventswoWeights"))->Fill(1, (int)(obj.IsTriggered));
    ((TH1F *)dir->FindObject("hNEventswoWeights"))->Fill(2, (int)(obj.ZVtxOK));
    ((TH1F *)dir->FindObject("hNEventswoWeights"))->Fill(3, (int)(obj.QualityOK));
    ((TH1F *)dir->FindObject("hNEventswoWeights"))->Fill(4, (int)(obj.Selected));
  }
}
