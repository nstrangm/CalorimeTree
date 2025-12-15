#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TKey.h"
#include "Logging.h"
#include "TFile.h"
#include "Utilities.h"
#include "PhysicsObjects.h"
#include "Cuts.h"
#include <bitset>

TDirectory *DefineEventHistograms(TFile *f, GlobalOptions optns)
{
  // TODO: add software triggers
  TDirectory *dir = f->mkdir("Events");
  dir->cd();
  TH1F *hNEvents = new TH1F("hNEvents", "hNEvents", 6, -0.5, 5.5);
  hNEvents->GetXaxis()->SetBinLabel(1, "All");
  hNEvents->GetXaxis()->SetBinLabel(2, "kTVXinEMC");
  hNEvents->GetXaxis()->SetBinLabel(3, "kEMC7");
  hNEvents->GetXaxis()->SetBinLabel(4, "kDMC7");
  hNEvents->GetXaxis()->SetBinLabel(5, "Fulfill Occupancy");
  hNEvents->GetXaxis()->SetBinLabel(6, "Selected");

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

  TH1F *hNEventswoWeights = new TH1F("hNEventswoWeights", "hNEventswoWeights", 6, -0.5, 5.5);
  hNEventswoWeights->GetXaxis()->SetBinLabel(1, "All");
  hNEventswoWeights->GetXaxis()->SetBinLabel(2, "kTVXinEMC");
  hNEventswoWeights->GetXaxis()->SetBinLabel(3, "kEMC7");
  hNEventswoWeights->GetXaxis()->SetBinLabel(4, "kDMC7");
  hNEventswoWeights->GetXaxis()->SetBinLabel(5, "Fulfill Occupancy");
  hNEventswoWeights->GetXaxis()->SetBinLabel(6, "Selected");

  THnSparseF* hCentralityRhoMultiplicity = new THnSparseF("hCentralityRhoMultiplicity", "hCentralityRhoMultiplicity", 3, new int[3]{100, 200, 500}, new double[3]{0, 0, 0}, new double[3]{100, 100, 20000});
  hCentralityRhoMultiplicity->GetAxis(0)->SetTitle("Centrality");
  hCentralityRhoMultiplicity->GetAxis(1)->SetTitle("Rho (GeV/c)");
  hCentralityRhoMultiplicity->GetAxis(2)->SetTitle("Multiplicity");
  hCentralityRhoMultiplicity->Sumw2();
  dir->Add(hCentralityRhoMultiplicity);

  TH1F* hOccupancy = new TH1F("hOccupancy", "hOccupancy", 500, 0, 20000);
  hOccupancy->GetXaxis()->SetTitle("Occupancy");
  dir->Add(hOccupancy);

  TH2F* hOccupancyCentrality = new TH2F("hOccupancyCentrality", "hOccupancyCentrality", 100, 0, 100, 500, 0, 20000);
  hOccupancyCentrality->GetXaxis()->SetTitle("Centrality");
  hOccupancyCentrality->GetYaxis()->SetTitle("Occupancy");
  dir->Add(hOccupancyCentrality);


  return dir;
}

TDirectory *DefineIsoGammaHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doIsoGamma)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  TH1F *hIsoGammaPt = new TH1F("hIsoGammaPt", "hPt", 1000, 0., 100.);
  TH1F *hIsoGammaE = new TH1F("hIsoGammaE", "hE", 1000, 0., 100.);

  if (optns.isMC)
  {
    TH1F *hGammaGenE = new TH1F("hGammaGenE", "hE", 1000, 0., 100.);
    TH1F *hGammaGenPt = new TH1F("hGammaGenPt", "hPt", 1000, 0., 100.);
    TH1F *hGammaGenESignal = new TH1F("hGammaGenESignal", "hESignal", 1000, 0., 100.);
    TH1F *hGammaGenPtSignal = new TH1F("hGammaGenPtSignal", "hPtSignal", 1000, 0., 100.);
    TH1F *hGammaGenEBackground = new TH1F("hGammaGenEBackground", "hEBackground", 1000, 0., 100.);
    TH1F *hGammaGenPtBackground = new TH1F("hGammaGenPtBackground", "hPtBackground", 1000, 0., 100.);

    TH1F *hGammaGenEAcceptanceCut = new TH1F("hGammaGenEAcceptanceCut", "hEAcceptanceCut", 1000, 0., 100.);
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

  if (optns.isMC)
  {
    // no-cuts signal
    THnSparseF *hIsoGammaIsovsM02vsPtSignal = new THnSparseF("hIsoGammaIsovsM02vsPtSignal", "hIsoGammaIsovsM02vsPtSignal", 3, nBins, xmin, xmax);
    hIsoGammaIsovsM02vsPtSignal->Sumw2();
    dir->Add(hIsoGammaIsovsM02vsPtSignal);
    // no-cuts background
    THnSparseF *hIsoGammaIsovsM02vsPtBackground = new THnSparseF("hIsoGammaIsovsM02vsPtBackground", "hIsoGammaIsovsM02vsPtBackground", 3, nBins, xmin, xmax);
    hIsoGammaIsovsM02vsPtBackground->Sumw2();
    dir->Add(hIsoGammaIsovsM02vsPtBackground);
  }

  return dir;
}

TDirectory *DefineIsoGammaQAHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doQA)
    return nullptr;

  const float M02Range[2] = {0, 5};
  const float M20Range[2] = {0, 5};
  const float pxyzRange[2] = {-50, 50};
  const float pTRange[2] = {0, 150};

  // QA plots for all IsoGammas
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  TH1F *hNIsoGamma = new TH1F("hNIsoGamma", "hNIsoGamma", 201, -0.5, 200.5);

  TH2F *hpTSpectraAfterSubsequentCuts = new TH2F("hpTSpectraAfterSubsequentCuts", "hpTSpectraAfterSubsequentCuts", 11, -0.5, 10.5, 300, pTRange[0], pTRange[1]);
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(1, "All");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(2, "Accepted");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(3, "E_{min}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(4, "Time");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(5, "N_{cells}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(6, "NLM");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(7, "Dist BC");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(8, "Track matching");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(9, "Exotic");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(10, "M_{02}");
  hpTSpectraAfterSubsequentCuts->GetXaxis()->SetBinLabel(11, "#it{p}_{T}^{iso}");

  TH2F *hpTSpectrumLossFromIndividualCuts = new TH2F("hpTSpectrumLossFromIndividualCuts", "hpTSpectrumLossFromIndividualCuts", 11, -0.5, 10.5, 300, pTRange[0], pTRange[1]);
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(1, "All");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(2, "Acceptance");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(3, "E_{min}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(4, "Time");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(5, "N_{cells}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(6, "NLM");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(7, "Dist BC");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(8, "Track matching");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(9, "Exotic");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(10, "M_{02}");
  hpTSpectrumLossFromIndividualCuts->GetXaxis()->SetBinLabel(11, "#it{p}_{T}^{iso}");

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
  TH2F *hIsoGammaRhoPt = new TH2F("hIsoGammaRhoPt", "hRhoPt", 600, 0, 300, 300, 0, 150);
  hIsoGammaRhoPt->GetXaxis()->SetTitle("#rho (GeV/#it{c})");
  hIsoGammaRhoPt->GetYaxis()->SetTitle("#it{p}_{T}^{trigger} (GeV/#it{c})");
  
  // replase EtaPhi later on with ThNsparseF. For now "double it and give it to the next person"
  THnSparseF *hIsoGammaEtaPhiPt = new THnSparseF("hIsoGammaEtaPhiPt", "hIsoGammaEtaPhiPt", 3, new int[3]{100, 100, 300}, new double[3]{-1, 0., 0}, new double[3]{1, 2 * TMath::Pi(), 150});
  hIsoGammaEtaPhiPt->Sumw2();
  hIsoGammaEtaPhiPt->GetAxis(0)->SetTitle("#eta");
  hIsoGammaEtaPhiPt->GetAxis(1)->SetTitle("#phi");
  hIsoGammaEtaPhiPt->GetAxis(2)->SetTitle("p_{T} (GeV/#it{c})");
  dir->Add(hIsoGammaEtaPhiPt);

  TH2F *hIsoGammaM02pT = new TH2F("hIsoGammaM02pT", "hIsoGammaM02pT", 100, M02Range[0], M02Range[1], 100, 0, 50);
  TH2F * hIsoGammaNCellsPt = new TH2F("hIsoGammaNCellsPt", "hIsoGammaNCellsPt", 100, 0, 100, 300, 0, 150);
  TH2F * hIsoGammaTimePt = new TH2F("hIsoGammaTimePt", "hIsoGammaTimePt", 300, -150, 150, 300, 0, 150);
  TH2F *hIsoGammaEBeforeAfterNL = new TH2F("hIsoGammaEBeforeAfterNL", "hIsoGammaEBeforeAfterNL;#bf{E_{cls}^{before} (GeV)};#bf{E_{cls}^{after}/E_{cls}^{before}}", 2000, 0, 200, 200, 0.8, 1.2);
  TH2F *hGGMinvDist = new TH2F("hGGMinvDist", "hGGMinvDist", 100, 0, 1, 500, 0, 50);
  TH2F *hIsoGammadEtadphi = new TH2F("hIsoGammadEtadphi", "hIsoGammadEtadphi", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2F *hIsoGammadEtadphiCut = new TH2F("hIsoGammadEtadphiCut", "hIsoGammadEtadphiCut", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2F *hDistanceToBadChannelvsPt = new TH2F("hDistanceToBadChannelvsPt", "hDistanceToBadChannelvsPt", 12, 0.5, 12.5, 100, 0, 100);
  THnSparseF * hIsoGammadEtadPhiPt = new THnSparseF("hIsoGammadEtadPhiPt", "hIsoGammadEtadPhiPt", 3, new int[3]{100, 100, 300}, new double[3]{-0.1, -0.1, 0}, new double[3]{0.1,0.1, 150});
  hIsoGammadEtadPhiPt->GetAxis(0)->SetTitle("d#eta");
  hIsoGammadEtadPhiPt->GetAxis(1)->SetTitle("d#phi");
  hIsoGammadEtadPhiPt->GetAxis(2)->SetTitle("p_{T} (GeV/#it{c})");
  hIsoGammadEtadPhiPt->Sumw2();
  dir->Add(hIsoGammadEtadPhiPt);

  if (optns.isMC)
  { // ------------------> MC IsoGammaQA histograms
    TH1F *hMinMassDiffToPi0Signal = new TH1F("hMinMassDiffToPi0Signal", "hMinMassDiffToPi0Signal", 50, 0., 0.1);
    TH1F *hMinMassDiffToPi0Background = new TH1F("hMinMassDiffToPi0Background", "hMinMassDiffToPi0Background", 50, 0., 0.1);

    TH2F *hNewVsOldSignalDef = new TH2F("hNewVsOldSignalDef", "hNewVsOldSignalDef", 2, -0.5, 1.5, 2, -0.5, 1.5);
    hNewVsOldSignalDef->GetXaxis()->SetBinLabel(2, "Signal Old");
    hNewVsOldSignalDef->GetXaxis()->SetBinLabel(1, "Background Old");
    hNewVsOldSignalDef->GetYaxis()->SetBinLabel(2, "Signal New");
    hNewVsOldSignalDef->GetYaxis()->SetBinLabel(1, "Background New");

    TH1F *hTrueIsoGammaMCTag = new TH1F("hTrueIsoGammaMCTag", "hMCTag", 20, 0.5, 20.5);
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(1, "Prompt #gamma");
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(2, "Frag #gamma");
    // Generated photons QA:
    if (dirname != "ClusterQA")
    {
      TH1F *hGammaGenIsoCharged = new TH1F("hGammaGenIsoCharged", "hGammaGenIsoCharged", 100, -10, 30);
      TH1F *hGammaGenBckPerb = new TH1F("hGammaGenBckPerb", "hGammaGenBckPerb", 100, 0, 200);
      TH2F *hGammaGenEtaPhi = new TH2F("hGammaGenEtaPhi", "hGammaGenEtaPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());
      TH2F *hGammaGenEtaPhiAcceptanceCut = new TH2F("hGammaGenEtaPhiAcceptanceCut", "hGammaGenEtaPhiAcceptanceCut", 100, -1., 1., 100, 0., 2 * TMath::Pi());
      //((TH2F *)dir->FindObject("hIsoGammaEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
    }

    // Signal IsoGammas
    TH1F *hNIsoGammaSignal = new TH1F("hNIsoGammaSignal", "hNIsoGammaSignal", 21, -0.5, 20.5);
    TH1F *hIsoGammaIsoChargedSignal = new TH1F("hIsoGammaIsoChargedSignal", "hIsoGammaIsoChargedSignal", 400, 0, 200);
    TH1F *hIsoGammaIsoChargedCorrectedSignal = new TH1F("hIsoGammaIsoChargedCorrectedSignal", "hIsoGammaIsoChargedCorrectedSignal", 400, -10, 30);
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
    TH2F *hDistanceToBadChannelvsPtSignal = new TH2F("hDistanceToBadChannelvsPtSignal", "hDistanceToBadChannelvsPtSignal", 12, 0.5, 12.5, 100, 0, 100);

    // NonSignal IsGammas
    TH1F *hNIsoGammaBackground = new TH1F("hNIsoGammaBackground", "hNIsoGammaBackground", 21, -0.5, 20.5);
    TH1F *hIsoGammaIsoChargedBackground = new TH1F("hIsoGammaIsoChargedBackground", "hIsoGammaIsoChargedBackground", 400, 0, 200);
    TH1F *hIsoGammaIsoChargedCorrectedBackground = new TH1F("hIsoGammaIsoChargedCorrectedBackground", "hIsoGammaIsoChargedCorrectedBackground", 400, -10, 30);
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
    TH2F *hDistanceToBadChannelvsPtBackground = new TH2F("hDistanceToBadChannelvsPtBackground", "hDistanceToBadChannelvsPtBackground", 12, 0.5, 12.5, 100, 0, 100);

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

TDirectory *DefineJetHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doJets)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
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

TDirectory *DefineJetQAHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doJets || !optns.doQA)
    return nullptr;
  const float pxyzRange[2] = {-50, 50};
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();
  TH1F *hNDLJets = new TH1F("hNDLJets", "hNDLJets", 21, -0.5, 20.5);
  TH1F *hJetPx = new TH1F("hJetPx", "hPx", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetPy = new TH1F("hJetPy", "hPy", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetPz = new TH1F("hJetPz", "hPz", 1000, pxyzRange[0], pxyzRange[1]);
  TH1F *hJetArea = new TH1F("hJetArea", "hArea", 600, 0., 3.);
  TH1F *hJetNch = new TH1F("hJetNch", "hNch", 51, -0.5, 50.5);
  TH1F *hJetNclus = new TH1F("hJetNclus", "hNclus", 51, -0.5, 50.5);
  TH1F *hJetNconstits = new TH1F("hJetNconstits", "hNconstits", 201, -0.5, 200.5);
  THnSparseF *hJetPtEtaPhi = new THnSparseF("hJetPtEtaPhi", "hJetPtEtaPhi", 3, new int[3]{350, 100, 100}, new double[3]{-50, -1, 0}, new double[3]{150, 1, 2 * TMath::Pi()});
  hJetPtEtaPhi->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  hJetPtEtaPhi->GetAxis(1)->SetTitle("#eta");
  hJetPtEtaPhi->GetAxis(2)->SetTitle("#phi");
  hJetPtEtaPhi->Sumw2();
  dir->Add(hJetPtEtaPhi);

  // define THNSparse that has three axes: jet pt, jet eta and z=leading hadron pt / jet pt
  THnSparseF *hJetPtEtaZ = new THnSparseF("hJetPtEtaZ", "hJetPtEtaZ", 3, new int[3]{300, 100, 100}, new double[3]{0, -1, 0}, new double[3]{150, 1, 2});
  hJetPtEtaZ->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  hJetPtEtaZ->GetAxis(1)->SetTitle("#eta");
  hJetPtEtaZ->GetAxis(2)->SetTitle("z = p_{T}^{leading hadron} / p_{T}^{jet}");
  hJetPtEtaZ->Sumw2();
  dir->Add(hJetPtEtaZ);

  THnSparseF *hJetPtEtaOccupancy = new THnSparseF("hJetPtEtaOccupancy", "hJetPtEtaOccupancy", 3, new int[3]{300, 100, 200}, new double[3]{0, -1, 0}, new double[3]{150, 1, 20000});
  hJetPtEtaOccupancy->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  hJetPtEtaOccupancy->GetAxis(1)->SetTitle("#eta");
  hJetPtEtaOccupancy->GetAxis(2)->SetTitle("Occupancy");
  hJetPtEtaOccupancy->Sumw2();
  dir->Add(hJetPtEtaOccupancy);

  // ------------------> MC Jet QA histograms
  if (optns.isMC)
  {
    TH1F *hNPLJets = new TH1F("hNPLJets", "hNPLJets", 21, -0.5, 20.5);
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

TDirectory *DefineGGPi0Histograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doGGPi0)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  TH2F *hGGMinvPt = new TH2F("hGGMinvPt", "hGGMinvPt", 400, 0., 0.4, 200, 0., 20.);

  return dir;
}

TDirectory *DefineGGPi0QAHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doGGPi0)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  TH2F *hGGMinvPt = new TH2F("hGGMinvPt", "hGGMinvPt", 400, 0., 0.4, 200, 0., 20.);

  return dir;
}

TDirectory *DefineGammaJetHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doJets)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  int nBinsPhotonPt = 300;
  int nBinsJetPt = 600;
  Double_t photonPtMinMax[2] = {0, 150};
  Double_t jetPtMinMax[2] = {-100, 200};

  TH2F *hpTImbalancevsDeltaPhi = new TH2F("hpTImbalancevsDeltaPhi", "hpTImbalancevsDeltaPhi;#bf{#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}};#bf{#Delta#phi = |#phi_{#gamma}-#phi_{jet}|}", 100, 0., 2., 314, 0, TMath::Pi());
  THnSparseF *hIsoGammaJetDeltaPhiJetPtGammaPt = new THnSparseF("hIsoGammaJetDeltaPhiJetPtGammaPt", "hIsoGammaJetDeltaPhiJetPtGammaPt", 3, new int[3]{100, nBinsJetPt, nBinsPhotonPt}, new double[3]{0, jetPtMinMax[0], photonPtMinMax[0]}, new double[3]{TMath::Pi(), jetPtMinMax[1], photonPtMinMax[1]});
  hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(0)->SetTitle("#Delta#phi = |#phi_{#gamma}-#phi_{jet}|");
  hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(1)->SetTitle("#it{p}_{T}^{jet} (GeV/c)");
  hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(2)->SetTitle("#it{p}_{T}^{#gamma} (GeV/c)");
  hIsoGammaJetDeltaPhiJetPtGammaPt->Sumw2();
  dir->Add(hIsoGammaJetDeltaPhiJetPtGammaPt);
  int nBinsArray[3] = {100, nBinsJetPt, nBinsPhotonPt};
  double minArray[3] = {0, jetPtMinMax[0], photonPtMinMax[0]};
  double maxArray[3] = {2*TMath::Pi(), jetPtMinMax[1], photonPtMinMax[1]};
  
  THnSparseF *hIsoGammaJetDeltaPhi2piJetPtGammaPt = new THnSparseF("hIsoGammaJetDeltaPhi2piJetPtGammaPt", 
                                                                    "hIsoGammaJetDeltaPhi2piJetPtGammaPt", 
                                                                    3, 
                                                                    nBinsArray,
                                                                    minArray,
                                                                    maxArray);
  hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(0)->SetTitle("#Delta#phi = |#phi_{#gamma}-#phi_{jet}|");
  hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(1)->SetTitle("#it{p}_{T}^{jet} (GeV/c)");
  hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(2)->SetTitle("#it{p}_{T}^{#gamma} (GeV/c)");
  hIsoGammaJetDeltaPhi2piJetPtGammaPt->Sumw2();
  dir->Add(hIsoGammaJetDeltaPhi2piJetPtGammaPt);

  int nBinsDeltaDeltaPhi[3] = {100, 100, nBinsJetPt};
  double minDeltaDeltaPhi[3] = {0, 0, jetPtMinMax[0]};
  double maxDeltaDeltaPhi[3] = {2*TMath::Pi(), 2,jetPtMinMax[1]};

  THnSparseF *hIsoGammaJetDeltaPhi2piDeltaEtaJetPt = new THnSparseF("hIsoGammaJetDeltaPhi2piDeltaEtaJetPt", "hIsoGammaJetDeltaPhi2piDeltaEtaJetPt", 3, nBinsDeltaDeltaPhi, minDeltaDeltaPhi, maxDeltaDeltaPhi);
  hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(0)->SetTitle("#Delta#phi = |#phi_{#gamma}-#phi_{jet}|");
  hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(1)->SetTitle("#Delta#eta = |#eta_{#gamma}-#eta_{jet}|");
  hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(2)->SetTitle("#it{p}_{T}^{jet} (GeV/c)");
  hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->Sumw2();
  dir->Add(hIsoGammaJetDeltaPhi2piDeltaEtaJetPt);
  return dir;
}

TDirectory *DefineGammaJetQAHistograms(TFile *f, std::string dirname, GlobalOptions optns)
{
  if (!optns.doIsoGamma || !optns.doJets || !optns.doQA)
    return nullptr;
  TDirectory *dir = f->mkdir(dirname.c_str());
  dir->cd();

  int nBinsJetPt = 600;
  Double_t jetPtMinMax[2] = {-100, 200};

  TH2F *hNIsoGammaJets = new TH2F("hNIsoGammaJets", "hNIsoGammaJets;#bf{N_{#gamma}};#bf{N_{jet}}", 11, -0.5, 10.5, 11, -0.5, 10.5);
  TH2F *hDeltaEtavsDeltaPhi = new TH2F("hDeltaEtavsDeltaPhi", "hDeltaEtavsDeltaPhi;#bf{#Delta#eta = |#eta_{#gamma}-#eta_{jet}|};#bf{#Delta#phi = |#phi_{#gamma}-#phi_{jet}|}", 100, 0., 2., 314, 0, TMath::Pi());
  TH2F *hIsoGammaFromCorrelationEtaPhiMap = new TH2F("hIsoGammaFromCorrelationEtaPhiMap", "hIsoGammaFromCorrelationEtaPhiMap;#bf{#eta_{iso #gamma}};#bf{#phi_{iso #gamma}}", 100, -1, 1., 628, 0, 2 * TMath::Pi());
  TH2F *hJetFromCorrelationEtaPhiMap = new TH2F("hJetFromCorrelationEtaPhiMap", "hJetFromCorrelationEtaPhiMap;#bf{#eta_{Jet}};#bf{#phi_{Jet}}", 100, -1, 1., 628, 0, 2 * TMath::Pi());

  TH1F *hpTImbalance = new TH1F("hpTImbalance", "hpTImbalance;#bf{#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}}", 100, 0., 2.);
  TH1F *hDeltaPhi = new TH1F("hDeltaPhi", "hDeltaPhi;#bf{#Delta#phi = |#phi_{#gamma}-#phi_{jet}|}", 314, 0, TMath::Pi());

  THnSparseF *hJetEtaPtDeltaPhi = new THnSparseF("hJetEtaPtDeltaPhi", "hJetEtaPtDeltaPhi", 3, new int[3]{100, nBinsJetPt, 100}, new double[3]{-1, jetPtMinMax[0], 0}, new double[3]{1, jetPtMinMax[1], 2*TMath::Pi()});
  hJetEtaPtDeltaPhi->GetAxis(0)->SetTitle("#eta");
  hJetEtaPtDeltaPhi->GetAxis(1)->SetTitle("#it{p}_{T}^{jet} (GeV/c)");
  hJetEtaPtDeltaPhi->GetAxis(2)->SetTitle("#Delta#phi = |#phi_{#gamma}-#phi_{jet}|");
  hJetEtaPtDeltaPhi->Sumw2();
  dir->Add(hJetEtaPtDeltaPhi);

  // jet pt vs eta vs jet constituents, only for back-to-back jets (delta 0.5pi < dphi < 1.5pi)
  THnSparseF *hJetEtaPtConstBackToBack = new THnSparseF("hJetEtaPtConstBackToBack", "hJetEtaPtConstBackToBack", 3, new int[3]{100, nBinsJetPt, 150}, new double[3]{-1, jetPtMinMax[0], 0}, new double[3]{1, jetPtMinMax[1], 150});
  hJetEtaPtConstBackToBack->GetAxis(0)->SetTitle("#eta");
  hJetEtaPtConstBackToBack->GetAxis(1)->SetTitle("#it{p}_{T}^{jet} (GeV/c)");
  hJetEtaPtConstBackToBack->GetAxis(2)->SetTitle("N jet constituens");
  hJetEtaPtConstBackToBack->Sumw2();
  dir->Add(hJetEtaPtConstBackToBack);

  return dir;
}

void fillHistograms(std::vector<Cluster> obj, Event evt,TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns, IsoGammaCuts IsoGammaCuts)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH1F *)dir->FindObject("hIsoGammaPt"))->Fill(obj.at(i).Pt(), eventWeight);
    ((TH1F *)dir->FindObject("hIsoGammaE"))->Fill(obj.at(i).E, eventWeight);
    ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPt"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
    if (optns.isMC)
    {
      if (IsoGammaCuts.isSignal(obj.at(i)))
      {
        ((TH1F *)dir->FindObject("hIsoGammaPtSignal"))->Fill(obj.at(i).Pt(), eventWeight);
        ((TH1F *)dir->FindObject("hIsoGammaESignal"))->Fill(obj.at(i).E, eventWeight);
        ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPtSignal"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
      }
      else
      {
        ((TH1F *)dir->FindObject("hIsoGammaPtBackground"))->Fill(obj.at(i).Pt(), eventWeight);
        ((TH1F *)dir->FindObject("hIsoGammaEBackground"))->Fill(obj.at(i).E, eventWeight);
        ((THnSparseF *)dir->FindObject("hIsoGammaIsovsM02vsPtBackground"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
      }
    }
  }
  if (optns.doQA)
  {
    ((TH1F *)QAdir->FindObject("hNIsoGamma"))->Fill((int)obj.size(), eventWeight);
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)QAdir->FindObject("hIsoGammaIsoCharged"))->Fill(obj.at(i).IsoCharged, eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaIsoChargedCorrected"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaPx"))->Fill(obj.at(i).Px(), eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaPy"))->Fill(obj.at(i).Py(), eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaPz"))->Fill(obj.at(i).Pz(), eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaE"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaM02"))->Fill(obj.at(i).M02, eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaM20"))->Fill(obj.at(i).M20, eventWeight);
      ((TH1F *)QAdir->FindObject("hIsoGammaNLM"))->Fill(obj.at(i).NLM, eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaNCellsPt"))->Fill(obj.at(i).NCells, obj.at(i).Pt(), eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaTimePt"))->Fill(obj.at(i).Time, obj.at(i).Pt(), eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      // hIsoGammaEtaPhi
      double coords[3] = {obj.at(i).Eta(), obj.at(i).Phi(), obj.at(i).Pt()};
      ((THnSparseF*) QAdir->FindObject("hIsoGammaEtaPhiPt"))->Fill(coords, eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaM02pT"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaEBeforeAfterNL"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
      ((TH2F *)QAdir->FindObject("hDistanceToBadChannelvsPt"))->Fill(obj.at(i).DistanceToBadChannel, obj.at(i).Pt(), eventWeight);
      if (obj.at(i).MatchedTrack.P > 0)
      {
        ((TH2F *)QAdir->FindObject("hIsoGammadEtadphi"))->Fill(obj.at(i).MatchedTrack.dPhi, obj.at(i).MatchedTrack.dEta, eventWeight);
        double coords[3] = {obj.at(i).MatchedTrack.dPhi, obj.at(i).MatchedTrack.dEta, obj.at(i).Pt()};
        ((THnSparseF *)QAdir->FindObject("hIsoGammadEtadPhiPt"))->Fill(coords, eventWeight);
      }
      ((TH2F *)QAdir->FindObject("hIsoGammaRhoPt"))->Fill(evt.Rho, obj.at(i).Pt(), eventWeight);

      if (optns.isMC)
      {
        double isSigOld = 0;
        double isSigNew = 0;
        if (IsoGammaCuts.isSignalClusterLevelIso(obj.at(i)))
        {
          isSigOld = 1;
        }
        if (IsoGammaCuts.isSignal(obj.at(i)))
        {
          isSigNew = 1;
        }
        ((TH2F *)QAdir->FindObject("hNewVsOldSignalDef"))->Fill(isSigOld, isSigNew);
        ((TH2F *)QAdir->FindObject("hTrueIsoGammaMCTag"))->Fill(obj.at(i).MCTag, eventWeight);
        if (IsoGammaCuts.isSignal(obj.at(i)))
        {
          ((TH1F *)QAdir->FindObject("hIsoGammaIsoChargedSignal"))->Fill(obj.at(i).IsoCharged, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaIsoChargedCorrectedSignal"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
          ((TH1F *)QAdir->FindObject("hNIsoGammaSignal"))->Fill((int)obj.size(), eventWeight);
          ((TH1F *)QAdir->FindObject("hMinMassDiffToPi0Signal"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPxSignal"))->Fill(obj.at(i).Px(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPySignal"))->Fill(obj.at(i).Py(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPzSignal"))->Fill(obj.at(i).Pz(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaESignal"))->Fill(obj.at(i).E, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaM02Signal"))->Fill(obj.at(i).M02, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaM20Signal"))->Fill(obj.at(i).M20, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaNLMSignal"))->Fill(obj.at(i).NLM, eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaEtaPhiSignal"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTSignal"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaEBeforeAfterNLSignal"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
          ((TH2F *)QAdir->FindObject("hDistanceToBadChannelvsPtSignal"))->Fill(obj.at(i).DistanceToBadChannel, obj.at(i).Pt(), eventWeight);
        }
        else
        {
          ((TH1F *)QAdir->FindObject("hIsoGammaIsoChargedBackground"))->Fill(obj.at(i).IsoCharged, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaIsoChargedCorrectedBackground"))->Fill(obj.at(i).IsoChargedCorrected, eventWeight);
          ((TH1F *)QAdir->FindObject("hNIsoGammaBackground"))->Fill((int)obj.size(), eventWeight);
          ((TH1F *)QAdir->FindObject("hMinMassDiffToPi0Background"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPxBackground"))->Fill(obj.at(i).Px(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPyBackground"))->Fill(obj.at(i).Py(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaPzBackground"))->Fill(obj.at(i).Pz(), eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaEBackground"))->Fill(obj.at(i).E, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaM02Background"))->Fill(obj.at(i).M02, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaM20Background"))->Fill(obj.at(i).M20, eventWeight);
          ((TH1F *)QAdir->FindObject("hIsoGammaNLMBackground"))->Fill(obj.at(i).NLM, eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaEtaPhiBackground"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTBackground"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
          ((TH2F *)QAdir->FindObject("hIsoGammaEBeforeAfterNLBackground"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
          ((TH2F *)QAdir->FindObject("hDistanceToBadChannelvsPtBackground"))->Fill(obj.at(i).DistanceToBadChannel, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isPhoton())
        {
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTAllGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isGammaFromPion())
        {
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTPionDecayGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isGammaFromEta())
        {
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTEtaDecayGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
        if (obj.at(i).isPion())
        {
          ((TH2F *)QAdir->FindObject("hIsoGammaM02pTMergedPionGamma"))->Fill(obj.at(i).M02, obj.at(i).Pt(), eventWeight);
        }
      }
      else
      {
        ((TH1F *)QAdir->FindObject("hMinMassDiffToPi0"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
      }
    }
  }
}

void fillHistograms(std::vector<DLJet> obj, Event evt, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH1F *)dir->FindObject("hJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
  }
  if (optns.doQA)
  {
    ((TH1F *)QAdir->FindObject("hNDLJets"))->Fill(obj.size(), eventWeight);
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)QAdir->FindObject("hJetPx"))->Fill(obj.at(i).Px(), eventWeight);
      ((TH1F *)QAdir->FindObject("hJetPy"))->Fill(obj.at(i).Py(), eventWeight);
      ((TH1F *)QAdir->FindObject("hJetPz"))->Fill(obj.at(i).Pz(), eventWeight);
      ((TH1F *)QAdir->FindObject("hJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F *)QAdir->FindObject("hJetNch"))->Fill(obj.at(i).Nch, eventWeight);
      ((TH1F *)QAdir->FindObject("hJetNclus"))->Fill(obj.at(i).Nclus, eventWeight);
      ((TH1F *)QAdir->FindObject("hJetNconstits"))->Fill(obj.at(i).Nconstits, eventWeight);
      ((THnSparseF *)QAdir->FindObject("hJetPtEtaPhi"))->Fill(obj.at(i).Pt(),obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);

      // define THNSparse that has three axes: jet pt, jet eta and z=leading hadron pt / jet pt
      if (obj.at(i).Pt() > 0) {
        double coords[3] = {obj.at(i).Pt(), obj.at(i).Eta(), obj.at(i).LeadingTrackPt / obj.at(i).Pt()};
        ((THnSparseF *)QAdir->FindObject("hJetPtEtaZ"))->Fill(coords, eventWeight);
      }
      // hJetPtEtaOccupancy
      double coords[3] = {obj.at(i).Pt(), obj.at(i).Eta(), static_cast<double>(evt.Occupancy)};
      ((THnSparseF *)QAdir->FindObject("hJetPtEtaOccupancy"))->Fill(coords, eventWeight);
    }
  }
}

void fillHistograms(std::vector<PLJet> obj, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH1F *)dir->FindObject("hPLJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
    if (obj.at(i).ClosestDLJet != nullptr)
    {
      ((TH2F *)dir->FindObject("hDLJetPtVsPLJetPt"))->Fill(obj.at(i).ClosestDLJet->Pt(), obj.at(i).Pt(), eventWeight);
      ((TH2F *)dir->FindObject("hDLSubPLDivPLJetPtVsPLJetPt"))->Fill((obj.at(i).ClosestDLJet->Pt() - obj.at(i).Pt()) / obj.at(i).Pt(), obj.at(i).Pt(), eventWeight);
    }
  }
  if (optns.doQA)
  {
    ((TH1F *)QAdir->FindObject("hNPLJets"))->Fill(obj.size(), eventWeight);
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)QAdir->FindObject("hPLJetPx"))->Fill(obj.at(i).Px(), eventWeight);
      ((TH1F *)QAdir->FindObject("hPLJetPy"))->Fill(obj.at(i).Py(), eventWeight);
      ((TH1F *)QAdir->FindObject("hPLJetPz"))->Fill(obj.at(i).Pz(), eventWeight);
      ((TH1F *)QAdir->FindObject("hPLJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F *)QAdir->FindObject("hPLJetNPart"))->Fill(obj.at(i).NPart, eventWeight);
      ((TH2F *)QAdir->FindObject("hPLJetEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      if (obj.at(i).ClosestDLJet != nullptr)
      {
        ((TH2F *)QAdir->FindObject("hDLJetPVsPLJetP"))->Fill(obj.at(i).ClosestDLJet->P(), obj.at(i).P(), eventWeight);
        ((TH2F *)QAdir->FindObject("hDLJetPxVsPLJetPx"))->Fill(fabs(obj.at(i).ClosestDLJet->Px()), fabs(obj.at(i).Px()), eventWeight);
        ((TH2F *)QAdir->FindObject("hDLJetPyVsPLJetPy"))->Fill(fabs(obj.at(i).ClosestDLJet->Py()), fabs(obj.at(i).Py()), eventWeight);
      }
    }
  }
}

void fillHistograms(std::vector<Pi0> obj, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH2F *)dir->FindObject("hGGMinvPt"))->Fill(obj.at(i).Mass, obj.at(i).Pt(), eventWeight);
  }
  if (optns.doQA)
  {
    if (((TString)dir->GetName()) == "IsoGammaQA")
    {
      for (unsigned long i = 0; i < obj.size(); i++)
      {
        ((TH2F *)dir->FindObject("hGGMinvDist"))->Fill(obj.at(i).Mass, obj.at(i).Pt(), eventWeight);
      }
    }
  }
}

void fillHistograms(std::vector<GammaJetPair> obj, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH2F *)dir->FindObject("hpTImbalancevsDeltaPhi"))->Fill(obj.at(i).pTImbalance, obj.at(i).DPhi, eventWeight);
    double coords[3] = {obj.at(i).DPhi, obj.at(i).jet->Pt(), obj.at(i).isoGamma->Pt()};
    ((THnSparseF *)dir->FindObject("hIsoGammaJetDeltaPhiJetPtGammaPt"))->Fill(coords, eventWeight);
    coords[0] = obj.at(i).DPhi2pi;
    ((THnSparseF *)dir->FindObject("hIsoGammaJetDeltaPhi2piJetPtGammaPt"))->Fill(coords, eventWeight);
    coords[0] = obj.at(i).DPhi2pi;
    coords[1] = TMath::Abs(obj.at(i).DEta);
    coords[2] = obj.at(i).jet->Pt();
    ((THnSparseF *)dir->FindObject("hIsoGammaJetDeltaPhi2piDeltaEtaJetPt"))->Fill(coords, eventWeight);
  }
  if (optns.doQA)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH2F *)QAdir->FindObject("hDeltaEtavsDeltaPhi"))->Fill(obj.at(i).DEta, obj.at(i).DPhi, eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaFromCorrelationEtaPhiMap"))->Fill(obj.at(i).isoGamma->Eta(), obj.at(i).isoGamma->Phi(), eventWeight);
      ((TH2F *)QAdir->FindObject("hJetFromCorrelationEtaPhiMap"))->Fill(obj.at(i).jet->Eta(), obj.at(i).jet->Phi(), eventWeight);
      ((TH1F *)QAdir->FindObject("hpTImbalance"))->Fill(obj.at(i).pTImbalance, eventWeight);
      ((TH1F *)QAdir->FindObject("hDeltaPhi"))->Fill(obj.at(i).DPhi, eventWeight);
      double coords[3] = {obj.at(i).jet->Eta(), obj.at(i).jet->Pt(), obj.at(i).DPhi2pi};
      ((THnSparseF *)QAdir->FindObject("hJetEtaPtDeltaPhi"))->Fill(coords, eventWeight);
      if(obj.at(i).DPhi2pi > 0.5 * TMath::Pi() && obj.at(i).DPhi2pi < 1.5 * TMath::Pi()){
        double coords[3] = {obj.at(i).jet->Eta(), obj.at(i).jet->Pt(), static_cast<double>(obj.at(i).jet->Nconstits)};
        ((THnSparseF *)QAdir->FindObject("hJetEtaPtConstBackToBack"))->Fill(coords, eventWeight);
      }
    }
  }
}

void fillHistograms(std::vector<GGPi0JetPair> obj, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH2F *)dir->FindObject("hpTImbalancevsDeltaPhi"))->Fill(obj.at(i).pTImbalance, obj.at(i).DPhi, eventWeight);
  }
  if (optns.doQA)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH2F *)QAdir->FindObject("hDeltaEtavsDeltaPhi"))->Fill(obj.at(i).DEta, obj.at(i).DPhi, eventWeight);
      ((TH2F *)QAdir->FindObject("hIsoGammaFromCorrelationEtaPhiMap"))->Fill(obj.at(i).ggPi0->Eta(), obj.at(i).ggPi0->Phi(), eventWeight);
      ((TH2F *)QAdir->FindObject("hJetFromCorrelationEtaPhiMap"))->Fill(obj.at(i).jet->Eta(), obj.at(i).jet->Phi(), eventWeight);
      ((TH1F *)QAdir->FindObject("hpTImbalance"))->Fill(obj.at(i).pTImbalance, eventWeight);
      ((TH1F *)QAdir->FindObject("hDeltaPhi"))->Fill(obj.at(i).DPhi, eventWeight);
    }
  }
}

void fillHistograms(Event obj, TDirectory *dir, TDirectory *QAdir, float eventWeight,EventCuts eventCuts, GlobalOptions optns)
{
  ((TH1F *)dir->FindObject("hNEvents"))->Fill(0., (Double_t)eventWeight);
  if(obj.IsTriggerAlias(kTVXinEMC))((TH1F *)dir->FindObject("hNEvents"))->Fill(1., eventWeight);
  if(obj.IsTriggerAlias(kEMC7))((TH1F *)dir->FindObject("hNEvents"))->Fill(2., eventWeight);
  if(obj.IsTriggerAlias(kDMC7))((TH1F *)dir->FindObject("hNEvents"))->Fill(3., eventWeight);
  if(obj.Occupancy >= eventCuts.GetOccupancyMin() && obj.Occupancy <= eventCuts.GetOccupancyMax()){
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(4., eventWeight);
  }
  if(eventCuts.PassedCuts(obj)){
    ((TH1F *)dir->FindObject("hNEvents"))->Fill(5., eventWeight);
  }
  if (optns.doQA)
  {
    ((TH1F *)QAdir->FindObject("hZVtx"))->Fill(obj.ZVtx, eventWeight);
    ((TH1F *)QAdir->FindObject("hWeights"))->Fill(obj.weight);
    ((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(0., 1);
    if(obj.IsTriggerAlias(kTVXinEMC))((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(1., 1);
    if(obj.IsTriggerAlias(kEMC7))((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(2., 1);
    if(obj.IsTriggerAlias(kDMC7))((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(3., 1);
    if(obj.Occupancy > eventCuts.GetOccupancyMin() && obj.Occupancy < eventCuts.GetOccupancyMax()){
      ((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(4., 1);
    }
    if(eventCuts.PassedCuts(obj)){
      ((TH1F *)QAdir->FindObject("hNEventswoWeights"))->Fill(5., 1);
    }
    ((THnSparseF *)QAdir->FindObject("hCentralityRhoMultiplicity"))->Fill(obj.Centrality, obj.Rho, obj.Multiplicity, eventWeight);
    ((TH1F *)QAdir->FindObject("hOccupancy"))->Fill(obj.Occupancy, eventWeight);
    ((TH2F *)QAdir->FindObject("hOccupancyCentrality"))->Fill(obj.Centrality, obj.Occupancy, eventWeight);
  }

  // loop over all triggerAliases and fill the histogram
  // for (int i = 0; i < triggerAliases::kNaliases; i++) {
  //   if (obj.IsTriggerAlias(static_cast<triggerAliases>(i))) {
  //     ((TH1F *)QAdir->FindObject("hNEventsTriggered"))->Fill(i);
  //   }
  // }
}

template <typename T>
void fillHistograms(T obj, TDirectory *dir, TDirectory *QAdir, float eventWeight, GlobalOptions optns, GammaGenCuts GammaGenCuts)
{
  for (unsigned long i = 0; i < obj.size(); i++)
  {
    ((TH1F *)dir->FindObject("hGammaGenE"))->Fill(obj.at(i).E, eventWeight);
    ((TH1F *)dir->FindObject("hGammaGenPt"))->Fill(obj.at(i).Pt(), eventWeight);
    if (GammaGenCuts.isSignal(obj.at(i)))
    {
      ((TH1F *)dir->FindObject("hGammaGenESignal"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenPtSignal"))->Fill(obj.at(i).Pt(), eventWeight);
    }
    else
    {
      ((TH1F *)dir->FindObject("hGammaGenEBackground"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenPtBackground"))->Fill(obj.at(i).Pt(), eventWeight);
    }
    if (GammaGenCuts.PassedGammaGenCuts(obj.at(i)))
    {
      ((TH1F *)dir->FindObject("hGammaGenEAcceptanceCut"))->Fill(obj.at(i).E, eventWeight);
      ((TH1F *)dir->FindObject("hGammaGenPtAcceptanceCut"))->Fill(obj.at(i).Pt(), eventWeight);
    }
  }
  if (optns.doQA)
  {
    for (unsigned long i = 0; i < obj.size(); i++)
    {
      ((TH1F *)QAdir->FindObject("hGammaGenIsoCharged"))->Fill(obj.at(i).IsoCharged, eventWeight);
      ((TH1F *)QAdir->FindObject("hGammaGenBckPerb"))->Fill(obj.at(i).IsoBckPerp, eventWeight);
      ((TH2F *)QAdir->FindObject("hGammaGenEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      if (GammaGenCuts.PassedGammaGenCuts(obj.at(i)))
      {
        ((TH2F *)QAdir->FindObject("hGammaGenEtaPhiAcceptanceCut"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      }
    }
  }
}
