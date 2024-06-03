#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TKey.h"
#include "Logging.h"
#include "TFile.h"
#include "Utilities.h"

TDirectory* DefineEventHistograms(TFile* f, GlobalOptions optns)
{
  TDirectory* dir = f->mkdir("Events");
  dir->cd();
  TH1F* hNEvents = new TH1F("hNEvents", "hNEvents", 100, -20, 20);
  return dir;
}

TDirectory* DefineEventQAHistograms(TFile* f, GlobalOptions optns)
{
  TDirectory* dir = f->mkdir("EventQA");
  dir->cd();
  TH1F* hZVtx = new TH1F("hZVtx", "hZVtx; Vtx_{z} (cm)", 100, -20, 20);

  const int nWeightBins = 100;
  Double_t weightBinning[nWeightBins + 1];
  setExpArray(weightBinning, nWeightBins + 1, 1E-9, 1E-2);
  TH1F* hWeights = new TH1F("hWeights", "hWeights;Weight", nWeightBins, weightBinning);

  return dir;
}

TDirectory* DefineIsoGammaHistograms(TFile* f, GlobalOptions optns)
{

  TDirectory* dir = f->mkdir("IsoGammas");
  dir->cd();

  TH1F* hIsoGammaPt = new TH1F("hIsoGammaPt", "hPt", 1000, 0., 100.);
  TH1F* hIsoGammaE = new TH1F("hIsoGammaE", "hE", 1000, 0., 100.);

  int nBins[3] = {1000, 200, 2000};
  double xmin[3] = {-20, 0, 0};
  double xmax[3] = {80, 2, 200};
  THnSparseF* hIsoGammaIsovsM02vsPt = new THnSparseF("hIsoGammaIsovsM02vsPt", "hIsoGammaIsovsM02vsPt", 3, nBins, xmin, xmax);
  hIsoGammaIsovsM02vsPt->Sumw2();
  dir->Add(hIsoGammaIsovsM02vsPt);

  return dir;
}

TDirectory* DefineIsoGammaQAHistograms(TFile* f, GlobalOptions optns)
{

  TDirectory* dir = f->mkdir("IsoGammaQA");
  dir->cd();

  TH1F* hNIsoGamma = new TH1F("hNIsoGamma", "hNIsoGamma", 21, -0.5, 20.5);

  TH1F* hIsoGammaPx = new TH1F("hIsoGammaPx", "hPx", 1000, -50., 50.);
  TH2F* hIsoGammaEtaPhi = new TH2F("hIsoGammaEtaPhi", "hEtaPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());

  TH2F* hIsoGammaEBeforeAfterNL = new TH2F("hIsoGammaEBeforeAfterNL", "hIsoGammaEBeforeAfterNL;#bf{E_{cls}^{before} (GeV)};#bf{E_{cls}^{after}/E_{cls}^{before}}", 2000, 0, 200, 200, 0.8, 1.2);

  TH2F* hGGMinvDist = new TH2F("hGGMinvDist", "hGGMinvDist", 100, 0, 1, 500, 0, 50);

  if (optns.isMC) { // ------------------> MC IsoGammaQA histograms
    TH1F* hMinMassDiffToPi0Signal = new TH1F("hMinMassDiffToPi0Signal", "hMinMassDiffToPi0Signal", 50, 0., 0.1);
    TH1F* hMinMassDiffToPi0Background = new TH1F("hMinMassDiffToPi0Background", "hMinMassDiffToPi0Background", 50, 0., 0.1);

    TH1F* hTrueIsoGammaMCTag = new TH1F("hTrueIsoGammaMCTag", "hMCTag", 20, 0.5, 20.5);
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(1, "Prompt #gamma");
    hTrueIsoGammaMCTag->GetXaxis()->SetBinLabel(2, "Frag #gamma");
  } else { // ------------------> Only in data Jet histograms
    TH1F* hMinMassDiffToPi0 = new TH1F("hMinMassDiffToPi0", "hMinMassDiffToPi0", 50, 0., 0.1);
  }

  return dir;
}

TDirectory* DefineJetHistograms(TFile* f, GlobalOptions optns)
{
  TDirectory* dir = f->mkdir("Jets");
  dir->cd();

  TH1F* hJetPt = new TH1F("hJetPt", "hPt", 1000, 0., 100.);

  // ------------------> MC Jet histograms
  if (optns.isMC) {
    TH1F* hPLJetPt = new TH1F("hPLJetPt", "hPt", 1000, 0., 100.);

    TH2F* hDLJetPtVsPLJetPt = new TH2F("hDLJetPtVsPLJetPt", "hPt;#it{p}_{T}^{DL Jet} (GeV/#it{c});#it{p}_{T}^{PL Jet} (GeV/#it{c})", 400, 0., 200., 400, 0., 200.);
  }

  return dir;
}

TDirectory* DefineJetQAHistograms(TFile* f, GlobalOptions optns)
{
  TDirectory* dir = f->mkdir("JetQA");
  dir->cd();
  TH1F* hJetPx = new TH1F("hJetPx", "hPx", 1000, -50., 50.);
  TH1F* hJetPy = new TH1F("hJetPy", "hPy", 1000, -50., 50.);
  TH1F* hJetPz = new TH1F("hJetPz", "hPz", 1000, -50., 50.);
  TH1F* hJetArea = new TH1F("hJetArea", "hArea", 100, 0., 1.);
  TH1F* hJetNch = new TH1F("hJetNch", "hNch", 51, -0.5, 50.5);
  TH1F* hJetNclus = new TH1F("hJetNclus", "hNclus", 51, -0.5, 50.5);
  TH2F* hJetEtaPhi = new TH2F("hJetEtaPhi", "hJetPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());

  // ------------------> MC Jet QA histograms
  if (optns.isMC) {
    TH1F* hPLJetPx = new TH1F("hPLJetPx", "hPx", 1000, -50., 50.);
    TH1F* hPLJetPy = new TH1F("hPLJetPy", "hPy", 1000, -50., 50.);
    TH1F* hPLJetPz = new TH1F("hPLJetPz", "hPz", 1000, -50., 50.);
    TH1F* hPLJetArea = new TH1F("hPLJetArea", "hArea", 100, 0., 1.);
    TH1F* hPLJetNPart = new TH1F("hPLJetNPart", "hNPart", 51, -0.5, 50.5);
    TH2F* hPLJetEtaPhi = new TH2F("hPLJetEtaPhi", "hPLJetPhi", 100, -1., 1., 100, 0., 2 * TMath::Pi());

    TH2F* hDLJetPVsPLJetP = new TH2F("hDLJetPVsPLJetP", "hPt;#it{p}_{T}^{DL Jet} (GeV/#it{c});#it{p}_{T}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
    TH2F* hDLJetPxVsPLJetPx = new TH2F("hDLJetPxVsPLJetPx", "hPx;#it{p}_{x}^{DL Jet} (GeV/#it{c});#it{p}_{x}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
    TH2F* hDLJetPyVsPLJetPy = new TH2F("hDLJetPyVsPLJetPy", "hPx;#it{p}_{y}^{DL Jet} (GeV/#it{c});#it{p}_{y}^{PL Jet} (GeV/#it{c})", 1000, 0., 100., 1000, 0., 100.);
  }

  return dir;
}

template <typename T>
void fillHistograms(T obj, TDirectory* dir, float eventWeight)
{
  if constexpr (std::is_same<T, std::vector<IsoGamma>>::value) {
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hIsoGammaPt"))->Fill(obj.at(i).Pt(), eventWeight);
      ((TH1F*)dir->FindObject("hIsoGammaE"))->Fill(obj.at(i).E, eventWeight);
      ((THnSparseF*)dir->FindObject("hIsoGammaIsovsM02vsPt"))->Fill(obj.at(i).IsoChargedCorrected, obj.at(i).M02, obj.at(i).Pt(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, std::vector<Jet>>::value) {
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, std::vector<PLJet>>::value) {
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hPLJetPt"))->Fill(obj.at(i).Pt(), eventWeight);
      if (obj.at(i).ClosestDLJet != nullptr)
        ((TH2F*)dir->FindObject("hDLJetPtVsPLJetPt"))->Fill(obj.at(i).ClosestDLJet->Pt(), obj.at(i).Pt(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, Event>::value) {
    ((TH1F*)dir->FindObject("hNEvents"))->Fill(1, eventWeight);
  }
}

template <typename T>
void fillQAHistograms(T obj, TDirectory* dir, float eventWeight, GlobalOptions optns)
{
  if constexpr (std::is_same<T, std::vector<IsoGamma>>::value) {
    ((TH1F*)dir->FindObject("hNIsoGamma"))->Fill((int)obj.size(), eventWeight);
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hIsoGammaPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH2F*)dir->FindObject("hIsoGammaEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      ((TH2F*)dir->FindObject("hIsoGammaEBeforeAfterNL"))->Fill(obj.at(i).EBeforeNL, obj.at(i).E / obj.at(i).EBeforeNL, eventWeight);
      if (optns.isMC) {
        ((TH2F*)dir->FindObject("hTrueIsoGammaMCTag"))->Fill(obj.at(i).MCTag, eventWeight);
        if (obj.at(i).isSignal())
          ((TH1F*)dir->FindObject("hMinMassDiffToPi0Signal"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
        else
          ((TH1F*)dir->FindObject("hMinMassDiffToPi0Background"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
      } else {
        ((TH1F*)dir->FindObject("hMinMassDiffToPi0"))->Fill(obj.at(i).MinMassDiffToPi0, eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<Pi0>>::value) {
    if (((TString)dir->GetName()) == "IsoGammaQA") {
      for (unsigned long i = 0; i < obj.size(); i++) {
        ((TH2F*)dir->FindObject("hGGMinvDist"))->Fill(obj.at(i).Mass, obj.at(i).Pt(), eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, std::vector<Jet>>::value) {
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hJetPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH1F*)dir->FindObject("hJetPy"))->Fill(obj.at(i).Py, eventWeight);
      ((TH1F*)dir->FindObject("hJetPz"))->Fill(obj.at(i).Pz, eventWeight);
      ((TH1F*)dir->FindObject("hJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F*)dir->FindObject("hJetNch"))->Fill(obj.at(i).Nch, eventWeight);
      ((TH1F*)dir->FindObject("hJetNclus"))->Fill(obj.at(i).Nclus, eventWeight);
      ((TH2F*)dir->FindObject("hJetEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
    }
  }
  if constexpr (std::is_same<T, std::vector<PLJet>>::value) {
    for (unsigned long i = 0; i < obj.size(); i++) {
      ((TH1F*)dir->FindObject("hPLJetPx"))->Fill(obj.at(i).Px, eventWeight);
      ((TH1F*)dir->FindObject("hPLJetPy"))->Fill(obj.at(i).Py, eventWeight);
      ((TH1F*)dir->FindObject("hPLJetPz"))->Fill(obj.at(i).Pz, eventWeight);
      ((TH1F*)dir->FindObject("hPLJetArea"))->Fill(obj.at(i).Area, eventWeight);
      ((TH1F*)dir->FindObject("hPLJetNPart"))->Fill(obj.at(i).NPart, eventWeight);
      ((TH2F*)dir->FindObject("hPLJetEtaPhi"))->Fill(obj.at(i).Eta(), obj.at(i).Phi(), eventWeight);
      if (obj.at(i).ClosestDLJet != nullptr) {
        ((TH2F*)dir->FindObject("hDLJetPVsPLJetP"))->Fill(obj.at(i).ClosestDLJet->P(), obj.at(i).P(), eventWeight);
        ((TH2F*)dir->FindObject("hDLJetPxVsPLJetPx"))->Fill(fabs(obj.at(i).ClosestDLJet->Px), fabs(obj.at(i).Px), eventWeight);
        ((TH2F*)dir->FindObject("hDLJetPyVsPLJetPy"))->Fill(fabs(obj.at(i).ClosestDLJet->Py), fabs(obj.at(i).Py), eventWeight);
      }
    }
  }
  if constexpr (std::is_same<T, Event>::value) {
    ((TH1F*)dir->FindObject("hZVtx"))->Fill(obj.ZVtx, eventWeight);
    if (dir->FindObject("hWeights")) {
      ((TH1F*)dir->FindObject("hWeights"))->Fill(obj.weight);
    }
  }
}