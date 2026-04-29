#include "Utilities.h"
#include "Binnings.h"
#include "TMath.h"
#include <filesystem>
#include <yaml-cpp/yaml.h>
#include "THnSparse.h"

// ============================================================
// ABCD cut boundaries (loaded from Cuts.yaml)
// ============================================================
struct ABCDCuts {
  double isoSigMax;  // axis 0: upper edge of signal ISO region
  double isoSbMin;   // axis 0: lower edge of ISO sideband
  double isoSbMax;   // axis 0: upper edge of ISO sideband
  double m02SigMin;  // axis 1: lower edge of signal M02 window
  double m02SigMax;  // axis 1: upper edge of signal M02 window
  double m02SbMin;   // axis 1: lower edge of background M02 sideband
  double m02SbMax;   // axis 1: upper edge of background M02 sideband
};



// ============================================================
// Load ABCD cuts from Cuts.yaml (inheriting from Standard)
// ============================================================
ABCDCuts loadABCDCuts(TString cutString)
{
  YAML::Node ycut    = YAML::LoadFile("Cuts.yaml");
  if (!ycut[(std::string)cutString])
    FATAL(Form("Cut string '%s' not found in Cuts.yaml", cutString.Data()))

  YAML::Node standard = ycut["Standard"];
  YAML::Node chosen   = ycut[(std::string)cutString];

  auto get = [&](const char *key) -> double {
    return chosen[key].IsDefined() ? chosen[key].as<double>()
                                   : standard[key].as<double>();
  };

  ABCDCuts c;
  c.isoSigMax = get("purity_iso_sig_max");
  c.isoSbMin  = get("purity_iso_sb_min");
  c.isoSbMax  = get("purity_iso_sb_max");
  c.m02SigMin = get("purity_m02_sig_min");
  c.m02SigMax = get("purity_m02_sig_max");
  c.m02SbMin  = get("purity_m02_sb_min");
  c.m02SbMax  = get("purity_m02_sb_max");
  return c;
}

// ============================================================
// ABCD region helpers
// ============================================================
// The hPurityABCD THnSparse axes are:
//   0: IsoCharged,  1: M02,  2: photonPt,  3: jetPt
// The hPurityABCDRg/Zg/Nsd add:
//   4: substructure observable
//
// Region definitions:
//              | M02 signal window  | M02 sideband window
// -------------+--------------------+--------------------
// ISO signal   |   A (signal-rich)  |   B
// ISO sideband |   C                |   D
//
// Purity = 1 - (N_B * N_C) / (N_A * N_D)

void setABCDRegion(THnSparseF *h, const ABCDCuts &cuts, char region)
{
  double isoMin = h->GetAxis(0)->GetXmin();
  double isoMax, m02Min, m02Max;
  switch (region) {
    case 'A': isoMax = cuts.isoSigMax; m02Min = cuts.m02SigMin; m02Max = cuts.m02SigMax; break;
    case 'B': isoMax = cuts.isoSigMax; m02Min = cuts.m02SbMin; m02Max = cuts.m02SbMax; break;
    case 'C': isoMin = cuts.isoSbMin; isoMax = cuts.isoSbMax; m02Min = cuts.m02SigMin; m02Max = cuts.m02SigMax; break;
    case 'D': isoMin = cuts.isoSbMin; isoMax = cuts.isoSbMax; m02Min = cuts.m02SbMin; m02Max = cuts.m02SbMax; break;
    default:  FATAL("setABCDRegion: unknown region (must be A, B, C, or D)")
  }
  h->GetAxis(0)->SetRangeUser(isoMin, isoMax);
  h->GetAxis(1)->SetRangeUser(m02Min, m02Max);
}

// Reset axes 0 and 1 to their full range (removes the ABCD restriction)
void resetABCDAxes(THnSparseF *h)
{
  h->GetAxis(0)->SetRange(0, 0);
  h->GetAxis(1)->SetRange(0, 0);
}

// ============================================================
// Purity computation from ABCD region histograms
// ============================================================
// Bin-by-bin: purity[i] = 1 - ((C/A) / (D/B))

// Core bin-by-bin ABCD purity on histograms with any binning (no rebinning)
static void fillPurityBinByBin(TH1D *hPurity,
                                const TH1D *hA, const TH1D *hB,
                                const TH1D *hC, const TH1D *hD)
{
  for (int i = 1; i <= hPurity->GetNbinsX(); i++) {
    double A = hA->GetBinContent(i), Aerr = hA->GetBinError(i);
    double B = hB->GetBinContent(i), Berr = hB->GetBinError(i);
    double C = hC->GetBinContent(i), Cerr = hC->GetBinError(i);
    double D = hD->GetBinContent(i), Derr = hD->GetBinError(i);

    if (A * D <= 0.)
      continue;

    double bgFrac = (B * C) / (A * D);
    double bgErr  = TMath::Sqrt(TMath::Power(Cerr * (B / (A * D)),           2) +
                                TMath::Power(Berr * (C / (A * D)),           2) +
                                TMath::Power(Aerr * (B * C) / (A * A * D),   2) +
                                TMath::Power(Derr * (B * C) / (A * D * D),   2));
    hPurity->SetBinContent(i, 1. - bgFrac);
    hPurity->SetBinError(i, bgErr);
  }
}

// Purity vs photon pT: rebin to ptBinsTrigger then compute
TH1D *purityFromABCD1D(TH1D *hA, TH1D *hB, TH1D *hC, TH1D *hD, TString name)
{
  TH1D *hAr = (TH1D *)hA->Rebin(nPtBinsTrigger, "_tmp_Ar", ptBinsTrigger);
  TH1D *hBr = (TH1D *)hB->Rebin(nPtBinsTrigger, "_tmp_Br", ptBinsTrigger);
  TH1D *hCr = (TH1D *)hC->Rebin(nPtBinsTrigger, "_tmp_Cr", ptBinsTrigger);
  TH1D *hDr = (TH1D *)hD->Rebin(nPtBinsTrigger, "_tmp_Dr", ptBinsTrigger);

  TH1D *hPurity = new TH1D(name.Data(), name.Data(), nPtBinsTrigger, ptBinsTrigger);
  hPurity->GetYaxis()->SetTitle("Purity");
  fillPurityBinByBin(hPurity, hAr, hBr, hCr, hDr);

  delete hAr;
  delete hBr;
  delete hCr;
  delete hDr;
  return hPurity;
}

// Purity vs substructure observable: use the native projection binning (no rebin)
TH1D *purityFromABCD1D_substructure(TH1D *hA, TH1D *hB, TH1D *hC, TH1D *hD, TString name)
{
  TH1D *hPurity = (TH1D *)hA->Clone(name.Data());
  hPurity->Reset();
  hPurity->GetYaxis()->SetTitle("Purity");
  fillPurityBinByBin(hPurity, hA, hB, hC, hD);
  return hPurity;
}

// ROOT has no variable-bin TH2::Rebin; accumulate source bins into target bins manually
static TH2D *rebin2D(const TH2D *h, int nx, const Double_t *xbins,
                                    int ny, const Double_t *ybins,
                                    const char *newName)
{
  TH2D *hr = new TH2D(newName, h->GetTitle(), nx, xbins, ny, ybins);
  hr->Sumw2();
  for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h->GetNbinsY(); iy++) {
      double val = h->GetBinContent(ix, iy);
      if (val == 0.) continue;
      double x  = h->GetXaxis()->GetBinCenter(ix);
      double y  = h->GetYaxis()->GetBinCenter(iy);
      int    jx = hr->GetXaxis()->FindBin(x);
      int    jy = hr->GetYaxis()->FindBin(y);
      if (jx < 1 || jx > nx || jy < 1 || jy > ny) continue;
      hr->SetBinContent(jx, jy, hr->GetBinContent(jx, jy) + val);
      hr->SetBinError(jx, jy, TMath::Sqrt(TMath::Power(hr->GetBinError(jx, jy), 2) +
                                           TMath::Power(h->GetBinError(ix, iy),   2)));
    }
  }
  return hr;
}

TH2D *purityFromABCD2D(TH2D *hA, TH2D *hB, TH2D *hC, TH2D *hD, TString name)
{
  TH2D *hAr = rebin2D(hA, nPtBinsTrigger, ptBinsTrigger, nPtBinsJetCoarse, ptBinsJetCoarse, "_tmp_Ar");
  TH2D *hBr = rebin2D(hB, nPtBinsTrigger, ptBinsTrigger, nPtBinsJetCoarse, ptBinsJetCoarse, "_tmp_Br");
  TH2D *hCr = rebin2D(hC, nPtBinsTrigger, ptBinsTrigger, nPtBinsJetCoarse, ptBinsJetCoarse, "_tmp_Cr");
  TH2D *hDr = rebin2D(hD, nPtBinsTrigger, ptBinsTrigger, nPtBinsJetCoarse, ptBinsJetCoarse, "_tmp_Dr");

  TH2D *hPurity = new TH2D(name.Data(), name.Data(), nPtBinsTrigger, ptBinsTrigger, nPtBinsJetCoarse, ptBinsJetCoarse);
  hPurity->GetZaxis()->SetTitle("Purity");

  for (int ix = 1; ix <= hPurity->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= hPurity->GetNbinsY(); iy++) {
      double A = hAr->GetBinContent(ix, iy), Aerr = hAr->GetBinError(ix, iy);
      double B = hBr->GetBinContent(ix, iy), Berr = hBr->GetBinError(ix, iy);
      double C = hCr->GetBinContent(ix, iy), Cerr = hCr->GetBinError(ix, iy);
      double D = hDr->GetBinContent(ix, iy), Derr = hDr->GetBinError(ix, iy);

      if (A * D <= 0.)
        continue;

      double bgFrac = (B * C) / (A * D);
      double bgErr  = TMath::Sqrt(TMath::Power(Cerr * (B / (A * D)),           2) +
                                  TMath::Power(Berr * (C / (A * D)),           2) +
                                  TMath::Power(Aerr * (B * C) / (A * A * D),   2) +
                                  TMath::Power(Derr * (B * C) / (A * D * D),   2));
      hPurity->SetBinContent(ix, iy, 1. - bgFrac);
      hPurity->SetBinError(ix, iy, bgErr);
    }
  }

  delete hAr;
  delete hBr;
  delete hCr;
  delete hDr;
  return hPurity;
}

// ============================================================
// Level 1: purity vs photon pT  (integrated over jet pT)
// Uses hPurityABCD (4D).
// ============================================================
TH1D *computeRawPurity1D(THnSparseF *h, const ABCDCuts &cuts,
                         TDirectory *outDir, TString name)
{
  h->GetAxis(3)->SetRange(0, 0); // integrate over all jetPt

  setABCDRegion(h, cuts, 'A');
  TH1D *hA = (TH1D *)h->Projection(2); hA->SetName(Form("%s_A", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'B');
  TH1D *hB = (TH1D *)h->Projection(2); hB->SetName(Form("%s_B", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'C');
  TH1D *hC = (TH1D *)h->Projection(2); hC->SetName(Form("%s_C", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'D');
  TH1D *hD = (TH1D *)h->Projection(2); hD->SetName(Form("%s_D", name.Data()));
  resetABCDAxes(h);

  TH1D *hPurity = purityFromABCD1D(hA, hB, hC, hD, name);

  if (outDir) {
    outDir->Add(hA);
    outDir->Add(hB);
    outDir->Add(hC);
    outDir->Add(hD);
    outDir->Add(hPurity);
  }
  return hPurity;
}

// ============================================================
// Level 2: purity vs (photon pT, jet pT)
// Uses hPurityABCD (4D).
// The returned TH2D has photon pT on x and jet pT on y.
// ============================================================
TH2D *computeRawPurity2D(THnSparseF *h, const ABCDCuts &cuts,
                         TDirectory *outDir, TString name)
{
  h->GetAxis(3)->SetRange(0, 0);

  // Projection(yDim, xDim): x = photonPt (axis 2), y = jetPt (axis 3)
  setABCDRegion(h, cuts, 'A');
  TH2D *hA = (TH2D *)h->Projection(3, 2); hA->SetName(Form("%s_A", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'B');
  TH2D *hB = (TH2D *)h->Projection(3, 2); hB->SetName(Form("%s_B", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'C');
  TH2D *hC = (TH2D *)h->Projection(3, 2); hC->SetName(Form("%s_C", name.Data()));
  resetABCDAxes(h);

  setABCDRegion(h, cuts, 'D');
  TH2D *hD = (TH2D *)h->Projection(3, 2); hD->SetName(Form("%s_D", name.Data()));
  resetABCDAxes(h);

  TH2D *hPurity = purityFromABCD2D(hA, hB, hC, hD, name);

  if (outDir) {
    outDir->Add(hA);
    outDir->Add(hB);
    outDir->Add(hC);
    outDir->Add(hD);
    outDir->Add(hPurity);
  }
  return hPurity;
}

// ============================================================
// Level 3: purity vs substructure observable
//          for a single (photon pT, jet pT) slice
// Uses hPurityABCDRg / hPurityABCDZg / hPurityABCDNsd (5D).
// Axes: 0=ISO, 1=M02, 2=photonPt, 3=jetPt, 4=observable
// ============================================================
TH1D *computeRawPurity1D_substructure(THnSparseF *h5D, const ABCDCuts &cuts,
                                      double photonPtMin, double photonPtMax,
                                      double jetPtMin,   double jetPtMax,
                                      TString name)
{
  h5D->GetAxis(2)->SetRangeUser(photonPtMin, photonPtMax);
  h5D->GetAxis(3)->SetRangeUser(jetPtMin,    jetPtMax);

  setABCDRegion(h5D, cuts, 'A');
  TH1D *hA = (TH1D *)h5D->Projection(4); hA->SetName(Form("%s_A", name.Data()));
  resetABCDAxes(h5D);

  setABCDRegion(h5D, cuts, 'B');
  TH1D *hB = (TH1D *)h5D->Projection(4); hB->SetName(Form("%s_B", name.Data()));
  resetABCDAxes(h5D);

  setABCDRegion(h5D, cuts, 'C');
  TH1D *hC = (TH1D *)h5D->Projection(4); hC->SetName(Form("%s_C", name.Data()));
  resetABCDAxes(h5D);

  setABCDRegion(h5D, cuts, 'D');
  TH1D *hD = (TH1D *)h5D->Projection(4); hD->SetName(Form("%s_D", name.Data()));
  resetABCDAxes(h5D);

  h5D->GetAxis(2)->SetRange(0, 0);
  h5D->GetAxis(3)->SetRange(0, 0);

  TH1D *hPurity = purityFromABCD1D_substructure(hA, hB, hC, hD, name);

  delete hA;
  delete hB;
  delete hC;
  delete hD;
  return hPurity;
}

// ============================================================
// MC Correction Factor (stub — implement when MC is available)
// ============================================================
// Returns a per-pT correction histogram to be multiplied with the raw purity.
// Currently returns nullptr, which signals "no correction" to computeCorrectedPurity.
//
// TODO: Implement when MC signal (e.g. GJ) and background (e.g. JJ) samples are available.
//
// Algorithm outline:
//   1. From MC signal sample: compute the true signal fraction in each ABCD region
//      using generator-level truth labels (is this cluster from a direct photon?).
//   2. From MC background sample: compute the background contamination in region A
//      and leakage of signal into regions B, C, D.
//   3. Derive a per-pT correction factor:
//      corrFactor(pT) = purityTrue_MC(pT) / purityABCD_MC(pT)
//   4. Apply corrFactor to the data raw purity.
//
// Expected inputs:
//   hSigMC : hPurityABCD (4D) from MC signal tune (e.g. GJ Pythia)
//   hBkgMC : hPurityABCD (4D) from MC background tune (e.g. JJ Pythia)

TH1D *computeMCCorrectionFactor1D(THnSparseF * /*hSigMC*/, THnSparseF * /*hBkgMC*/,
                                  const ABCDCuts & /*cuts*/)
{
  return nullptr;
}

TH2D *computeMCCorrectionFactor2D(THnSparseF * /*hSigMC*/, THnSparseF * /*hBkgMC*/,
                                  const ABCDCuts & /*cuts*/)
{
  return nullptr;
}

TH1D *computeMCCorrectionFactor1D_substructure(THnSparseF * /*hSigMC5D*/, THnSparseF * /*hBkgMC5D*/,
                                               const ABCDCuts & /*cuts*/,
                                               double /*photonPtMin*/, double /*photonPtMax*/,
                                               double /*jetPtMin*/,    double /*jetPtMax*/)
{
  return nullptr;
}

// ============================================================
// Apply MC correction. Pass-through when mcCorr == nullptr.
// ============================================================
TH1D *computeCorrectedPurity(TH1D *rawPurity, TH1D *mcCorr)
{
  if (!mcCorr)
    return rawPurity;
  TH1D *corrected = (TH1D *)rawPurity->Clone(Form("%s_MCcorrected", rawPurity->GetName()));
  corrected->Multiply(mcCorr);
  return corrected;
}

TH2D *computeCorrectedPurity(TH2D *rawPurity, TH2D *mcCorr)
{
  if (!mcCorr)
    return rawPurity;
  TH2D *corrected = (TH2D *)rawPurity->Clone(Form("%s_MCcorrected", rawPurity->GetName()));
  corrected->Multiply(mcCorr);
  return corrected;
}

// ============================================================
// Main entry point
// ============================================================
void calculatePurity(TString AnalysisDirectory, bool isDebugRun = false,
                     TString configPath = "RunConfig.yaml")
{
  ENTER

  // Prevent ROOT from auto-registering histograms into gDirectory on construction.
  // Without this, temporary histograms (e.g. Rebin intermediates) get added to the
  // output TFile's directory list and are double-deleted when TFile::Close() runs.
  TH1::AddDirectory(false);

  GlobalOptions optns(AnalysisDirectory, isDebugRun ? 0 : 1, configPath);


  ABCDCuts cuts = loadABCDCuts(optns.cutString);
  INFO(Form("ABCD cuts: ISO signal < %.1f | ISO sideband [%.1f, %.1f] | M02 signal [%.2f, %.2f] | M02 sideband [%.2f, %.2f]",
    cuts.isoSigMax, cuts.isoSbMin, cuts.isoSbMax,
    cuts.m02SigMin, cuts.m02SigMax, cuts.m02SbMin, cuts.m02SbMax))

  // ---- Open input file ----
  TString inputPath = Form("%s/HistosFromTree.root", AnalysisDirectory.Data());
  if (!std::filesystem::exists(inputPath.Data()))
    inputPath = Form("%s/HistosFromTree_0.root", AnalysisDirectory.Data());

  TFile *fIn = TFile::Open(inputPath, "READ");
  if (!fIn || fIn->IsZombie())
    FATAL(Form("Input file not found: %s", inputPath.Data()))

  TDirectory *dPurity = (TDirectory *)fIn->Get("PurityHistograms");
  if (!dPurity)
    FATAL("PurityHistograms directory not found in input file")

  THnSparseF *hABCD = (THnSparseF *)dPurity->Get("hPurityABCD");
  if (!hABCD)
    FATAL("hPurityABCD not found in PurityHistograms")

  THnSparseF *hABCDRg  = optns.doSubstructure ? (THnSparseF *)dPurity->Get("hPurityABCDRg")  : nullptr;
  THnSparseF *hABCDZg  = optns.doSubstructure ? (THnSparseF *)dPurity->Get("hPurityABCDZg")  : nullptr;
  THnSparseF *hABCDNsd = optns.doSubstructure ? (THnSparseF *)dPurity->Get("hPurityABCDNsd") : nullptr;

  // ---- Create output file ----
  TFile *fOut = TFile::Open(Form("%s/Purity.root", AnalysisDirectory.Data()), "RECREATE");

  // ============================================================
  // Level 1: purity vs photon pT  (integrated over all jet pT)
  // ============================================================
  INFO("Computing Level 1: purity vs photon pT")
  TDirectory *dL1 = fOut->mkdir("PurityVsPhotonPt");

  TH1D *hRawPurity1D = computeRawPurity1D(hABCD, cuts, dL1, "hPurity_vs_photonPt");
  // Placeholder for MC correction (currently a no-op)
  TH1D *mcCorr1D = computeMCCorrectionFactor1D(nullptr, nullptr, cuts);
  TH1D *hPurity1D = computeCorrectedPurity(hRawPurity1D, mcCorr1D);
  if (hPurity1D != hRawPurity1D)
    dL1->Add(hPurity1D);

  // ============================================================
  // Level 2: purity vs (photon pT, jet pT)
  // ============================================================
  INFO("Computing Level 2: purity vs photon pT and jet pT")
  TDirectory *dL2 = fOut->mkdir("PurityVsPhotonPtJetPt");

  TH2D *hRawPurity2D = computeRawPurity2D(hABCD, cuts, dL2, "hPurity_vs_photonPt_jetPt");

  TH2D *mcCorr2D = computeMCCorrectionFactor2D(nullptr, nullptr, cuts);
  TH2D *hPurity2D = computeCorrectedPurity(hRawPurity2D, mcCorr2D);
  if (hPurity2D != hRawPurity2D)
    dL2->Add(hPurity2D);

  // ============================================================
  // Level 3: purity vs substructure observable
  //          in slices of (photon pT, jet pT)
  // ============================================================
  if (optns.doSubstructure) {
    if (!hABCDRg || !hABCDZg || !hABCDNsd) {
      INFO("doSubstructure is enabled but one or more 5D purity histograms are missing — skipping Level 3.")
    } else {
      INFO("Computing Level 3: purity vs substructure observables in (photon pT, jet pT) slices")

      TDirectory *dRg  = fOut->mkdir("PurityVsRg");
      TDirectory *dZg  = fOut->mkdir("PurityVsZg");
      TDirectory *dNsd = fOut->mkdir("PurityVsNsd");

      for (int ig = 0; ig < (int)kPhotonPtBins.size() - 1; ig++) {
        double gPtLo = kPhotonPtBins[ig];
        double gPtHi = kPhotonPtBins[ig + 1];
        for (int ij = 0; ij < (int)kJetPtBins.size() - 1; ij++) {
          double jPtLo = kJetPtBins[ij];
          double jPtHi = kJetPtBins[ij + 1];

          TString sliceLabel = Form("photonPt_%.0f_%.0f_jetPt_%.0f_%.0f",
                                    gPtLo, gPtHi, jPtLo, jPtHi);

          TH1D *hRg = computeRawPurity1D_substructure(
            hABCDRg, cuts, gPtLo, gPtHi, jPtLo, jPtHi,
            Form("hPurity_Rg_%s", sliceLabel.Data()));
          TH1D *mcCorrRg = computeMCCorrectionFactor1D_substructure(
            nullptr, nullptr, cuts, gPtLo, gPtHi, jPtLo, jPtHi);
          TH1D *hRgCorr = computeCorrectedPurity(hRg, mcCorrRg);
          dRg->Add(hRgCorr);
          if (hRgCorr != hRg) delete hRg;

          TH1D *hZg = computeRawPurity1D_substructure(
            hABCDZg, cuts, gPtLo, gPtHi, jPtLo, jPtHi,
            Form("hPurity_Zg_%s", sliceLabel.Data()));
          TH1D *mcCorrZg = computeMCCorrectionFactor1D_substructure(
            nullptr, nullptr, cuts, gPtLo, gPtHi, jPtLo, jPtHi);
          TH1D *hZgCorr = computeCorrectedPurity(hZg, mcCorrZg);
          dZg->Add(hZgCorr);
          if (hZgCorr != hZg) delete hZg;

          TH1D *hNsd = computeRawPurity1D_substructure(
            hABCDNsd, cuts, gPtLo, gPtHi, jPtLo, jPtHi,
            Form("hPurity_Nsd_%s", sliceLabel.Data()));
          TH1D *mcCorrNsd = computeMCCorrectionFactor1D_substructure(
            nullptr, nullptr, cuts, gPtLo, gPtHi, jPtLo, jPtHi);
          TH1D *hNsdCorr = computeCorrectedPurity(hNsd, mcCorrNsd);
          dNsd->Add(hNsdCorr);
          if (hNsdCorr != hNsd) delete hNsd;
        }
      }
    }
  }

  fOut->Write();
  fIn->Close();
  fOut->Close();

  EXIT
}
