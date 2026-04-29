#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "Binnings.h"
#include "THnSparse.h"

const char *suffix = "pdf";

Int_t GoetheColor[4][2];

TString legendHeader;
TString legendHeaderJets;
void BuildLegendHeader(GlobalOptions optns)
{
  // load relevant cuts
  EventCuts eventCuts(optns);
  YAML::Node ycut = YAML::LoadFile("Cuts.yaml");
  if (!ycut[(std::string)optns.cutString])
    FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))
  YAML::Node standardcut = ycut["Standard"];
  YAML::Node chosencut = ycut[(std::string)optns.cutString];
  double EMin = chosencut["gamma_min_E"].IsDefined() ? chosencut["gamma_min_E"].as<float>() : standardcut["gamma_min_E"].as<float>();
  
  TString trainConfig = optns.trainConfig;
  // extract centrality min and centrality max. Trainconfig has always form Run3-$centMin-$centMax
  // Tokenize the string
  std::vector<std::string> tokens;
  std::stringstream ss(trainConfig.Data());
  std::string token;
  while (std::getline(ss, token, '-')) {
    tokens.push_back(token);
  }
  int centMin = 0;
  int centMax = 100;
  if(tokens.size() == 3)
  {
    centMin = std::stoi(tokens[1]);
    centMax = std::stoi(tokens[2]);
  }
  if (optns.isPP)
    legendHeader = Form("#bf{ALICE work-in-progress};%s", optns.dataSet.Data());
  else
    legendHeader = Form("#bf{ALICE work-in-progress};%s;Centrality: %d-%d%%", optns.dataSet.Data(), centMin, centMax);

  // get jet radius
  double jetRadius = chosencut["jet_Radius"].IsDefined() ? chosencut["jet_Radius"].as<float>() : standardcut["jet_Radius"].as<float>();
  TString jetAcceptance = chosencut["jet_Acceptance"].IsDefined() ? chosencut["jet_Acceptance"].as<std::string>() : standardcut["jet_Acceptance"].as<std::string>();
  legendHeaderJets = Form("%s;Jet Radius: %.2f;Jet Acceptance: %s",legendHeader.Data(),jetRadius,jetAcceptance.Data());
}

void plotPurityInclusive(TDirectory *dPurity, TString dirname, GlobalOptions optns)
{
  ENTER

  TString outputDir = Form("%s/%s", optns.analysisDirPath.Data(), dirname.Data());
  createDirectory(outputDir);

  // from dPurity, get hPurity_vs_photonPt
  TH1F* hPurity_vs_photonPt = (TH1F*)dPurity->Get("hPurity_vs_photonPt");
  if (!hPurity_vs_photonPt)
    FATAL("hPurity_vs_photonPt not found in Purity")

  // plot hPurity_vs_photonPt
  Plotting1D PurityVsPhotonPt;
  PurityVsPhotonPt.New(hPurity_vs_photonPt);
  PurityVsPhotonPt.SetAxisLabel("#it{p}_{T}^{#gamma} (GeV/c)", "#bf{Purity}");
  PurityVsPhotonPt.NewLatex(0.9,0.4,legendHeaderJets);
  PurityVsPhotonPt.Plot(Form("%s/PurityVsPhotonPt.%s", outputDir.Data(), suffix));  
  EXIT
}

void plotPurityVsPhotonPtJetPt(TDirectory *dPurity, TString dirname, GlobalOptions optns)
{
  ENTER

  TString outputDir = Form("%s/%s", optns.analysisDirPath.Data(), dirname.Data());
  createDirectory(outputDir);

  TH2F* hPurity_vs_PhotonPt_jetPt = (TH2F*)dPurity->Get("hPurity_vs_photonPt_jetPt");
  if (!hPurity_vs_PhotonPt_jetPt)
    FATAL("hPurity_vs_PhotonPt_jetPt not found in PurityVsPhotonPtJetPt")

  gStyle->SetPaintTextFormat(".2f");
  hPurity_vs_PhotonPt_jetPt->SetMarkerColor(kRed);

  Plotting2D PurityVsPhotonPtJetPt;
  PurityVsPhotonPtJetPt.SetMargins(0.12, 0.1, 0.025, 0.17);
  PurityVsPhotonPtJetPt.New(hPurity_vs_PhotonPt_jetPt, "COLZ TEXT");
  PurityVsPhotonPtJetPt.SetAxisRange(hPurity_vs_PhotonPt_jetPt->GetXaxis()->GetXmin(), hPurity_vs_PhotonPt_jetPt->GetXaxis()->GetXmax(),
                                     hPurity_vs_PhotonPt_jetPt->GetYaxis()->GetXmin(), hPurity_vs_PhotonPt_jetPt->GetYaxis()->GetXmax());
  PurityVsPhotonPtJetPt.SetAxisLabel("#it{p}_{T}^{#gamma} (GeV/c)", "#it{p}_{T}^{jet} (GeV/c)");
  // PurityVsPhotonPtJetPt.NewLatex(0.9, 0.4, legendHeaderJets);
  PurityVsPhotonPtJetPt.Plot(Form("%s/PurityVsPhotonPtJetPt.%s", outputDir.Data(), suffix));

  EXIT
}

void plotPuritySubstructure(TDirectory* dRg, TDirectory* dZg, TDirectory* dNsd, TString dirname, GlobalOptions optns)
{
  ENTER

  TH1::AddDirectory(false);

  TString outputDir = Form("%s/%s/SubstructurePurity", optns.analysisDirPath.Data(), dirname.Data());
  createDirectory(outputDir);

  int nGammaBins = (int)kPhotonPtBins.size() - 1;
  int nJetBins   = (int)kJetPtBins.size()    - 1;

  const int kNObs = 3;
  TDirectory* obsDirs[kNObs]    = {dRg,               dZg,               dNsd};
  TString     obsPrefixes[kNObs] = {"Rg",             "Zg",              "Nsd"};
  TString     obsXLabels[kNObs]  = {"#it{R}_{g}",     "#it{z}_{g}",      "#it{n}_{SD}"};
  TString     obsNames[kNObs]    = {"Rg",             "Zg",              "Nsd"};
  Double_t    xmin[kNObs] = {0., 0., 0.};
  Double_t    xmax[kNObs] = {0.6, 0.6, 15.};
  for (int iObs = 0; iObs < kNObs; iObs++) {
    if (!obsDirs[iObs]) {
      INFO(Form("Directory for observable %s not available, skipping", obsNames[iObs].Data()))
      continue;
    }

    TString obsOutputDir = Form("%s/%s", outputDir.Data(), obsNames[iObs].Data());
    createDirectory(obsOutputDir);

    // Grid canvas: columns = photon pT bins, rows = jet pT bins
    TCanvas* cGrid = new TCanvas(
      Form("cGrid_Purity_%s", obsNames[iObs].Data()),
      Form("Purity vs %s", obsNames[iObs].Data()),
      300 * nGammaBins, 300 * nJetBins
    );
    cGrid->Divide(nGammaBins, nJetBins, 0, 0);

    for (int ig = 0; ig < nGammaBins; ig++) {
      double gPtLo = kPhotonPtBins[ig];
      double gPtHi = kPhotonPtBins[ig + 1];

      for (int ij = 0; ij < nJetBins; ij++) {
        double jPtLo = kJetPtBins[ij];
        double jPtHi = kJetPtBins[ij + 1];

        TString sliceLabel = Form("photonPt_%.0f_%.0f_jetPt_%.0f_%.0f",
                                   gPtLo, gPtHi, jPtLo, jPtHi);
        TString histName   = Form("hPurity_%s_%s", obsPrefixes[iObs].Data(), sliceLabel.Data());

        TH1D* hPurity = (TH1D*)obsDirs[iObs]->Get(histName);
        if (!hPurity) {
          INFO(Form("Histogram %s not found, skipping", histName.Data()))
          continue;
        }

        // ---- Individual 1D plot ----
          TString legLabel = Form("%.0f < #it{p}_{T}^{#gamma} < %.0f GeV/#it{c};%.0f < #it{p}_{T}^{jet} < %.0f GeV/#it{c}",
                                  gPtLo, gPtHi, jPtLo, jPtHi);
        Plotting1D pPlot;
        pPlot.New(hPurity);
        pPlot.SetAxisRange(xmin[iObs], xmax[iObs], 0., 1.1);
        pPlot.SetAxisLabel(obsXLabels[iObs], "#bf{Purity}");
        pPlot.SetLegend(0.15, 0.62, 0.72, 0.92);
        TString legendText = Form("%s;%.0f < #it{p}_{T}^{#gamma} < %.0f GeV/#it{c};%.0f < #it{p}_{T}^{jet} < %.0f GeV/#it{c}",
                                  legendHeaderJets.Data(), gPtLo, gPtHi, jPtLo, jPtHi);
        pPlot.NewLatex(0.9, 0.9, legendText);
        pPlot.Plot(Form("%s/Purity_%s_%s.%s",
                        obsOutputDir.Data(), obsNames[iObs].Data(), sliceLabel.Data(), suffix));

        // ---- Grid pad ----
        // rows = jet pT (top row = lowest jet pT), cols = photon pT (leftmost = lowest)
        int padIdx = ij * nGammaBins + ig + 1;
        cGrid->cd(padIdx);
        gPad->SetTopMargin(0.14);
        gPad->SetBottomMargin(0.14);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.03);

        TH1D* hDraw = (TH1D*)hPurity->Clone(Form("hDraw_%s_%s", obsNames[iObs].Data(), sliceLabel.Data()));
        hDraw->SetTitle(Form("#it{p}_{T}^{#gamma}: %.0f-%.0f  #it{p}_{T}^{jet}: %.0f-%.0f GeV/#it{c}",
                              gPtLo, gPtHi, jPtLo, jPtHi));
        hDraw->GetXaxis()->SetTitle(obsXLabels[iObs]);
        hDraw->GetXaxis()->SetTitleSize(0.07);
        hDraw->GetXaxis()->SetLabelSize(0.06);
        hDraw->GetXaxis()->SetTitleOffset(0.9);
        hDraw->GetYaxis()->SetTitle("Purity");
        hDraw->GetYaxis()->SetTitleSize(0.07);
        hDraw->GetYaxis()->SetLabelSize(0.06);
        hDraw->GetYaxis()->SetTitleOffset(1.0);
        hDraw->GetXaxis()->SetRangeUser(xmin[iObs], xmax[iObs]);
        hDraw->GetYaxis()->SetRangeUser(0., 1.1);
        hDraw->SetMarkerStyle(20);
        hDraw->SetMarkerSize(0.8);
        hDraw->SetMarkerColor(kBlack);
        hDraw->SetLineColor(kBlack);
        hDraw->Draw("e p");
      }
    }

    cGrid->SaveAs(Form("%s/PurityGrid_%s.%s", outputDir.Data(), obsNames[iObs].Data(), suffix));
    delete cGrid;
  }

  EXIT
}

void plotPurity(TString AnalysisDirectory, bool isDebugRun = false, TString configPath = "RunConfig.yaml")
{
  ENTER

  GlobalOptions optns(AnalysisDirectory, isDebugRun ? 0 : 1, configPath);

  BuildLegendHeader(optns);

  TString inputFilePath = Form("%s/Purity.root", AnalysisDirectory.Data());

  TFile *fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  TDirectory *dPurityVsPhotonPt = (TDirectory *)fIn->Get("PurityVsPhotonPt");
  if (!dPurityVsPhotonPt)
    FATAL("PurityVsPhotonPt directory not found in input file")

  plotPurityInclusive(dPurityVsPhotonPt, "Purity", optns);

  TDirectory *dPurityVsPhotonPtJetPt = (TDirectory *)fIn->Get("PurityVsPhotonPtJetPt");
  if (!dPurityVsPhotonPtJetPt)
    FATAL("PurityVsPhotonPtJetPt directory not found in input file")

  plotPurityVsPhotonPtJetPt(dPurityVsPhotonPtJetPt, "Purity", optns);

  TDirectory *dPurityVsRg  = (TDirectory *)fIn->Get("PurityVsRg");
  TDirectory *dPurityVsZg  = (TDirectory *)fIn->Get("PurityVsZg");
  TDirectory *dPurityVsNsd = (TDirectory *)fIn->Get("PurityVsNsd");

  if (dPurityVsRg || dPurityVsZg || dPurityVsNsd) {
    plotPuritySubstructure(dPurityVsRg, dPurityVsZg, dPurityVsNsd, "Purity", optns);
  } else {
    INFO("No substructure purity directories found in input file — skipping plotPuritySubstructure")
  }

  EXIT
}