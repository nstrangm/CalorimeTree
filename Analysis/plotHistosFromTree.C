#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "plotHistosFromTree.h"

const char *suffix = "pdf";

Int_t GoetheColor[4][2];

// Add to plot class?
void AvoidLogYconflict(TH1 *h, float logymin = 1e-10)
{
  double min_nonzero = std::numeric_limits<double>::max();
  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    float content = h->GetBinContent(i);
    if (content > 1e-10 && content < min_nonzero)
    {
      min_nonzero = content;
    }
    if (h->GetBinContent(i) == 0)
    {
      h->SetBinContent(i, logymin);
    }
  }
  h->SetMinimum(min_nonzero);
  return;
}

void SetTitleAndLabels(TH1 *h, const char *title)
{
  h->SetTitle(title);
  return;
}

void plotEventQA(TDirectory* dEventQA, GlobalOptions optns)
{
  ENTER

  const char *outputDir = Form("%s/EventQA", optns.analysisDirPath.Data());
  createDirectory(outputDir);

  THnSparseF *hCentralityRhoMultiplicity = (THnSparseF *)dEventQA->Get("hCentralityRhoMultiplicity");
  TH1F* hCentrality = (TH1F*)hCentralityRhoMultiplicity->Projection(0);
  TH1F* hRho = (TH1F*)hCentralityRhoMultiplicity->Projection(1);
  TH1F* hMultiplicity = (TH1F*)hCentralityRhoMultiplicity->Projection(2);



  Plotting1D PCentrality;
  // set log y
  AvoidLogYconflict(hCentrality);
  

  PCentrality.New(hCentrality);
  PCentrality.SetAxisLabel("#bf{Centrality}", "#bf{Counts}");
  PCentrality.Plot(Form("%s/Centrality.%s", outputDir, suffix),false, true);

  Plotting1D PRho;
  PRho.New(hRho);
  PRho.SetAxisLabel("#bf{#rho} (GeV/#it{c})", "#bf{Counts}");
  PRho.Plot(Form("%s/Rho.%s", outputDir, suffix));

  Plotting1D PMultiplicity;
  PMultiplicity.New(hMultiplicity);
  PMultiplicity.SetAxisLabel("#bf{Multiplicity}", "#bf{Counts}");
  PMultiplicity.Plot(Form("%s/Multiplicity.%s", outputDir, suffix));
  
  TH2F* hRhoVsCentrality = (TH2F*)hCentralityRhoMultiplicity->Projection(1, 0);
  Plotting2D PRhoVsCentrality;
  PRhoVsCentrality.New(hRhoVsCentrality);
  PRhoVsCentrality.SetAxisLabel("#bf{Centrality (%)}", "#bf{#rho (GeV/#it{c})}");
  PRhoVsCentrality.Plot(Form("%s/RhoVsCentrality.%s", outputDir, suffix),false,false,true);

  TH2F* hRhoVsMultiplicity = (TH2F*)hCentralityRhoMultiplicity->Projection(1, 2);
  Plotting2D PRhoVsMultiplicity;
  PRhoVsMultiplicity.New(hRhoVsMultiplicity);
  PRhoVsMultiplicity.SetAxisLabel("#bf{Multiplicity}", "#bf{#rho} (GeV/#it{c})");
  PRhoVsMultiplicity.Plot(Form("%s/RhoVsMultiplicity.%s", outputDir, suffix), true);

  // Build an event histogram that has one bin called 0-10% and another called 10-30%,30-50%,50-90% and 0-90%
  TH1F* hCentralityBinned = new TH1F("hCentralityBinned", "Centrality", 5, 0, 5);
  hCentralityBinned->GetXaxis()->SetBinLabel(1, "0-10%");
  hCentralityBinned->GetXaxis()->SetBinLabel(2, "10-30%");
  hCentralityBinned->GetXaxis()->SetBinLabel(3, "30-50%");
  hCentralityBinned->GetXaxis()->SetBinLabel(4, "50-90%");
  hCentralityBinned->GetXaxis()->SetBinLabel(5, "0-90%");


  for(Int_t bin = 0; bin < hCentrality->GetNbinsX(); bin++)
  {
     double cent = hCentrality->GetBinCenter(bin);
     double content = hCentrality->GetBinContent(bin);
      if(cent <= 10)
      {
        hCentralityBinned->Fill(0., content);
      }
      if(cent > 10 && cent <= 30)
      {
        hCentralityBinned->Fill(1., content);
      }
      if(cent > 30 && cent <= 50)
      {
        hCentralityBinned->Fill(2., content);
      }
      if(cent > 50 && cent <= 90)
      {
        hCentralityBinned->Fill(3., content);
      } 
      if(cent >= 0 && cent <= 90)
      {
        hCentralityBinned->Fill(4., content);
      }
  }
  

  TCanvas* cCentralityBinned = new TCanvas("cCentralityBinned", "cCentralityBinned", 800, 800);
  hCentralityBinned->SetStats(0);
  hCentralityBinned->SetLineColor(kBlack);
  hCentralityBinned->SetLineWidth(3);
  hCentralityBinned->Draw("hist");
  // print number of events in each bin
  for (int i = 1; i <= hCentralityBinned->GetNbinsX(); i++)
  {
    float content = hCentralityBinned->GetBinContent(i);
    float error = hCentralityBinned->GetBinError(i);
    TLatex l;
    l.SetTextSize(0.03);
    // draw labels with the number of events rotated by 90 degrees
    l.DrawLatex(i-0.5,  1.1 * content, Form("%.0f", content, error));
  }
  cCentralityBinned->SaveAs(Form("%s/CentralityBinned.%s", outputDir, suffix));

  EXIT
}

void plotJets(TDirectory *dJets, GlobalOptions optns)
{
  ENTER

  const char *outputDir = Form("%s/Jets", optns.analysisDirPath.Data());
  createDirectory(outputDir);

  if (optns.isMC)
  {
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
    float yMax = 1;
    for (int iProj = 0; iProj < NProjections; iProj++)
    {
      int binProjFrom = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][0] + 0.001);
      int binProjTo = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][1] - 0.001);
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
      // gRatioDLJetPtToPLJetPtProjection[iProj]->GetXaxis()->SetRangeUser(0.5, 1.5);

      for (int iP = 0; iP < gRatioDLJetPtToPLJetPtProjection[iProj]->GetN(); iP++)
      {
        if (gRatioDLJetPtToPLJetPtProjection[iProj]->GetPointX(iP) > 0.6)
        {
          float thisYMax = gRatioDLJetPtToPLJetPtProjection[iProj]->GetPointY(iP);
          yMax = ((thisYMax > yMax) ? thisYMax : yMax);
          break;
        }
      }
      PResponseMatrixProjections.New(gRatioDLJetPtToPLJetPtProjection[iProj], Form("#it{p}_{T}^{PL} = %.0f #pm 1 GeV/#it{c}", PLJetpT));
    }
    PResponseMatrixProjections.SetAxisRange(0, 2.5, 0, 1.5 * yMax);
    PResponseMatrixProjections.SetLegend(0.5, 0.9, 0.5, 0.9);
    PResponseMatrixProjections.SetAxisLabel("#bf{#it{p}_{T}^{DL}/#it{p}_{T}^{PL}} #hat{=} #bf{#it{p}_{T}^{rec}/#it{p}_{T}^{gen}}");
    PResponseMatrixProjections.Plot(Form("%s/ResponseMatrixProjections.%s", outputDir, suffix));
  }

  EXIT
}

void plotJetQA(TDirectory *dJetQA,TString dirname,GlobalOptions optns)
{
  ENTER

  DLJetCuts jetCuts(optns);

  const char *outputDir = Form("%s/%s", optns.analysisDirPath.Data(), dirname.Data());
  createDirectory(outputDir);

  int const nPtBins = 7;
  Double_t jetPtBins[nPtBins+1] = {10,15,20,30,40,50,150,200};
  
  // Plot eta phi map of reconstructed jets (DL)
  THnSparseF *hJetPtEtaPhi = (THnSparseF *)dJetQA->Get("hJetPtEtaPhi");
  // axis 0: pT, axis 1: eta, axis 2: phi
  TH2F* hDLJetEtaPhi = (TH2F*)hJetPtEtaPhi->Projection(2, 1);
  hDLJetEtaPhi->SetName("hDLJetEtaPhi");
  Plotting2D PDLEtaPhiMap;
  PDLEtaPhiMap.SetMargins(0.12, 0.12, 0.05, 0.125);
  PDLEtaPhiMap.New(hDLJetEtaPhi);
  PDLEtaPhiMap.AddEMCalOutline();
  PDLEtaPhiMap.AddPHOSOutline();
  PDLEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PDLEtaPhiMap.Plot(Form("%s/DLEtaPhiMap.%s", outputDir, suffix));

  TH2F* hDLPtvsEta = (TH2F*)hJetPtEtaPhi->Projection(0, 1);
  hDLPtvsEta->SetName("hDLPtvsEta");
  Plotting2D PDLPtvsEta;
  PDLPtvsEta.SetMargins(0.12, 0.12, 0.05, 0.125);
  PDLPtvsEta.New(hDLPtvsEta);
  PDLPtvsEta.SetAxisLabel("#bf{#eta}", "#bf{#it{p}_{T} (GeV/#it{c})}");
  PDLPtvsEta.Plot(Form("%s/DLPtvsEta.%s", outputDir, suffix));
  
  // create eta projections
  // PlottingGrid PDLEtaProjections;
  // PDLEtaProjections.SetMargins(0.12, 0.12, 0.05, 0.125);
  // PDLEtaProjections.SetAxisLabel("#bf{#eta}", "#bf{Counts}");
  // TH1F* hDLEta[nPtBins];
  // for (int i = 0; i < nPtBins; i++)
  // {
  //   // set range of hPtEtaPhi to the current pT bin
  //   int minbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i]);
  //   int maxbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i + 1]);
  //   hJetPtEtaPhi->GetAxis(0)->SetRange(minbin, maxbin);
  //   hDLEta[i] = (TH1F*)hJetPtEtaPhi->Projection(1);
  //   hDLEta[i]->SetName(Form("hDLEta_%d", i));
  //   PDLEtaProjections.New(hDLEta[i],"",-1,1,kBlack,"hist");
  //   PDLEtaProjections.NextPad(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", jetPtBins[i], jetPtBins[i + 1]));
  // }
  // PDLEtaProjections.Plot(Form("%s/DLEtaProjections.%s", outputDir, suffix));

  TCanvas* cDLEtaProjections = new TCanvas("cDLEtaProjections", "cDLEtaProjections", 800, 800);
  cDLEtaProjections->Divide(2, 4);
  for (int i = 0; i < nPtBins; i++)
  {
    // set range of hPtEtaPhi to the current pT bin
    int minbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i]);
    int maxbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i + 1]);
    hJetPtEtaPhi->GetAxis(0)->SetRange(minbin, maxbin);
    TH1F* hDLEta = (TH1F*)hJetPtEtaPhi->Projection(1);
    hDLEta->SetName(Form("hDLEta_%d", i));
    cDLEtaProjections->cd(i+1);
    hDLEta->Draw();
  }
  cDLEtaProjections->SaveAs(Form("%s/DLEtaProjections.%s", outputDir, suffix));
  

  hJetPtEtaPhi->GetAxis(0)->SetRange(0, hJetPtEtaPhi->GetAxis(0)->GetNbins()-1);
  TH2F* hDLPtvsPhi = (TH2F*)hJetPtEtaPhi->Projection(0, 2);
  hDLPtvsPhi->SetName("hDLPtvsPhi");
  Plotting2D PDLPtvsPhi;
  PDLPtvsPhi.New(hDLPtvsPhi);
  PDLPtvsPhi.SetAxisLabel("#bf{#phi}", "#bf{#it{p}_{T} (GeV/#it{c})}");
  PDLPtvsPhi.Plot(Form("%s/DLPtvsPhi.%s", outputDir, suffix),false,false,true);

  // create phi projections
  // {
  // PlottingGrid PDLPhiProjections;
  // PDLPhiProjections.SetMargins(0.12, 0.12, 0.05, 0.125);
  // PDLPhiProjections.SetAxisLabel("#bf{#phi}", "#bf{Counts}");
  // TH1F* hDLPhi[nPtBins];
  // for (int i = 0; i < nPtBins; i++)
  // {
  //   // set range of hPtEtaPhi to the current pT bin
  //   int minbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i]);
  //   int maxbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i + 1]);
  //   hJetPtEtaPhi->GetAxis(0)->SetRange(minbin, maxbin);
  //   hDLPhi[i] = (TH1F*)hJetPtEtaPhi->Projection(2);
  //   hDLPhi[i]->SetName(Form("hDLPhi_%d", i));
  //   PDLPhiProjections.New(hDLPhi[i],"",-1,1,kBlack,"hist");
  //   PDLPhiProjections.NextPad(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", jetPtBins[i], jetPtBins[i + 1]));
  // }
  // PDLPhiProjections.Plot(Form("%s/DLPhiProjections.%s", outputDir, suffix));
  // }

  TCanvas* cDLPhiProjections = new TCanvas("cDLPhiProjections", "cDLPhiProjections", 800, 800);
  cDLPhiProjections->Divide(2, 4);
  for (int i = 0; i < nPtBins; i++)
  {
    // set range of hPtEtaPhi to the current pT bin
    int minbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i]);
    int maxbin = hJetPtEtaPhi->GetAxis(0)->FindBin(jetPtBins[i + 1]);
    hJetPtEtaPhi->GetAxis(0)->SetRange(minbin, maxbin);
    TH1F* hDLPhi = (TH1F*)hJetPtEtaPhi->Projection(2);
    hDLPhi->SetName(Form("hDLPhi_%d", i));
    cDLPhiProjections->cd(i+1);
    hDLPhi->Draw();
  }
  cDLPhiProjections->SaveAs(Form("%s/DLPhiProjectionsTest.%s", outputDir, suffix));


  
  // reset all ranges
  hJetPtEtaPhi->GetAxis(0)->SetRange(0,hJetPtEtaPhi->GetAxis(0)->GetNbins());
  hJetPtEtaPhi->GetAxis(1)->SetRange(0,hJetPtEtaPhi->GetAxis(1)->GetNbins());
  hJetPtEtaPhi->GetAxis(2)->SetRange(0,hJetPtEtaPhi->GetAxis(2)->GetNbins());
  
  TH1F* hDLPt = (TH1F*)hJetPtEtaPhi->Projection(0);
  Plotting1D PDLPt;
  PDLPt.SetMargins(0.12, 0.12, 0.05, 0.125);
  AvoidLogYconflict(hDLPt);
  PDLPt.New(hDLPt);
  PDLPt.SetAxisLabel("#bf{#it{p}_{T} (GeV/#it{c})}", "#bf{Counts}");
  PDLPt.Plot(Form("%s/DLPt.%s", outputDir, suffix),false,true);
  // Plot N constituents
  TH1F* hDLConst = (TH1F*)dJetQA->Get("hJetNconstits");
  Plotting1D PDLConst;
  PDLConst.New(hDLConst, "", -1, 1, kBlack, "hist");
  PDLConst.SetAxisLabel("#bf{N constituents}", "#bf{Counts}");
  PDLConst.Plot(Form("%s/DLConst.%s", outputDir, suffix));

  // Jet Area
  TH1F* hDLArea = (TH1F*)dJetQA->Get("hJetArea");
  Plotting1D PDLArea;
  PDLArea.New(hDLArea, "", -1, 1, kBlack, "hist");
  PDLArea.SetAxisLabel("#bf{Area}", "#bf{Counts}");
  PDLArea.Plot(Form("%s/DLArea.%s", outputDir, suffix));

  // Figures related to leading hadron pt

  // get THnSparseF hJetPtEtaPhi from dJetQA
  THnSparseF *hJetPtEtaZ = (THnSparseF *)dJetQA->Get("hJetPtEtaZ");

  // do 2D projection of Z vs pT
  TH2F* hDLZvsPt = (TH2F*)hJetPtEtaZ->Projection(2, 0);
  hDLZvsPt->SetName("hDLZvsPt");
  Plotting2D PDLZvsPt;
  PDLZvsPt.New(hDLZvsPt);
  PDLZvsPt.SetAxisLabel("#bf{#it{p}_{T} (GeV/#it{c})}", "#bf{z = #it{p}_{T}^{leading hadron} / #it{p}_{T}^{jet}}");
  PDLZvsPt.Plot(Form("%s/DLZvsPt.%s", outputDir, suffix),false,false,true);

  EXIT
}

void plotGammas(TDirectory *dIsoGammas, TString GammaType, GlobalOptions optns)
{
  ENTER

  TString outputDir = Form("%s/%s", optns.analysisDirPath.Data(), GammaType.Data());

  createDirectory(outputDir);

  // get THnSparseF for M02, Iso and Pt QA
  // 0: Iso, 1: M02, 2: Pt
  THnSparseF *hIsoGammaIsovsM02vsPt = (THnSparseF *)dIsoGammas->Get("hIsoGammaIsovsM02vsPt");
   

  // do two dimensional plots
  TH2F* hIsoVsM02 = (TH2F*)hIsoGammaIsovsM02vsPt->Projection(0, 1);
  hIsoVsM02->Sumw2();
  hIsoVsM02->Rebin2D(1, 2);
  Plotting2D PIsoVsM02;
  PIsoVsM02.SetAxisLabel("#bf{#sigma^{2}_{long}}", "#bf{Isolation (GeV/#it{c})}");
  PIsoVsM02.New(hIsoVsM02);
  PIsoVsM02.Plot(Form("%s/IsoVsM02.%s", outputDir.Data(), suffix), false, false, true);

  TH2F* hM02VsPt = (TH2F*)hIsoGammaIsovsM02vsPt->Projection(2, 1);
  hM02VsPt->Sumw2();
  hM02VsPt->Rebin2D(2, 1);
  Plotting2D PM02VsPt;
  PM02VsPt.SetAxisLabel("#bf{#sigma^{2}_{long}}", "#bf{#it{p}_{T} (GeV/#it{c})}");
  PM02VsPt.New(hM02VsPt);
  PM02VsPt.Plot(Form("%s/M02VsPt.%s", outputDir.Data(), suffix), false, false, true);

  PlottingGrid PM02Projections;
  PM02Projections.SetMargins(0.12, 0.12, 0.05, 0.125);
  PM02Projections.SetAxisLabel("#bf{#sigma^{2}_{long}}", "#bf{Norm. counts}");
  TH1F* hM02[nPtBins];

  for (int i = 0; i < nPtBins; i++)
  {
    // set range of hIsoGammaIsovsM02vsPt to the current pT bin
    int minbin = hIsoGammaIsovsM02vsPt->GetAxis(2)->FindBin(ptBins[i]);
    int maxbin = hIsoGammaIsovsM02vsPt->GetAxis(2)->FindBin(ptBins[i + 1]);
    hIsoGammaIsovsM02vsPt->GetAxis(2)->SetRange(minbin, maxbin);
    hM02[i] = (TH1F*)hIsoGammaIsovsM02vsPt->Projection(1);
    hM02[i]->SetName(Form("hM02_%d", i));
    hM02[i]->Scale(1. / hM02[i]->Integral());
    PM02Projections.New(hM02[i],"",-1,1,kBlack,"hist e");
    PM02Projections.NextPad(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", ptBins[i], ptBins[i + 1]));
  }
  PM02Projections.Plot(Form("%s/M02Projections.%s", outputDir.Data(), suffix));

}

void plotGammaQA(TDirectory *dIsoGammaQA, TString GammaType, GlobalOptions optns)
{
  ENTER

  TString outputDir = Form("%s/%s", optns.analysisDirPath.Data(), GammaType.Data());

  createDirectory(outputDir.Data());

  // Plot response matrix
  Plotting1D PMinMassDiffToPi0;
  TH1F *hMinMassDiffToPi0Signal;
  TH1F *hMinMassDiffToPi0Background;
  TH1F *hMinMassDiffToPi0;

  Plotting1D PIsoGammaE;
  TH1F *hIsoGammaESignal;
  TH1F *hIsoGammaEBackground;
  TH1F *hIsoGammaE;

  Plotting1D PIsoGammaM02;
  TH1F *hIsoGammaM02Signal;
  TH1F *hIsoGammaM02Background;
  TH1F *hIsoGammaM02;

  Plotting1D PIsoGammaM20;
  TH1F *hIsoGammaM20Signal;
  TH1F *hIsoGammaM20Background;
  TH1F *hIsoGammaM20;

  Plotting1D PIsoGammaPx;
  TH1F *hIsoGammaPxSignal;
  TH1F *hIsoGammaPxBackground;
  TH1F *hIsoGammaPx;

  Plotting1D PIsoGammaPy;
  TH1F *hIsoGammaPySignal;
  TH1F *hIsoGammaPyBackground;
  TH1F *hIsoGammaPy;

  Plotting1D PIsoGammaPz;
  TH1F *hIsoGammaPzSignal;
  TH1F *hIsoGammaPzBackground;
  TH1F *hIsoGammaPz;

  Plotting1D PIsoGammaIsoCharged;
  TH1F *hIsoGammaIsoChargedSignal;
  TH1F *hIsoGammaIsoChargedBackground;
  TH1F *hIsoGammaIsoCharged;

  Plotting1D PIsoGammaIsoChargedCorrected;
  TH1F *hIsoGammaIsoChargedCorrectedSignal;
  TH1F *hIsoGammaIsoChargedCorrectedBackground;
  TH1F *hIsoGammaIsoChargedCorrected;

  if (optns.isMC)
  {
    // Diff to pi0 mass
    hMinMassDiffToPi0Signal = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0Signal")->Clone("hMinMassDiffToPi0Signal");
    hMinMassDiffToPi0Background = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0Background")->Clone("hMinMassDiffToPi0Background");
    float integralSignal = hMinMassDiffToPi0Signal->Integral(1, hMinMassDiffToPi0Signal->GetNbinsX());
    float integralBG = hMinMassDiffToPi0Background->Integral(1, hMinMassDiffToPi0Background->GetNbinsX());
    hMinMassDiffToPi0Signal->Scale(1. / integralSignal);
    hMinMassDiffToPi0Background->Scale(1. / integralBG);
    AvoidLogYconflict(hMinMassDiffToPi0Signal);
    AvoidLogYconflict(hMinMassDiffToPi0Background);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Signal, "Signal");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Background, "Background");

    // E
    hIsoGammaESignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaESignal")->Clone("hIsoGammaESignal");
    hIsoGammaEBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaEBackground")->Clone("hIsoGammaEBackground");
    float integralESignal = hIsoGammaESignal->Integral(1, hIsoGammaESignal->GetNbinsX());
    float integralEBG = hIsoGammaEBackground->Integral(1, hIsoGammaEBackground->GetNbinsX());
    hIsoGammaESignal->Scale(1. / integralESignal);
    hIsoGammaEBackground->Scale(1. / integralEBG);
    AvoidLogYconflict(hIsoGammaESignal);
    AvoidLogYconflict(hIsoGammaEBackground);
    PIsoGammaE.New(hIsoGammaESignal, "Signal");
    PIsoGammaE.New(hIsoGammaEBackground, "Background");

    // M02
    hIsoGammaM02Signal = (TH1F *)dIsoGammaQA->Get("hIsoGammaM02Signal")->Clone("hIsoGammaM02Signal");
    hIsoGammaM02Background = (TH1F *)dIsoGammaQA->Get("hIsoGammaM02Background")->Clone("hIsoGammaM02Background");
    float integralM02Signal = hIsoGammaM02Signal->Integral(1, hIsoGammaM02Signal->GetNbinsX());
    float integralM02BG = hIsoGammaM02Background->Integral(1, hIsoGammaM02Background->GetNbinsX());
    hIsoGammaM02Signal->Scale(1. / integralM02Signal);
    hIsoGammaM02Background->Scale(1. / integralM02BG);
    AvoidLogYconflict(hIsoGammaM02Signal);
    AvoidLogYconflict(hIsoGammaM02Background);
    PIsoGammaM02.New(hIsoGammaM02Signal, "Signal");
    PIsoGammaM02.New(hIsoGammaM02Background, "Background");

    // M20
    hIsoGammaM20Signal = (TH1F *)dIsoGammaQA->Get("hIsoGammaM20Signal")->Clone("hIsoGammaM20Signal");
    hIsoGammaM20Background = (TH1F *)dIsoGammaQA->Get("hIsoGammaM20Background")->Clone("hIsoGammaM20Background");
    float integralM20Signal = hIsoGammaM20Signal->Integral(1, hIsoGammaM20Signal->GetNbinsX());
    float integralM20BG = hIsoGammaM20Background->Integral(1, hIsoGammaM20Background->GetNbinsX());
    hIsoGammaM20Signal->Scale(1. / integralM20Signal);
    hIsoGammaM20Background->Scale(1. / integralM20BG);
    AvoidLogYconflict(hIsoGammaM20Signal);
    AvoidLogYconflict(hIsoGammaM20Background);
    PIsoGammaM20.New(hIsoGammaM20Signal, "Signal");
    PIsoGammaM20.New(hIsoGammaM20Background, "Background");

    // Px
    hIsoGammaPxSignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaPxSignal")->Clone("hIsoGammaPxSignal");
    hIsoGammaPxBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaPxBackground")->Clone("hIsoGammaPxBackground");
    float integralPxSignal = hIsoGammaPxSignal->Integral(1, hIsoGammaPxSignal->GetNbinsX());
    float integralPxBG = hIsoGammaPxBackground->Integral(1, hIsoGammaPxBackground->GetNbinsX());
    hIsoGammaPxSignal->Scale(1. / integralPxSignal);
    hIsoGammaPxBackground->Scale(1. / integralPxBG);
    AvoidLogYconflict(hIsoGammaPxSignal);
    AvoidLogYconflict(hIsoGammaPxBackground);
    PIsoGammaPx.New(hIsoGammaPxSignal, "Signal");
    PIsoGammaPx.New(hIsoGammaPxBackground, "Background");

    // Py
    hIsoGammaPySignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaPySignal")->Clone("hIsoGammaPySignal");
    hIsoGammaPyBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaPyBackground")->Clone("hIsoGammaPyBackground");
    float integralPySignal = hIsoGammaPySignal->Integral(1, hIsoGammaPySignal->GetNbinsX());
    float integralPyBG = hIsoGammaPyBackground->Integral(1, hIsoGammaPyBackground->GetNbinsX());
    hIsoGammaPySignal->Scale(1. / integralPySignal);
    hIsoGammaPyBackground->Scale(1. / integralPyBG);
    AvoidLogYconflict(hIsoGammaPySignal);
    AvoidLogYconflict(hIsoGammaPyBackground);
    PIsoGammaPy.New(hIsoGammaPySignal, "Signal");
    PIsoGammaPy.New(hIsoGammaPyBackground, "Background");

    // Pz
    hIsoGammaPzSignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaPzSignal")->Clone("hIsoGammaPzSignal");
    hIsoGammaPzBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaPzBackground")->Clone("hIsoGammaPzBackground");
    float integralPzSignal = hIsoGammaPzSignal->Integral(1, hIsoGammaPzSignal->GetNbinsX());
    float integralPzBG = hIsoGammaPzBackground->Integral(1, hIsoGammaPzBackground->GetNbinsX());
    hIsoGammaPzSignal->Scale(1. / integralPzSignal);
    hIsoGammaPzBackground->Scale(1. / integralPzBG);
    AvoidLogYconflict(hIsoGammaPzSignal);
    AvoidLogYconflict(hIsoGammaPzBackground);
    PIsoGammaPz.New(hIsoGammaPzSignal, "Signal");
    PIsoGammaPz.New(hIsoGammaPzBackground, "Background");

    // IsoCharged
    hIsoGammaIsoChargedSignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoChargedSignal")->Clone("hIsoGammaIsoChargedSignal");
    hIsoGammaIsoChargedBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoChargedBackground")->Clone("hIsoGammaIsoChargedBackground");
    float integralIsolationSignal = hIsoGammaIsoChargedSignal->Integral(1, hIsoGammaIsoChargedSignal->GetNbinsX());
    float integralIsolationBG = hIsoGammaIsoChargedBackground->Integral(1, hIsoGammaIsoChargedBackground->GetNbinsX());
    hIsoGammaIsoChargedSignal->Scale(1. / integralIsolationSignal);
    hIsoGammaIsoChargedBackground->Scale(1. / integralIsolationBG);
    AvoidLogYconflict(hIsoGammaIsoChargedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedBackground);
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedSignal, "Signal");
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedBackground, "Background");

    // IsoChargedCorrected
    hIsoGammaIsoChargedCorrectedSignal = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoChargedCorrectedSignal")->Clone("hIsoGammaIsoChargedCorrectedSignal");
    hIsoGammaIsoChargedCorrectedBackground = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoChargedCorrectedBackground")->Clone("hIsoGammaIsoChargedCorrectedBackground");
    float integralIsolationCorrectedSignal = hIsoGammaIsoChargedCorrectedSignal->Integral(1, hIsoGammaIsoChargedCorrectedSignal->GetNbinsX());
    float integralIsolationCorrectedBG = hIsoGammaIsoChargedCorrectedBackground->Integral(1, hIsoGammaIsoChargedCorrectedBackground->GetNbinsX());
    hIsoGammaIsoChargedCorrectedSignal->Scale(1. / integralIsolationSignal);
    hIsoGammaIsoChargedCorrectedBackground->Scale(1. / integralIsolationBG);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedBackground);
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedSignal, "Signal");
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedBackground, "Background");

    // M02vsPt (Slicing 2D hist pr. specified pT-bins)
    TH2F *h2M02vspT;
    TH2F *h2M02vspTSignal;
    TH2F *h2M02vspTBackground;
    TH2F *h2M02vspTAllGamma;
    TH2F *h2M02vspTPionDecayGamma;
    TH2F *h2M02vspTEtaDecayGamma;
    TH2F *h2M02vspTMergedPionGamma;

    h2M02vspT = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pT")->Clone("hIsoGammaM02pT");
    h2M02vspTSignal = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTSignal")->Clone("hIsoGammaM02pTSignal");
    h2M02vspTBackground = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTBackground")->Clone("hIsoGammaM02pTBackground");
    h2M02vspTAllGamma = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTAllGamma")->Clone("hIsoGammaM02pTAllGamma");
    h2M02vspTPionDecayGamma = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTPionDecayGamma")->Clone("hIsoGammaM02pTPionDecayGamma");
    h2M02vspTEtaDecayGamma = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTEtaDecayGamma")->Clone("hIsoGammaM02pTEtaDecayGamma");
    h2M02vspTMergedPionGamma = (TH2F *)dIsoGammaQA->Get("hIsoGammaM02pTMergedPionGamma")->Clone("hIsoGammaM02pTMergedPionGamma");

    int NSlicesPtBins = SlicesPt.size();

    // Initiate arrays for storing Plotting1D's for loop
    Plotting1D PM02vspTSignalBackground[NSlicesPtBins];
    Plotting1D PM02vspTSignalContributions[NSlicesPtBins];

    TH1F *hM02vspTSlice[NSlicesPtBins];
    TH1F *hM02vspTSliceSignal[NSlicesPtBins];
    TH1F *hM02vspTSliceBackground[NSlicesPtBins];
    TH1F *hM02vspTSliceAllGamma[NSlicesPtBins];
    TH1F *hM02vspTSlicePionDecayGamma[NSlicesPtBins];
    TH1F *hM02vspTSliceEtaDecayGamma[NSlicesPtBins];
    TH1F *hM02vspTSliceMergedPionGamma[NSlicesPtBins];

    PlottingGrid PM02Grid;
    PM02Grid.SetAxisLabel("#bf{#it{M}_{02}}", "#bf{#it{N}}");

    for (int i = 0; i < NSlicesPtBins - 1; i++)
    {
      // Project TH2F's
      string titleM02vspTsplice = "M02:" + to_string(SlicesPt.at(i)).substr(0, 4) + "GeV/c<p_{T}<" + to_string(SlicesPt.at(i + 1)).substr(0, 4) + "GeV/c";
      int pTminbin = h2M02vspT->GetYaxis()->FindBin(SlicesPt.at(i));
      int pTmaxbin = h2M02vspT->GetYaxis()->FindBin(SlicesPt.at(i + 1));

      hM02vspTSlice[i] = (TH1F *)h2M02vspT->ProjectionX(titleM02vspTsplice.c_str(), pTminbin, pTmaxbin);
      hM02vspTSlice[i]->SetTitle(titleM02vspTsplice.c_str());
      hM02vspTSliceSignal[i] = (TH1F *)h2M02vspTSignal->ProjectionX(Form("a%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceBackground[i] = (TH1F *)h2M02vspTBackground->ProjectionX(Form("b%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceAllGamma[i] = (TH1F *)h2M02vspTAllGamma->ProjectionX(Form("c%i", i), pTminbin, pTmaxbin);
      hM02vspTSlicePionDecayGamma[i] = (TH1F *)h2M02vspTPionDecayGamma->ProjectionX(Form("d%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceEtaDecayGamma[i] = (TH1F *)h2M02vspTEtaDecayGamma->ProjectionX(Form("e%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceMergedPionGamma[i] = (TH1F *)h2M02vspTMergedPionGamma->ProjectionX(Form("f%i", i), pTminbin, pTmaxbin);
      // Unique hist name.
      // Try array.
      AvoidLogYconflict(hM02vspTSlice[i]);
      AvoidLogYconflict(hM02vspTSliceSignal[i]);
      AvoidLogYconflict(hM02vspTSliceBackground[i]);
      AvoidLogYconflict(hM02vspTSliceAllGamma[i]);
      AvoidLogYconflict(hM02vspTSlicePionDecayGamma[i]);
      AvoidLogYconflict(hM02vspTSliceMergedPionGamma[i]);

      // Add to figure and plot:
      // All, Signal and Background
      PM02vspTSignalBackground[i].New(hM02vspTSlice[i], "All");
      PM02vspTSignalBackground[i].New(hM02vspTSliceSignal[i], "Signal");
      PM02vspTSignalBackground[i].New(hM02vspTSliceBackground[i], "Background");
      PM02vspTSignalBackground[i].SetAxisLabel(("M02"));
      PM02vspTSignalBackground[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalBackgroundSig.%s", outputDir.Data(), SlicesPt.at(i), SlicesPt.at(i + 1), suffix), kFALSE, kTRUE);
      //
      ////All, Signal and Contributions
      PM02vspTSignalContributions[i].New(hM02vspTSlice[i], "All");
      PM02vspTSignalContributions[i].New(hM02vspTSliceSignal[i], "Signal");
      PM02vspTSignalContributions[i].New(hM02vspTSlicePionDecayGamma[i], "#pi^{0}-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceEtaDecayGamma[i], "#eta-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceMergedPionGamma[i], "#pi^{0}");
      PM02vspTSignalContributions[i].SetAxisLabel(("M02"));
      PM02vspTSignalContributions[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalContributions.%s", outputDir.Data(), SlicesPt.at(i), SlicesPt.at(i + 1), suffix), kFALSE, kTRUE);

      PM02Grid.New(hM02vspTSlice[i], "All");
      PM02Grid.New(hM02vspTSliceSignal[i], "Signal");
      PM02Grid.New(hM02vspTSliceBackground[i], "Background");
      PM02Grid.NextPad(Form("%.1f < #it{p}_{T} (GeV/c) < %.1f", SlicesPt.at(i), SlicesPt.at(i + 1)));
    }
    PM02Grid.Plot(Form("%s/M02_SignalBackgroundSig_Grid.%s", outputDir.Data(), suffix), kFALSE, kFALSE);

    //
  }
  else
  {
    hMinMassDiffToPi0 = (TH1F *)dIsoGammaQA->Get("hMinMassDiffToPi0");
    hMinMassDiffToPi0->SetTitle("");
    float integralDiffToPi0 = hMinMassDiffToPi0->Integral(1, hMinMassDiffToPi0->GetNbinsX());
    hMinMassDiffToPi0->Scale(1. / integralDiffToPi0);
    AvoidLogYconflict(hMinMassDiffToPi0);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0);

    // E
    hIsoGammaE = (TH1F *)dIsoGammaQA->Get("hIsoGammaE")->Clone("hIsoGammaE");
    float integralE = hIsoGammaE->Integral(1, hIsoGammaE->GetNbinsX());
    hIsoGammaE->Scale(1. / integralE);
    AvoidLogYconflict(hIsoGammaE);
    PIsoGammaE.New(hIsoGammaE, "Data");

    // M02
    hIsoGammaM02 = (TH1F *)dIsoGammaQA->Get("hIsoGammaM02")->Clone("hIsoGammaM02");
    float integralM02 = hIsoGammaM02->Integral(1, hIsoGammaM02->GetNbinsX());
    hIsoGammaM02->Scale(1. / integralM02);
    AvoidLogYconflict(hIsoGammaM02);
    if (GammaType == "Clusters")
    {
      PIsoGammaM02.New(hIsoGammaM02, "Data", 1, 2, kBlack, "h");

      TH1F *hIsogammaSelectedM02 = (TH1F *)hIsoGammaM02->Clone("hIsogammaSelectedM02");
      TH1F *hPi0SelectedM02 = (TH1F *)hIsoGammaM02->Clone("hPi0SelectedM02");
      for (int ibin = 0; ibin < hIsoGammaM02->GetNbinsX(); ibin++)
      {
        if (hIsogammaSelectedM02->GetBinCenter(ibin) < 0.1 || hIsogammaSelectedM02->GetBinCenter(ibin) > 0.3)
        {
          hIsogammaSelectedM02->SetBinContent(ibin, 0.);
          hIsogammaSelectedM02->SetBinError(ibin, 0.);
        }
        if (hPi0SelectedM02->GetBinCenter(ibin) < 0.4)
        {
          hPi0SelectedM02->SetBinContent(ibin, 0.);
          hPi0SelectedM02->SetBinError(ibin, 0.);
        }
      }

      PIsoGammaM02.New(hIsogammaSelectedM02, "#gamma candidates", 1, 4, GoetheColor[1][0], "hf");
      PIsoGammaM02.New(hPi0SelectedM02, "Merged #pi^{0} candidates", 1, 4, GoetheColor[3][0], "hf");
      PIsoGammaM02.New(hIsoGammaM02, "", 1, 4, kBlack, "h");
    }
    else
    {
      PIsoGammaM02.New(hIsoGammaM02, "Data");
    }

    // M20
    hIsoGammaM20 = (TH1F *)dIsoGammaQA->Get("hIsoGammaM20")->Clone("hIsoGammaM20");
    float integralM20 = hIsoGammaM20->Integral(1, hIsoGammaM20->GetNbinsX());
    hIsoGammaM20->Scale(1. / integralM20);
    AvoidLogYconflict(hIsoGammaM20);
    PIsoGammaM20.New(hIsoGammaM20, "Data");

    // Px
    hIsoGammaPx = (TH1F *)dIsoGammaQA->Get("hIsoGammaPx")->Clone("hIsoGammaPx");
    float integralPx = hIsoGammaPx->Integral(1, hIsoGammaPx->GetNbinsX());
    hIsoGammaPx->Scale(1. / integralPx);
    AvoidLogYconflict(hIsoGammaPx);
    PIsoGammaPx.New(hIsoGammaPx, "Data");

    // Py
    hIsoGammaPy = (TH1F *)dIsoGammaQA->Get("hIsoGammaPy")->Clone("hIsoGammaPy");
    float integralPy = hIsoGammaPy->Integral(1, hIsoGammaPy->GetNbinsX());
    hIsoGammaPy->Scale(1. / integralPy);
    AvoidLogYconflict(hIsoGammaPy);
    PIsoGammaPy.New(hIsoGammaPy, "Data");

    // Pz
    hIsoGammaPz = (TH1F *)dIsoGammaQA->Get("hIsoGammaPz")->Clone("hIsoGammaPz");
    float integralPz = hIsoGammaPz->Integral(1, hIsoGammaPz->GetNbinsX());
    hIsoGammaPz->Scale(1. / integralPz);
    AvoidLogYconflict(hIsoGammaPz);
    PIsoGammaPz.New(hIsoGammaPz, "Data");

    // IsoCharged
    hIsoGammaIsoCharged = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoCharged")->Clone("hIsoGammaIsoCharged");
    float integralIsoCharged = hIsoGammaIsoCharged->Integral(1, hIsoGammaIsoCharged->GetNbinsX());
    hIsoGammaIsoCharged->Scale(1. / integralIsoCharged);
    AvoidLogYconflict(hIsoGammaIsoCharged);
    PIsoGammaIsoCharged.New(hIsoGammaIsoCharged, "Data");

    // IsoChargedCorrected
    hIsoGammaIsoChargedCorrected = (TH1F *)dIsoGammaQA->Get("hIsoGammaIsoChargedCorrected")->Clone("hIsoGammaIsoChargedCorrected");
    float integralIsoChargedCorrected = hIsoGammaIsoChargedCorrected->Integral(1, hIsoGammaIsoChargedCorrected->GetNbinsX());
    // hIsoGammaIsoChargedCorrected->Scale(1. / integralIsoChargedCorrected);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrected);
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrected, "");
  }
  PMinMassDiffToPi0.SetAxisLabel("#bf{#Delta#it{m}_{inv} (GeV/#it{c}^{2})}", "#bf{dN/d#Delta#it{m}_{inv} (GeV/#it{c}^{2})}");
  PMinMassDiffToPi0.SetMargins(0.12, 0.1, 0.08, 0.025);
  PMinMassDiffToPi0.Plot(Form("%s/MinMassDiffToPi0.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaE.SetAxisLabel("#bf{#it{E} (GeV)}", "#bf{dN/dE GeV}");
  PIsoGammaE.Plot(Form("%s/IsoGammaE.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaM02.SetAxisLabel("#bf{#it{#sigma}_{long}^{2}}", "#bf{Normalized counts}", 1., 1.3);
  PIsoGammaM02.SetAxisRange(0., 2.,1E-3,0.3);
  PIsoGammaM02.SetLegend(0.4, 0.7, 0.65, 0.85, true);
  PIsoGammaM02.NewLatex(0.9, 0.9, "ALICE work in progress;2018 data, #sqrt{s} = 13 TeV", StdTextSize, 0.04);
  PIsoGammaM02.SetMargins(0.11, 0.11, 0.025, 0.025, 1750, 1750);
  PIsoGammaM02.Plot(Form("%s/IsoGammaM02.%s", outputDir.Data(), suffix), kFALSE, kTRUE, true);

  PIsoGammaM20.SetAxisLabel("#bf{#it{M20}}", "#bf{log(dN/dM20)}");
  PIsoGammaM20.Plot(Form("%s/IsoGammaM20.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaPx.SetAxisLabel("#bf{#it{p_x}(GeV/#it{c})}", "#bf{log(dN/Px GeV/c)}");
  PIsoGammaPx.Plot(Form("%s/IsoGammaPx.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaPy.SetAxisLabel("#bf{#it{p_y}(GeV/#it{c})}", "#bf{log(dN/Py GeV/c)}");
  PIsoGammaPy.Plot(Form("%s/IsoGammaPy.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaIsoChargedCorrected.SetAxisLabel("#bf{#it{p}^{iso, corrected}_{T} (GeV/#it{c})}", "#bf{N_{iso #gamma}}", 0.9);
  PIsoGammaIsoChargedCorrected.Plot(Form("%s/IsoGammaIsoChargedCorrected.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaPz.SetAxisLabel("#bf{#it{p_z}(GeV/#it{c})}", "#bf{dN/p_{z} (GeV/c)}");
  PIsoGammaPz.Plot(Form("%s/IsoGammaPz.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  PIsoGammaIsoCharged.SetAxisLabel("#bf{#it{p^{iso}_{T}}(GeV/#it{c})}", "#bf{dN/p_{T}^{iso} (GeV/c)}");
  PIsoGammaIsoCharged.Plot(Form("%s/IsoGammaIsoCharged.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

  TH2F *hIsoGammaEtaPhi = (TH2F *)dIsoGammaQA->Get("hIsoGammaEtaPhi");
  Plotting2D PEtaPhiMap;
  PEtaPhiMap.SetMargins(0.12, 0.12, 0.05, 0.125);
  PEtaPhiMap.New(hIsoGammaEtaPhi);
  PEtaPhiMap.AddEMCalOutline();
  PEtaPhiMap.AddPHOSOutline();
  PEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PEtaPhiMap.Plot(Form("%s/EtaPhiMap.%s", outputDir.Data(), suffix));

  // get ThNSparse to get energy dependence
  THnSparseF *hIsoGammaEtaPhiPt = (THnSparseF *)dIsoGammaQA->Get("hIsoGammaEtaPhiPt");
  TH1F* hIsoGammaEtaProjections[nPtBins];
  TH1F* hIsoGammaPhiProjections[nPtBins];
  for (int iPtBin = 0; iPtBin < nPtBins; iPtBin++)
  {
    int minBin = hIsoGammaEtaPhiPt->GetAxis(2)->FindBin(ptBins[iPtBin]);
    int maxBin = hIsoGammaEtaPhiPt->GetAxis(2)->FindBin(ptBins[iPtBin + 1]);
    
    hIsoGammaEtaPhiPt->GetAxis(2)->SetRange(minBin, maxBin);
    hIsoGammaEtaProjections[iPtBin] = (TH1F*)hIsoGammaEtaPhiPt->Projection(0);
    hIsoGammaEtaProjections[iPtBin]->SetName(Form("hIsoGammaEtaProjection_%d", iPtBin));
    hIsoGammaPhiProjections[iPtBin] = (TH1F*)hIsoGammaEtaPhiPt->Projection(1);
    hIsoGammaPhiProjections[iPtBin]->SetName(Form("hIsoGammaPhiProjection_%d", iPtBin));
  }

  PlottingGrid PGammaEtaProjections;
  PGammaEtaProjections.SetAxisLabel("#bf{#eta}", "#bf{N_{clusters}}");
  PGammaEtaProjections.SetMargins(0.12, 0.12, 0.05, 0.125);
  DEBUG
  PlottingGrid PGammaPhiProjections;
  PGammaPhiProjections.SetAxisLabel("#bf{#phi}", "#bf{N_{clusters}}");
  PGammaPhiProjections.SetMargins(0.12, 0.12, 0.05, 0.125);

  for (int iPtBin = 0; iPtBin < nPtBins; iPtBin++)
  {
    PGammaEtaProjections.New(hIsoGammaEtaProjections[iPtBin], "", -1, 1, kBlack, "hist e");
    PGammaEtaProjections.NextPad(Form("%.1f < #it{p}_{T}^{cluster} < %.1f GeV/#it{c}", ptBins[iPtBin], ptBins[iPtBin + 1]));

    PGammaPhiProjections.New(hIsoGammaPhiProjections[iPtBin], "", -1, 1, kBlack, "hist e");
    PGammaPhiProjections.NextPad(Form("%.1f < #it{p}_{T}^{cluster} < %.1f GeV/#it{c}", ptBins[iPtBin], ptBins[iPtBin + 1]));
  }
  PGammaEtaProjections.Plot(Form("%s/EtaProjections.%s", outputDir.Data(), suffix), false, false);
  PGammaPhiProjections.Plot(Form("%s/PhiProjections.%s", outputDir.Data(), suffix), false, false);

  // plot cluster time
  TH2F* hIsoGammaTimePt = (TH2F*)dIsoGammaQA->Get("hIsoGammaTimePt");
  Plotting2D PTimePt;
  PTimePt.SetMargins(0.12, 0.12, 0.05, 0.125);
  PTimePt.New(hIsoGammaTimePt);
  PTimePt.SetAxisLabel("#bf{t_{cluster} (ns)}", "#bf{#it{p}_{T}^{cluster} (GeV/#it{c})}");
  PTimePt.Plot(Form("%s/TimePt.%s", outputDir.Data(), suffix),false,false,true);

  PlottingGrid PTimeProjections;
  PTimeProjections.SetAxisLabel("#bf{t_{cluster} (ns)}", "#bf{N_{clusters}}");
  PTimeProjections.SetMargins(0.12, 0.12, 0.05, 0.125);

  TH1F* hIsoGammaTime[nPtBins];
  for (int iPtBin = 0; iPtBin < nPtBins; iPtBin++)
  {
    int minBin = hIsoGammaTimePt->GetYaxis()->FindBin(ptBins[iPtBin]);
    int maxBin = hIsoGammaTimePt->GetYaxis()->FindBin(ptBins[iPtBin + 1]);
    hIsoGammaTime[iPtBin] = (TH1F*)hIsoGammaTimePt->ProjectionX(Form("hIsoGammaTime_%d", iPtBin), minBin, maxBin);
    hIsoGammaTime[iPtBin]->GetXaxis()->SetRangeUser(-105, 105);

    // get mean and sigma of time distribution
    double mean = hIsoGammaTime[iPtBin]->GetMean();
    double sigma = hIsoGammaTime[iPtBin]->GetRMS();
  
    PTimeProjections.New(hIsoGammaTime[iPtBin], "", -1, 1, kBlack, "hist e");
    // put label on plot with mean and sigma of time distribution
    PTimeProjections.NextPad(Form("%.1f < #it{p}_{T}^{cluster} < %.1f GeV/#it{c}", ptBins[iPtBin], ptBins[iPtBin + 1]));
  }
  PTimeProjections.Plot(Form("%s/TimeProjections.%s", outputDir.Data(), suffix), false, false);

  // NCells
  TH2F *hIsoGammaNCellsPt = (TH2F *)dIsoGammaQA->Get("hIsoGammaNCellsPt");
  Plotting2D PNCellsPt;
  PNCellsPt.SetMargins(0.12, 0.12, 0.05, 0.125);
  PNCellsPt.New(hIsoGammaNCellsPt);
  PNCellsPt.SetAxisLabel("#bf{N_{cells}}", "#bf{#it{p}_{T}^{cluster} (GeV/#it{c})}");
  PNCellsPt.Plot(Form("%s/NCellsPt.%s", outputDir.Data(), suffix));

  // projections
  PlottingGrid PNCellsProjections;
  PNCellsProjections.SetAxisLabel("#bf{N_{cells}}", "#bf{N_{clusters}}");
  PNCellsProjections.SetMargins(0.12, 0.12, 0.05, 0.125);

  TH1F *hIsoGammaNCells[nPtBins];
  for (int iPtBin = 0; iPtBin < nPtBins; iPtBin++)
  {
    int binMin = hIsoGammaNCellsPt->GetYaxis()->FindBin(ptBins[iPtBin]);
    int binMax = hIsoGammaNCellsPt->GetYaxis()->FindBin(ptBins[iPtBin + 1]);
    hIsoGammaNCells[iPtBin] = (TH1F *)hIsoGammaNCellsPt->ProjectionX(Form("hIsoGammaNCells_%d", iPtBin), binMin, binMax);
    hIsoGammaNCells[iPtBin]->GetXaxis()->SetRangeUser(0, 40);
    PNCellsProjections.New(hIsoGammaNCells[iPtBin], "", -1, 1, kBlack, "hist e");
    PNCellsProjections.NextPad(Form("%.1f < #it{p}_{T}^{cluster} < %.1f GeV/#it{c}", ptBins[iPtBin], ptBins[iPtBin + 1]));
  }
  PNCellsProjections.Plot(Form("%s/NCellsProjections.%s", outputDir.Data(), suffix), false, false);

  // plot individual and subsequen rejections of cuts
  // get hpTSpectrumLossFromIndividualCuts
  TH2F *hpTSpectrumLossFromIndividualCuts = (TH2F *)dIsoGammaQA->Get("hpTSpectrumLossFromIndividualCuts");
  Plotting2D PPtSpectrumLossFromIndividualCuts;
  PPtSpectrumLossFromIndividualCuts.New(hpTSpectrumLossFromIndividualCuts);
  PPtSpectrumLossFromIndividualCuts.Plot(Form("%s/pTSpectrumLossFromIndividualCuts.%s", outputDir.Data(), suffix),false,false,true);

  // hpTSpectraAfterSubsequentCuts
  TH2F *hpTSpectraAfterSubsequentCuts = (TH2F *)dIsoGammaQA->Get("hpTSpectraAfterSubsequentCuts");
  Plotting2D PPtSpectraAfterSubsequentCuts;
  PPtSpectraAfterSubsequentCuts.New(hpTSpectraAfterSubsequentCuts);
  PPtSpectraAfterSubsequentCuts.Plot(Form("%s/pTSpectraAfterSubsequentCuts.%s", outputDir.Data(), suffix),false,false,true);


  



  EXIT
}

void plotExclusiveSelectionGammaJetCorrelations(TDirectory *dGammaJetSignal, TDirectory *dGammaJetReference, TDirectory *dGammaSignal, TDirectory *dGammaSignalQA, TDirectory *dGammaReference, TDirectory *dGammaReferenceQA, TString GammaType, GlobalOptions optns)
{
  ENTER
  ExclusiveTriggerParticleSelection exclusiveTriggerParticleSelection(optns);

  TString outputDir = Form("%s/%s", optns.analysisDirPath.Data(), GammaType.Data());
  createDirectory(outputDir.Data());
  // Check if input histograms exist
  if (!dGammaJetSignal || !dGammaJetReference || !dGammaSignal || !dGammaReference) {
    printf("Error: One or more input directories are null\n");
    return;
  }

  // Get signal THnSparseF and check if valid
  THnSparseF* hGammaJetCorrelationsSignal = (THnSparseF*)dGammaJetSignal->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
  THnSparseF* hGammaJetCorrelationsReference = (THnSparseF*)dGammaJetReference->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
  
  if (!hGammaJetCorrelationsSignal || !hGammaJetCorrelationsReference) {
    printf("Error: Could not get correlation histograms\n");
    return;
  }

  // Get signal pt distributions and check if valid
  TH1F* hGammaPtSignal = (TH1F*)dGammaSignal->Get("hIsoGammaPt");
  hGammaPtSignal->SetName("hGammaPtSignal");
  TH1F* hGammaPtReference = (TH1F*)dGammaReference->Get("hIsoGammaPt"); 
  hGammaPtReference->SetName("hGammaPtReference");
  
  if (!hGammaPtSignal || !hGammaPtReference) {
    printf("Error: Could not get pt histograms\n");
    return;
  }

  double nSignal = hGammaPtSignal->Integral();
  double nReference = hGammaPtReference->Integral();

  if (nSignal > 0) hGammaJetCorrelationsSignal->Scale(1./nSignal);
  if (nReference > 0) hGammaJetCorrelationsReference->Scale(1./nReference);

  // Create projections of jet pt vs deltaPhi
  TH2F* hGammaJetDeltaPhiJetPtSignal = (TH2F*)hGammaJetCorrelationsSignal->Projection(1,0);
  hGammaJetDeltaPhiJetPtSignal->SetName("hGammaJetDeltaPhiJetPtSignal");
  TH2F* hGammaJetDeltaPhiJetPtReference = (TH2F*)hGammaJetCorrelationsReference->Projection(1,0);
  hGammaJetDeltaPhiJetPtReference->SetName("hGammaJetDeltaPhiJetPtReference");

  if (!hGammaJetDeltaPhiJetPtSignal || !hGammaJetDeltaPhiJetPtReference) {
    printf("Error: Projection failed\n");
    return;
  }

  // Do 2D plots
  Plotting2D PGammaJetDeltaPhiJetPtSignal;
  PGammaJetDeltaPhiJetPtSignal.New(hGammaJetDeltaPhiJetPtSignal);
  PGammaJetDeltaPhiJetPtSignal.SetAxisLabel("#bf{#Delta#phi}", "#bf{#it{p}_{T}^{jet} (GeV/#it{c})}");
  PGammaJetDeltaPhiJetPtSignal.NewLatex(0.9,0.9,Form("ALICE work-in-progress;TT{%.0f,%.0f} GeV/#it{c}", 
                                              exclusiveTriggerParticleSelection.getPtMinSignal(),
                                              exclusiveTriggerParticleSelection.getPtMaxSignal()));
  PGammaJetDeltaPhiJetPtSignal.Plot(Form("%s/GammaJetDeltaPhiJetPtSignal.%s", outputDir.Data(), suffix),false,false,true);

  Plotting2D PGammaJetDeltaPhiJetPtReference;
  PGammaJetDeltaPhiJetPtReference.New(hGammaJetDeltaPhiJetPtReference);
  PGammaJetDeltaPhiJetPtReference.SetAxisLabel("#bf{#Delta#phi}", "#bf{#it{p}_{T}^{jet} (GeV/#it{c})}");
  PGammaJetDeltaPhiJetPtReference.NewLatex(0.9,0.9,Form("TT{%.0f,%.0f} GeV/#it{c}", 
                                                 exclusiveTriggerParticleSelection.getPtMinReference(),
                                                 exclusiveTriggerParticleSelection.getPtMaxReference()));
  PGammaJetDeltaPhiJetPtReference.Plot(Form("%s/GammaJetDeltaPhiJetPtReference.%s", outputDir.Data(), suffix),false,false,true);

  // Do a projection of the distribution for back-to-back recoil jets
  hGammaJetCorrelationsSignal->GetAxis(0)->SetRangeUser(TMath::Pi()-0.6,TMath::Pi()+0.6);
  TH1F* hRecoilJetPtSignal = (TH1F*)hGammaJetCorrelationsSignal->Projection(1);
  hRecoilJetPtSignal->SetName("hRecoilJetPtSignal");
  hRecoilJetPtSignal->Rebin(2);
  hGammaJetCorrelationsReference->GetAxis(0)->SetRangeUser(TMath::Pi()-0.6,TMath::Pi()+0.6);
  TH1F* hRecoilJetPtReference = (TH1F*)hGammaJetCorrelationsReference->Projection(1);
  hRecoilJetPtReference->SetName("hRecoilJetPtReference");
  hRecoilJetPtReference->Rebin(2);

  if (!hRecoilJetPtSignal || !hRecoilJetPtReference) {
    printf("Error: Recoil jet projection failed\n");
    return;
  }

  // Divide histograms by bin width
  if (hRecoilJetPtSignal->GetBinWidth(1) > 0) hRecoilJetPtSignal->Scale(1./hRecoilJetPtSignal->GetBinWidth(1));
  if (hRecoilJetPtReference->GetBinWidth(1) > 0) hRecoilJetPtReference->Scale(1./hRecoilJetPtReference->GetBinWidth(1));

  // Plot signal and reference distribution in same plot
  Plotting1D PRecoilJetPt;
  AvoidLogYconflict(hRecoilJetPtSignal);
  AvoidLogYconflict(hRecoilJetPtReference);

  PRecoilJetPt.SetAxisLabel("#bf{#it{p}_{T, ch jet}^{reco} (GeV/#it{c})}", "#frac{1}{N_{trig}} dN/d#it{p}_{T, ch jet}^{reco}");
  PRecoilJetPt.SetLegend(0.6, 0.9, 0.70, 0.95, true);
  PRecoilJetPt.New(hRecoilJetPtSignal,Form("TT{%.0f,%.0f} GeV/#it{c}", 
                                           exclusiveTriggerParticleSelection.getPtMinSignal(),
                                           exclusiveTriggerParticleSelection.getPtMaxSignal()),1,1,kRed+1,"hist");
  PRecoilJetPt.New(hRecoilJetPtReference,Form("TT{%.0f,%.0f} GeV/#it{c}", 
                                                       exclusiveTriggerParticleSelection.getPtMinReference(),
                                                       exclusiveTriggerParticleSelection.getPtMaxReference()),1,1,kBlue+1,"hist");
  PRecoilJetPt.NewLatex(0.9,0.9,"ALICE work-in-progress"); 
  PRecoilJetPt.Plot(Form("%s/RecoilJetPt.%s", outputDir.Data(), suffix),false,true);

  // plot ratio of signal to reference
  TH1F* hRecoilJetPtRatio = (TH1F*)hRecoilJetPtSignal->Clone("hRecoilJetPtRatio");
  hRecoilJetPtRatio->Divide(hRecoilJetPtSignal,hRecoilJetPtReference);
  hRecoilJetPtRatio->SetName("hRecoilJetPtRatio");

  Plotting1D PRecoilJetPtRatio;
  hRecoilJetPtRatio->GetYaxis()->SetRangeUser(0.01,10);
  PRecoilJetPtRatio.New(hRecoilJetPtRatio);
  PRecoilJetPtRatio.SetAxisLabel("#bf{#it{p}_{T, ch jet}^{reco} (GeV/#it{c})}", "sig / ref");
  PRecoilJetPtRatio.Plot(Form("%s/RecoilJetPtRatio.%s", outputDir.Data(), suffix),false,false);


  // plot rho distribution for signal and reference
  if (!dGammaSignalQA || !dGammaReferenceQA){
    printf("Error: QA histograms not found\n");
    printf("dGammaSignalQA: %p\n", dGammaSignalQA);
    printf("dGammaReferenceQA: %p\n", dGammaReferenceQA);
  }
    TH2F* hIsoGammaRhoPtSignal = (TH2F*)dGammaSignalQA->Get("hIsoGammaRhoPt");
    TH2F* hIsoGammaRhoPtReference = (TH2F*)dGammaReferenceQA->Get("hIsoGammaRhoPt");
    hIsoGammaRhoPtSignal->Scale(1./hGammaPtSignal->Integral());
    hIsoGammaRhoPtReference->Scale(1./hGammaPtReference->Integral());
    // project onto x axis
  TH1F* hIsoGammaRhoPtSignalProjection = (TH1F*)hIsoGammaRhoPtSignal->ProjectionX("hIsoGammaRhoPtSignalProjection");
  TH1F* hIsoGammaRhoPtReferenceProjection = (TH1F*)hIsoGammaRhoPtReference->ProjectionX("hIsoGammaRhoPtReferenceProjection");
  AvoidLogYconflict(hIsoGammaRhoPtSignalProjection);
  AvoidLogYconflict(hIsoGammaRhoPtReferenceProjection);
  Plotting1D PRhoDistribution;
  PRhoDistribution.New(hIsoGammaRhoPtSignalProjection,Form("TT{%.0f,%.0f} GeV/#it{c}", 
                                              exclusiveTriggerParticleSelection.getPtMinSignal(),
                                              exclusiveTriggerParticleSelection.getPtMaxSignal()),1,1,kRed+1,"hist");
  PRhoDistribution.New(hIsoGammaRhoPtReferenceProjection,Form("TT{%.0f,%.0f} GeV/#it{c}", 
                                              exclusiveTriggerParticleSelection.getPtMinReference(),
                                              exclusiveTriggerParticleSelection.getPtMaxReference()),1,1,kBlue+1,"hist");
    PRhoDistribution.Plot(Form("%s/RhoDistribution.%s", outputDir.Data(), suffix),false,true);

  // plot ratio of signal to reference
  TH1F* hIsoGammaRhoPtRatio = (TH1F*)hIsoGammaRhoPtSignalProjection->Clone("hIsoGammaRhoPtRatio");
  hIsoGammaRhoPtRatio->Divide(hIsoGammaRhoPtSignalProjection,hIsoGammaRhoPtReferenceProjection);
  hIsoGammaRhoPtRatio->SetName("hIsoGammaRhoPtRatio");

  Plotting1D PRhoPtRatio;
  hIsoGammaRhoPtRatio->GetYaxis()->SetRangeUser(0.01,10);
  PRhoPtRatio.New(hIsoGammaRhoPtRatio);
  PRhoPtRatio.SetAxisLabel("#bf{#it{p}_{T}^{#gamma} (GeV/#it{c})}", "sig / ref");
  PRhoPtRatio.Plot(Form("%s/RhoPtRatio.%s", outputDir.Data(), suffix),false,false);
}

void plotIsoGammaJetCorrelations(TDirectory *dGammaJetCorrelations, TString GammaType, GlobalOptions optns)
{
  TString outputDir = Form("%s/%sJetCorrelations", optns.analysisDirPath.Data(), GammaType.Data());
  const char* gammaLatex = GammaType.Contains("Pi0") ? "#pi^{0}" : "#gamma";

  createDirectory(outputDir.Data());

  TH2F *hpTImbalancevsDeltaPhi = (TH2F *)dGammaJetCorrelations->Get("hpTImbalancevsDeltaPhi")->Clone("hpTImbalancevsDeltaPhi");
  hpTImbalancevsDeltaPhi->Rebin2D(2, 2);

  Plotting2D P2DC; // Plot 2D Correlation
  P2DC.New(hpTImbalancevsDeltaPhi);
  P2DC.SetMargins(0.15, 0.1, 0.025, 0.15, 2000, 1750);
  P2DC.SetAxisLabel(Form("#bf{#it{p}_{T}^{jet}/#it{p}_{T}^{%s}}",gammaLatex), Form("#bf{#Delta#phi = |#phi_{%s}-#phi_{jet}|}",gammaLatex));
  P2DC.Plot(Form("%s/pTImbalancevsDeltaPhi.%s", outputDir.Data(), suffix), 0, 0, 0);

  const int NProjections = 8;

  TH1F *hpTImbalance[NProjections];
  TH1F *hDeltaPhi[NProjections];
  int rebinFactorspTImbalance[NProjections] = {1, 1, 1, 1, 1, 1, 1, 1};
  int rebinFactorsDeltaPhi[NProjections] = {1, 1, 1, 1, 1, 1, 1, 1};

  float imbalanceMax = hpTImbalancevsDeltaPhi->GetXaxis()->GetXmax();
  float deltaPhiMax = TMath::Pi();

  PlottingGrid PImbalanceGrid;
  PlottingGrid PDeltaPhiGrid;

  for (int iProj = 0; iProj < NProjections; iProj++)
  {
    float ProjFracFrom = ((float)iProj) / ((float)NProjections);
    float ProjFracTo = ((float)iProj + 1) / ((float)NProjections);
    int ProjBinFromX = hpTImbalancevsDeltaPhi->GetXaxis()->FindBin(ProjFracFrom * imbalanceMax);
    int ProjBinToX = hpTImbalancevsDeltaPhi->GetXaxis()->FindBin(ProjFracTo * imbalanceMax) - 1;
    int ProjBinFromY = hpTImbalancevsDeltaPhi->GetYaxis()->FindBin(ProjFracFrom * deltaPhiMax);
    int ProjBinToY = hpTImbalancevsDeltaPhi->GetYaxis()->FindBin(ProjFracTo * deltaPhiMax) - 1;
    INFO(Form("Correlation projection %d/%d: %d - %d | %d - %d", iProj, NProjections, ProjBinFromX, ProjBinToX, ProjBinFromY, ProjBinToY))

    hpTImbalance[iProj] = (TH1F *)hpTImbalancevsDeltaPhi->ProjectionX(Form("hpTImbalanceProjection_%d", iProj), ProjBinFromY, ProjBinToY);
    hpTImbalance[iProj]->Rebin(rebinFactorspTImbalance[iProj]);
    PImbalanceGrid.New(hpTImbalance[iProj]);
    PImbalanceGrid.NextPad(Form("%d/%d#pi < #Delta#phi < %d/%d#pi", iProj, NProjections, iProj + 1, NProjections));

    hDeltaPhi[iProj] = (TH1F *)hpTImbalancevsDeltaPhi->ProjectionY(Form("hDeltaPhiProjection_%d", iProj), ProjBinFromX, ProjBinToX);
    hDeltaPhi[iProj]->Rebin(rebinFactorsDeltaPhi[iProj]);
    PDeltaPhiGrid.New(hDeltaPhi[iProj]);
    PDeltaPhiGrid.NextPad(Form("%.1f < #it{p}_{T}^{jet}/#it{p}_{T}^{#gamma} < %.1f", ProjFracFrom * imbalanceMax, ProjFracTo * imbalanceMax));
  }

  PImbalanceGrid.Plot(Form("%s/pTImbalanceGrid.%s", outputDir.Data(), suffix));
  PDeltaPhiGrid.Plot(Form("%s/deltaPhiGrid.%s", outputDir.Data(), suffix));
  

  // do correlation plots in bins of photon pt and jet pt
  THnSparseF* hIsoGammaJetDeltaPhiJetPtGammaPt = (THnSparseF*)dGammaJetCorrelations->Get("hIsoGammaJetDeltaPhiJetPtGammaPt");

  TH1F* hJetGammaDeltaPhiProjections[nPtBinsCoarse][nPtBinsCoarse];
  for (int iJetPtBin = 0; iJetPtBin < nPtBinsCoarse; iJetPtBin++)
  {
    for (int iGammaPtBin = 0; iGammaPtBin < nPtBinsCoarse; iGammaPtBin++)
    {
      int minBinJetPt = hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(1)->FindBin(ptBinsCoarse[iJetPtBin]);
      int maxBinJetPt = hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(1)->FindBin(ptBinsCoarse[iJetPtBin + 1]);
      int minBinGammaPt = hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(2)->FindBin(ptBinsCoarse[iGammaPtBin]);
      int maxBinGammaPt = hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(2)->FindBin(ptBinsCoarse[iGammaPtBin + 1]);
      
      hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(1)->SetRange(minBinJetPt, maxBinJetPt);
      hIsoGammaJetDeltaPhiJetPtGammaPt->GetAxis(2)->SetRange(minBinGammaPt, maxBinGammaPt);

      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin] = (TH1F*)hIsoGammaJetDeltaPhiJetPtGammaPt->Projection(0);
      
    }
  }

  TCanvas* cJetGammaDeltaPhi = new TCanvas("cJetGammaDeltaPhi", "cJetGammaDeltaPhi", 800, 800);
  cJetGammaDeltaPhi->Divide(nPtBinsCoarse, nPtBinsCoarse,0,0);
  for (int iJetPtBin = 0; iJetPtBin < nPtBinsCoarse; iJetPtBin++)
  {
    for (int iGammaPtBin = 0; iGammaPtBin < nPtBinsCoarse; iGammaPtBin++)
    {
      cJetGammaDeltaPhi->cd(iJetPtBin * nPtBinsCoarse + iGammaPtBin + 1);
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->Rebin(4);
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->SetTitle(Form("%.1f < #it{p}_{T}^{jet} < %.1f, %.1f < #it{p}_{T}^{#gamma} < %.1f", ptBinsCoarse[iJetPtBin], ptBinsCoarse[iJetPtBin + 1], ptBinsCoarse[iGammaPtBin], ptBinsCoarse[iGammaPtBin + 1]));
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->GetXaxis()->SetTitle("#Delta#phi");
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->GetYaxis()->SetTitle("dN/d#Delta#phi");
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->SetLineColor(kBlack);
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->SetMarkerColor(kBlack);
      hJetGammaDeltaPhiProjections[iJetPtBin][iGammaPtBin]->Draw("hist e");
    }
  }
  cJetGammaDeltaPhi->SaveAs(Form("%s/JetGammaDeltaPhiPtProjections.%s", outputDir.Data(), suffix));

  // do same correlation plot but for 0 to 2pi range
  // do correlation plots in bins of photon pt and jet pt
  THnSparseF* hIsoGammaJetDeltaPhi2piJetPtGammaPt = (THnSparseF*)dGammaJetCorrelations->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");

  TH1F* hJetGammaDeltaPhi2piProjections[nPtBinsCoarse][nPtBinsCoarse];
  TGraphPolar* gJetGammaDeltaPhi2piProjections[nPtBinsCoarse][nPtBinsCoarse];
  for (int iJetPtBin = 0; iJetPtBin < nPtBinsCoarse; iJetPtBin++)
  {
    for (int iGammaPtBin = 0; iGammaPtBin < nPtBinsCoarse; iGammaPtBin++)
    {
      int minBinJetPt = hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(1)->FindBin(ptBinsCoarse[iJetPtBin]);
      int maxBinJetPt = hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(1)->FindBin(ptBinsCoarse[iJetPtBin + 1]);
      int minBinGammaPt = hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(2)->FindBin(ptBinsCoarse[iGammaPtBin]);
      int maxBinGammaPt = hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(2)->FindBin(ptBinsCoarse[iGammaPtBin + 1]);
      
      hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(1)->SetRange(minBinJetPt, maxBinJetPt);
      hIsoGammaJetDeltaPhi2piJetPtGammaPt->GetAxis(2)->SetRange(minBinGammaPt, maxBinGammaPt);

      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin] = (TH1F*)hIsoGammaJetDeltaPhi2piJetPtGammaPt->Projection(0);
      
      
    }
  }

  TCanvas* cJetGammaDeltaPhi2pi = new TCanvas("cJetGammaDeltaPhi2pi", "cJetGammaDeltaPhi2pi", 800, 800);
  cJetGammaDeltaPhi2pi->Divide(nPtBinsCoarse, nPtBinsCoarse,0,0);
  TLine* piLineDashed[nPtBinsCoarse][nPtBinsCoarse];
  for (int iJetPtBin = 0; iJetPtBin < nPtBinsCoarse; iJetPtBin++)
  {
    for (int iGammaPtBin = 0; iGammaPtBin < nPtBinsCoarse; iGammaPtBin++)
    {
      cJetGammaDeltaPhi2pi->cd(iJetPtBin * nPtBinsCoarse + iGammaPtBin + 1);
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->Rebin(4);
            gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin] = TH1ToTGraphPolar(hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin],Form("gJetGammaDeltaPhi2piProjections_%d_%d",iJetPtBin,iGammaPtBin));
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetTitle(Form("%.1f < #it{p}_{T}^{jet} < %.1f, %.1f < #it{p}_{T}^{#gamma} < %.1f", ptBinsCoarse[iJetPtBin], ptBinsCoarse[iJetPtBin + 1], ptBinsCoarse[iGammaPtBin], ptBinsCoarse[iGammaPtBin + 1]));
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetXaxis()->SetTitle("#Delta#phi");
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetYaxis()->SetTitle("dN/d#Delta#phi");
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetLineColor(kBlack);
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetMarkerColor(kBlack);
      hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->Draw("hist e");

      // draw a dashed gray line at pi
      piLineDashed[iJetPtBin][iGammaPtBin] = new TLine(TMath::Pi(), hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetMinimum(), TMath::Pi(), hJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetMaximum());
      piLineDashed[iJetPtBin][iGammaPtBin]->SetLineStyle(2);
      piLineDashed[iJetPtBin][iGammaPtBin]->SetLineColor(kGray);
      piLineDashed[iJetPtBin][iGammaPtBin]->Draw("same");
    }
  }
  cJetGammaDeltaPhi2pi->SaveAs(Form("%s/JetGammaDeltaPhiPtProjections_2pi.%s", outputDir.Data(), suffix));

  TCanvas* cJetGammaDeltaPhi2piPolar = new TCanvas("cJetGammaDeltaPhi2piPolar", "cJetGammaDeltaPhi2piPolar", 800, 800);
  cJetGammaDeltaPhi2piPolar->Divide(nPtBinsCoarse, nPtBinsCoarse,0,0);

  for (int iJetPtBin = 0; iJetPtBin < nPtBinsCoarse; iJetPtBin++)
  {
    for (int iGammaPtBin = 0; iGammaPtBin < nPtBinsCoarse; iGammaPtBin++)
    {
      cJetGammaDeltaPhi2piPolar->cd(iJetPtBin * nPtBinsCoarse + iGammaPtBin + 1);
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetTitle(Form("%.1f < #it{p}_{T}^{jet} < %.1f, %.1f < #it{p}_{T}^{#gamma} < %.1f", ptBinsCoarse[iJetPtBin], ptBinsCoarse[iJetPtBin + 1], ptBinsCoarse[iGammaPtBin], ptBinsCoarse[iGammaPtBin + 1]));
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetXaxis()->SetLabelSize(0.05);
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetLineColor(kRed);
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetMarkerSize(0);
      // make less divisions for y axis
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetYaxis()->SetNdivisions(400);
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->SetFillColorAlpha(kRed, 0.5);
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->Draw("ALF");
      cJetGammaDeltaPhi2piPolar->Update();
      gPad->Update();
      gJetGammaDeltaPhi2piProjections[iJetPtBin][iGammaPtBin]->GetPolargram()->SetToRadian();
    }
  }

  cJetGammaDeltaPhi2piPolar->SaveAs(Form("%s/JetGammaDeltaPhiPtProjections_2piPolar.%s", outputDir.Data(), suffix));
  // now do the same but for gJetGammaDeltaPhi2piProjections

  // try only to plot one bin
  TCanvas* cJetGammaDeltaPhi2piPolarExampleBin = new TCanvas("cJetGammaDeltaPhi2piPolarExampleBin", "cJetGammaDeltaPhi2piPolarExampleBin", 800, 800);

  cJetGammaDeltaPhi2piPolarExampleBin->cd();
  gJetGammaDeltaPhi2piProjections[0][0]->SetTitle(Form("%.1f < #it{p}_{T}^{jet} < %.1f, %.1f < #it{p}_{T}^{#gamma} < %.1f", ptBinsCoarse[0], ptBinsCoarse[0 + 1], ptBinsCoarse[0], ptBinsCoarse[0 + 1]));
  gJetGammaDeltaPhi2piProjections[0][0]->SetLineColor(kRed);
  gJetGammaDeltaPhi2piProjections[0][0]->SetMarkerColor(kRed);
  gJetGammaDeltaPhi2piProjections[0][0]->SetFillColorAlpha(kRed, 0.5);
  gJetGammaDeltaPhi2piProjections[0][0]->Draw("ALF");
  cJetGammaDeltaPhi2piPolarExampleBin->Update();
  gJetGammaDeltaPhi2piProjections[0][0]->GetPolargram()->SetToRadian();
  cJetGammaDeltaPhi2piPolarExampleBin->SaveAs(Form("%s/JetGammaDeltaPhiPtProjections_2piPolarExampleBin.%s", outputDir.Data(), suffix));


}

void plotIsoGammaJetCorrelationQA(TDirectory *dGammaJetCorrelationQA, TString GammaType, GlobalOptions optns)
{
  TString outputDir = Form("%s/%sJetCorrelations", optns.analysisDirPath.Data(), GammaType.Data());

  createDirectory(outputDir.Data());

  TH2F *hIsoGammaFromCorrEtaPhi = (TH2F *)dGammaJetCorrelationQA->Get("hIsoGammaFromCorrelationEtaPhiMap");
  Plotting2D PIEtaPhiMap;
  PIEtaPhiMap.SetMargins(0.12, 0.12, 0.05, 0.125);
  PIEtaPhiMap.New(hIsoGammaFromCorrEtaPhi);
  PIEtaPhiMap.AddEMCalOutline();
  PIEtaPhiMap.AddPHOSOutline();
  PIEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PIEtaPhiMap.Plot(Form("%s/IsoGammaEtaPhiMap.%s", outputDir.Data(), suffix));

  TH2F *hJetFromCorrEtaPhi = (TH2F *)dGammaJetCorrelationQA->Get("hJetFromCorrelationEtaPhiMap");
  Plotting2D PJEtaPhiMap;
  PJEtaPhiMap.SetMargins(0.12, 0.12, 0.05, 0.125);
  PJEtaPhiMap.New(hJetFromCorrEtaPhi);
  PJEtaPhiMap.AddEMCalOutline();
  PJEtaPhiMap.AddPHOSOutline();
  PJEtaPhiMap.SetAxisLabel("#bf{#eta}", "#bf{#phi}");
  PJEtaPhiMap.Plot(Form("%s/JetEtaPhiMap.%s", outputDir.Data(), suffix));
}

void plotHistosFromTree(TString AnalysisDirectory, bool isDebugRun = false)
{

  Float_t r[4][2] = {{77, 228}, {115, 165}, {0, 72}, {134, 173}};
    Float_t g[4][2] = {{75, 227}, {124, 171}, {97, 143}, {0, 59}};
    Float_t b[4][2] = {{70, 221}, {69, 82}, {143, 218}, {71, 118}};
    for (Int_t iC = 0; iC < 4; iC++)
    {
        for (Int_t iH = 0; iH < 2; iH++)
        {
            Int_t ColorCode = TColor::GetFreeColorIndex();
            TColor *color = new TColor(ColorCode, r[iC][iH] / 255., g[iC][iH] / 255., b[iC][iH] / 255.);
            GoetheColor[iC][iH] = ColorCode;
        }
    }

  GlobalOptions optns(AnalysisDirectory, isDebugRun);

  TString inputFilePath = Form("%s/HistosFromTree.root", AnalysisDirectory.Data());

  if (!std::filesystem::exists(inputFilePath.Data()))
  {
    inputFilePath = Form("%s/HistosFromTree_0.root", AnalysisDirectory.Data());
    INFO(Form("Did not find HistosFromTree.root, looking for %s now.", inputFilePath.Data()))
  }

  TFile *fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  // -------------------------
  //        Event Plots
  // -------------------------
  TDirectory *dEventQA = (TDirectory *)fIn->Get("EventQA");
  if (!dEventQA)
    FATAL("Dir EventQA not found")
  plotEventQA(dEventQA, optns);
  
  // -------------------------
  //        Jet Plots
  // -------------------------
  TDirectory *dJets = (TDirectory *)fIn->Get("Jets");   
  if(dJets) plotJets(dJets, optns);
  else WARN("Dir Jets not found. Skipping");
  if (optns.doQA)
  {
    TDirectory *dJetQA = (TDirectory *)fIn->Get("JetQA");
    TDirectory *dJetQARaw = (TDirectory *)fIn->Get("JetQARaw");
    if(dJetQA) plotJetQA(dJetQA, "Jets" ,optns);
    else WARN("Dir JetQA not found. Skipping");
    if(dJetQARaw) plotJetQA(dJetQARaw, "JetsRaw" ,optns);
    else WARN("Dir JetQARaw not found. Skipping");

  }
  
  // ------------------------------------------
  //        Cluster / Gamma / Iso Gamma Plots
  // ------------------------------------------
  TDirectory *dClusters = (TDirectory *)fIn->Get("Clusters");
  TDirectory *dGammas = (TDirectory *)fIn->Get("Gammas");
  TDirectory *dIsoGammas = (TDirectory *)fIn->Get("IsoGammas");
  TDirectory *dTriggerPhotonsSignal = (TDirectory *)fIn->Get("TriggerPhotonsSignal");
  TDirectory *dTriggerPhotonsReference = (TDirectory *)fIn->Get("TriggerPhotonsReference");
  if(dClusters)  plotGammas(dClusters, "Clusters", optns);
  else WARN("Dir Clusters not found. Skipping");
  if(dGammas) plotGammas(dGammas, "Gammas", optns);
  else WARN("Dir Gammas not found. Skipping");
  if(dIsoGammas) plotGammas(dIsoGammas, "IsoGammas", optns);
  else WARN("Dir IsoGammas not found. Skipping");
  if(dTriggerPhotonsSignal) plotGammas(dTriggerPhotonsSignal, "TriggerPhotonsSignal", optns);
  else WARN("Dir TriggerPhotonsSignal not found. Skipping");
  if(dTriggerPhotonsReference) plotGammas(dTriggerPhotonsReference, "TriggerPhotonsReference", optns);
  else WARN("Dir TriggerPhotonsReference not found. Skipping");

  if (optns.doQA)
  {
    TDirectory *dClusterQA = (TDirectory *)fIn->Get("ClusterQA");
    TDirectory *dGammaQA = (TDirectory *)fIn->Get("GammaQA");
    TDirectory *dIsoGammaQA = (TDirectory *)fIn->Get("IsoGammaQA");
    TDirectory *dTriggerPhotonSignalQA = (TDirectory *)fIn->Get("TriggerPhotonSignalQA");
    TDirectory *dTriggerPhotonReferenceQA = (TDirectory *)fIn->Get("TriggerPhotonReferenceQA");

    if(dClusterQA) plotGammaQA(dClusterQA, "Clusters", optns);
    else WARN("Dir ClusterQA not found. Skipping");
    if(dGammaQA) plotGammaQA(dGammaQA, "Gammas", optns);
    else WARN("Dir GammaQA not found. Skipping");
    if(dIsoGammaQA) plotGammaQA(dIsoGammaQA, "IsoGammas", optns);
    else WARN("Dir IsoGammaQA not found. Skipping");
    if(dTriggerPhotonSignalQA) plotGammaQA(dTriggerPhotonSignalQA, "TriggerPhotonsSignal", optns);
    else WARN("Dir TriggerPhotonSignalQA not found. Skipping");
    if(dTriggerPhotonReferenceQA) plotGammaQA(dTriggerPhotonReferenceQA, "TriggerPhotonsReference", optns);
    else WARN("Dir TriggerPhotonReferenceQA not found. Skipping");
  }
  // -------------------------------
  //        Merged Pi0 exclusive selection
  // -------------------------------
  TDirectory *dTriggerMergedPi0sSignal = (TDirectory *)fIn->Get("TriggerMergedPi0sSignal");
  TDirectory *dTriggerMergedPi0sReference = (TDirectory *)fIn->Get("TriggerMergedPi0sReference");
  if(dTriggerMergedPi0sSignal) plotGammas(dTriggerMergedPi0sSignal, "TriggerMergedPi0sSignal", optns);
  else WARN("Dir TriggerMergedPi0sSignal not found. Skipping");
  if(dTriggerMergedPi0sReference) plotGammas(dTriggerMergedPi0sReference, "TriggerMergedPi0sReference", optns);
  else WARN("Dir TriggerMergedPi0sReference not found. Skipping");

  if (optns.doQA)
  {
    TDirectory *dTriggerMergedPi0sSignalQA = (TDirectory *)fIn->Get("TriggerMergedPi0sSignalQA");
    TDirectory *dTriggerMergedPi0sReferenceQA = (TDirectory *)fIn->Get("TriggerMergedPi0sReferenceQA");
    if(dTriggerMergedPi0sSignalQA) plotGammaQA(dTriggerMergedPi0sSignalQA, "TriggerMergedPi0sSignal", optns);
    else WARN("Dir TriggerMergedPi0sSignalQA not found. Skipping");
    if(dTriggerMergedPi0sReferenceQA) plotGammaQA(dTriggerMergedPi0sReferenceQA, "TriggerMergedPi0sReference", optns);
    else WARN("Dir TriggerMergedPi0sReferenceQA not found. Skipping");
  }

  // -------------------------------
  //        Gamma-Jet Correlations
  // -------------------------------
  TDirectory *dGammaJetCorrelations = (TDirectory *)fIn->Get("GammaJetCorrelations");
  if(dGammaJetCorrelations) plotIsoGammaJetCorrelations(dGammaJetCorrelations, "Gamma", optns);
  else WARN("Dir GammaJetCorrelations not found. Skipping");

  if (optns.doQA)
  {
    TDirectory *dGammaJetCorrelationQA = (TDirectory *)fIn->Get("GammaJetCorrelationQA");
    if(dGammaJetCorrelationQA) plotIsoGammaJetCorrelationQA(dGammaJetCorrelationQA, "Gamma", optns);
    else WARN("Dir GammaJetCorrelationQA not found. Skipping");
  }
  
  //------------------------------------
  //       merged Pi0-Jet Correlations
  //------------------------------------
  TDirectory *dmPi0JetCorrelations = (TDirectory *)fIn->Get("mPi0JetCorrelations");
  if (dmPi0JetCorrelations) plotIsoGammaJetCorrelations(dmPi0JetCorrelations, "mPi0", optns);
  else WARN("Dir mPi0JetCorrelations not found. Skipping");

  if (optns.doQA)
  {
    TDirectory *dmPi0JetCorrelationQA = (TDirectory *)fIn->Get("mPi0JetCorrelationQA");
    if (dmPi0JetCorrelationQA) plotIsoGammaJetCorrelationQA(dmPi0JetCorrelationQA, "mPi0", optns);
    else WARN("Dir mPi0JetCorrelationQA not found. Skipping");
  }

  //------------------------------------
  //       Pi0-Jet Correlations
  //------------------------------------
  TDirectory *dGGPi0JetCorrelations = (TDirectory *)fIn->Get("GGPi0JetCorrelations");
  if(dGGPi0JetCorrelations) plotIsoGammaJetCorrelations(dGGPi0JetCorrelations, "GGPi0", optns);
  else WARN("Dir GGPi0JetCorrelations not found. Skipping");

  if (optns.doQA)
  {
    TDirectory *dGGPi0JetCorrelationQA = (TDirectory *)fIn->Get("GGPi0JetCorrelationQA");
    if(dGGPi0JetCorrelationQA) plotIsoGammaJetCorrelationQA(dGGPi0JetCorrelationQA, "GGPi0", optns);
    else WARN("Dir GGPi0JetCorrelationQA not found. Skipping");
  }

  //----------------------------------------------------------------
  //       Trigger Photon - Jet Correlations (Exclusive selections)
  //----------------------------------------------------------------
  TDirectory *dTriggerGammaJetCorrelationsSignal = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationsSignal");
  TDirectory *dTriggerGammaJetCorrelationsReference = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationsReference");
  TDirectory *dTriggerGammaJetCorrelationsSignalQA = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationSignalQA");
  TDirectory *dTriggerGammaJetCorrelationsReferenceQA = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationReferenceQA");
  TDirectory *dTriggerPhotonSignalQA = (TDirectory *)fIn->Get("TriggerPhotonSignalQA");
  TDirectory *dTriggerPhotonReferenceQA = (TDirectory *)fIn->Get("TriggerPhotonReferenceQA");
  if(dTriggerGammaJetCorrelationsSignal) plotIsoGammaJetCorrelations(dTriggerGammaJetCorrelationsSignal, "TriggerPhotonsSignal", optns);
  else WARN("Dir TriggerGammaJetCorrelationsSignal not found. Skipping");
  if(dTriggerGammaJetCorrelationsReference) plotIsoGammaJetCorrelations(dTriggerGammaJetCorrelationsReference, "TriggerPhotonsReference", optns);
  else WARN("Dir TriggerGammaJetCorrelationsReference not found. Skipping");
  if(dTriggerGammaJetCorrelationsSignalQA) plotIsoGammaJetCorrelationQA(dTriggerGammaJetCorrelationsSignalQA, "TriggerPhotonsSignal", optns);
  else WARN("Dir TriggerGammaJetCorrelationSignalQA not found. Skipping");
  if(dTriggerGammaJetCorrelationsReferenceQA) plotIsoGammaJetCorrelationQA(dTriggerGammaJetCorrelationsReferenceQA, "TriggerPhotonsReference", optns);
  else WARN("Dir TriggerGammaJetCorrelationReferenceQA not found. Skipping");
  if(dTriggerGammaJetCorrelationsSignal && dTriggerGammaJetCorrelationsReference && dTriggerPhotonSignalQA && dTriggerPhotonReferenceQA) plotExclusiveSelectionGammaJetCorrelations(dTriggerGammaJetCorrelationsSignal, dTriggerGammaJetCorrelationsReference, dTriggerPhotonsSignal, dTriggerPhotonSignalQA, dTriggerPhotonsReference, dTriggerPhotonReferenceQA, "TriggerPhotonCorrelationsCombined", optns);
  else WARN("Dir TriggerGammaJetCorrelationsSignal or TriggerGammaJetCorrelationsReference not found. Skipping");

  //----------------------------------------------------------------
  //       Trigger Merged Pi0 - Jet Correlations (Exclusive selections)
  //----------------------------------------------------------------
  TDirectory *dTriggerMergedPi0JetCorrelationsSignal = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationsSignal");
  TDirectory *dTriggerMergedPi0JetCorrelationsReference = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationsReference");
  TDirectory *dTriggerMergedPi0JetCorrelationsSignalQA = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationSignalQA");
  TDirectory *dTriggerMergedPi0JetCorrelationsReferenceQA = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationReferenceQA");
  TDirectory *dTriggerMergedPi0sSignalQA = (TDirectory *)fIn->Get("TriggerMergedPi0sSignalQA");
  TDirectory *dTriggerMergedPi0sReferenceQA = (TDirectory *)fIn->Get("TriggerMergedPi0sReferenceQA");
  if(dTriggerMergedPi0JetCorrelationsSignal) plotIsoGammaJetCorrelations(dTriggerMergedPi0JetCorrelationsSignal, "TriggerMergedPi0sSignal", optns);
  else WARN("Dir TriggerMergedPi0JetCorrelationsSignal not found. Skipping");
  if(dTriggerMergedPi0JetCorrelationsReference) plotIsoGammaJetCorrelations(dTriggerMergedPi0JetCorrelationsReference, "TriggerMergedPi0sReference", optns);
  else WARN("Dir TriggerMergedPi0JetCorrelationsReference not found. Skipping");
  if(dTriggerMergedPi0JetCorrelationsSignalQA) plotIsoGammaJetCorrelationQA(dTriggerMergedPi0JetCorrelationsSignalQA, "TriggerMergedPi0sSignal", optns);
  else WARN("Dir TriggerMergedPi0JetCorrelationSignalQA not found. Skipping");
  if(dTriggerMergedPi0JetCorrelationsReferenceQA) plotIsoGammaJetCorrelationQA(dTriggerMergedPi0JetCorrelationsReferenceQA, "TriggerMergedPi0sReference", optns);
  else WARN("Dir TriggerMergedPi0JetCorrelationReferenceQA not found. Skipping");
  if(dTriggerMergedPi0JetCorrelationsSignal && dTriggerMergedPi0JetCorrelationsReference && dTriggerMergedPi0JetCorrelationsSignalQA && dTriggerMergedPi0JetCorrelationsReferenceQA) plotExclusiveSelectionGammaJetCorrelations(dTriggerMergedPi0JetCorrelationsSignal, dTriggerMergedPi0JetCorrelationsReference, dTriggerMergedPi0sSignal, dTriggerMergedPi0sSignalQA, dTriggerMergedPi0sReference, dTriggerMergedPi0sReferenceQA, "TriggerMergedPi0CorrelationsCombined", optns);
  else WARN("Dir TriggerMergedPi0JetCorrelationsSignal or TriggerMergedPi0JetCorrelationsReference not found. Skipping");

  //------------------------------------
  //       Good Night! zzz
  //------------------------------------
  fIn->Close();
}