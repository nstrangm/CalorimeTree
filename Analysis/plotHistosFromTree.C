#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <filesystem>
#include "Utilities.h"
#include "Cuts.h"
#include "plotHistosFromTree.h"

const char* suffix = "png";

void createDirectory(const char* path)
{
  struct stat info;
  if (stat(path, &info) != 0) {
    if (mkdir(path, 0755) == -1) // Directory does not exist, create it
      std::cerr << "Error creating directory: " << path << std::endl;
  }
  return;
}
//Add to plot class?
void AvoidLogYconflict(TH1* h, float logymin=1e-10)
{
  double min_nonzero = std::numeric_limits<double>::max();
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
      float content = h->GetBinContent(i);
      if (content > 1e-10 && content < min_nonzero) {
          min_nonzero = content;
      }
      if (h->GetBinContent(i) == 0) {
          h->SetBinContent(i, logymin);
      }
  }
  h->SetMinimum(min_nonzero);
  return;
}

void SetTitleAndLabels(TH1* h, const char* title)
{
  h->SetTitle(title);
  return;
}

void plotJets(TDirectory* dJets, GlobalOptions optns)
{
  ENTER

  const char* outputDir = Form("%s/Jets", optns.analysisDirPath.Data());
  createDirectory(outputDir);

  if (optns.isMC) {
    // Plot response matrix
    TH2F* hDLJetPtVsPLJetPt = (TH2F*)dJets->Get("hDLJetPtVsPLJetPt");
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
    TH1F* hDLJetPtVsPLJetPtProjection[NProjections];
    TGraphAsymmErrors* gRatioDLJetPtToPLJetPtProjection[NProjections];
    Plotting1D PResponseMatrixProjections;
    for (int iProj = 0; iProj < NProjections; iProj++) {
      int binProjFrom = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][0] + 0.001);
      int binProjTo = hDLJetPtVsPLJetPt->GetYaxis()->FindBin(projectionRanges[iProj][1] + 0.001);
      hDLJetPtVsPLJetPtProjection[iProj] = (TH1F*)hDLJetPtVsPLJetPt->ProjectionX(Form("hDLJetPtVsPLJetPtProjection_%d", iProj), binProjFrom, binProjTo);
      hDLJetPtVsPLJetPtProjection[iProj]->Rebin(rebinFactors[iProj]);
      gRatioDLJetPtToPLJetPtProjection[iProj] = new TGraphAsymmErrors();
      float PLJetpT = (projectionRanges[iProj][0] + projectionRanges[iProj][1]) / 2.;
      for (int iDLJetPtBin = 0; iDLJetPtBin < hDLJetPtVsPLJetPtProjection[iProj]->GetNbinsX(); iDLJetPtBin++) {
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
    PResponseMatrixProjections.SetLegend(0.5, 0.9, 0.5, 0.9);
    PResponseMatrixProjections.SetAxisLabel("#bf{#it{p}_{T}^{DL}/#it{p}_{T}^{PL}} #hat{=} #bf{#it{p}_{T}^{rec}/#it{p}_{T}^{gen}}");
    PResponseMatrixProjections.Plot(Form("%s/ResponseMatrixProjections.%s", outputDir, suffix));
  }

  EXIT
}

void plotJetQA(TDirectory* dJetQA, GlobalOptions optns)
{
  ENTER

  JetCuts jetCuts(optns);

  const char* outputDir = Form("%s/Jets", optns.analysisDirPath.Data());
  createDirectory(outputDir);

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

void plotIsoGammaQA(TDirectory* dIsoGammaQA, GlobalOptions optns)
{
  ENTER

  //const char* outputDir = Form("%s/IsoGammas", optns.analysisDirPath.Data());
  //const char* outputDirChar = Form("%s/IsoGammas", optns.analysisDirPath.Data());
  TString outputDir = Form("%s/IsoGammas", optns.analysisDirPath.Data());

  createDirectory(outputDir.Data());

  // Plot response matrix
  Plotting1D PMinMassDiffToPi0;
  TH1F* hMinMassDiffToPi0Signal;
  TH1F* hMinMassDiffToPi0Background;
  TH1F* hMinMassDiffToPi0;

  Plotting1D PIsoGammaE;
  TH1F* hIsoGammaESignal;
  TH1F* hIsoGammaEBackground;  
  TH1F* hIsoGammaE;

  Plotting1D PIsoGammaM02;
  TH1F* hIsoGammaM02Signal;
  TH1F* hIsoGammaM02Background;
  TH1F* hIsoGammaM02;

  Plotting1D PIsoGammaM20;
  TH1F* hIsoGammaM20Signal;
  TH1F* hIsoGammaM20Background;
  TH1F* hIsoGammaM20;

  Plotting1D PIsoGammaPx;
  TH1F* hIsoGammaPxSignal;
  TH1F* hIsoGammaPxBackground;
  TH1F* hIsoGammaPx;

  Plotting1D PIsoGammaPy;
  TH1F* hIsoGammaPySignal;
  TH1F* hIsoGammaPyBackground;
  TH1F* hIsoGammaPy;

  Plotting1D PIsoGammaPz;
  TH1F* hIsoGammaPzSignal;
  TH1F* hIsoGammaPzBackground;
  TH1F* hIsoGammaPz;

  Plotting1D PIsoGammaIsoCharged;
  TH1F* hIsoGammaIsoChargedSignal;
  TH1F* hIsoGammaIsoChargedBackground;
  TH1F* hIsoGammaIsoCharged;

  Plotting1D PIsoGammaIsoChargedCorrected;
  TH1F* hIsoGammaIsoChargedCorrectedSignal;
  TH1F* hIsoGammaIsoChargedCorrectedBackground;
  TH1F* hIsoGammaIsoChargedCorrected;

  if (optns.isMC) {
    //Diff to pi0 mass
    hMinMassDiffToPi0Signal = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0Signal")->Clone("hMinMassDiffToPi0Signal");
    hMinMassDiffToPi0Background = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0Background")->Clone("hMinMassDiffToPi0Background");
    float integralSignal = hMinMassDiffToPi0Signal->Integral(1, hMinMassDiffToPi0Signal->GetNbinsX());
    float integralBG = hMinMassDiffToPi0Background->Integral(1, hMinMassDiffToPi0Background->GetNbinsX());
    hMinMassDiffToPi0Signal->Scale(1. / integralSignal);
    hMinMassDiffToPi0Background->Scale(1. / integralBG);
    AvoidLogYconflict(hMinMassDiffToPi0Signal);
    AvoidLogYconflict(hMinMassDiffToPi0Background);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Signal, "Signal");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Background, "Background");

    //E
    hIsoGammaESignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaESignal")->Clone("hIsoGammaESignal");
    hIsoGammaEBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaEBackground")->Clone("hIsoGammaEBackground");
    float integralESignal = hIsoGammaESignal->Integral(1, hIsoGammaESignal->GetNbinsX());
    float integralEBG = hIsoGammaEBackground->Integral(1, hIsoGammaEBackground->GetNbinsX());
    hIsoGammaESignal->Scale(1./integralESignal);
    hIsoGammaEBackground->Scale(1./integralEBG);
    AvoidLogYconflict(hIsoGammaESignal);
    AvoidLogYconflict(hIsoGammaEBackground);
    PIsoGammaE.New(hIsoGammaESignal,"Signal");
    PIsoGammaE.New(hIsoGammaEBackground,"Background");

    //M02
    hIsoGammaM02Signal = (TH1F*)dIsoGammaQA->Get("hIsoGammaM02Signal")->Clone("hIsoGammaM02Signal");
    hIsoGammaM02Background = (TH1F*)dIsoGammaQA->Get("hIsoGammaM02Background")->Clone("hIsoGammaM02Background");
    float integralM02Signal = hIsoGammaM02Signal->Integral(1, hIsoGammaM02Signal->GetNbinsX());
    float integralM02BG = hIsoGammaM02Background->Integral(1, hIsoGammaM02Background->GetNbinsX());
    hIsoGammaM02Signal->Scale(1./integralM02Signal);
    hIsoGammaM02Background->Scale(1./integralM02BG);
    AvoidLogYconflict(hIsoGammaM02Signal);
    AvoidLogYconflict(hIsoGammaM02Background);
    PIsoGammaM02.New(hIsoGammaM02Signal,"Signal");
    PIsoGammaM02.New(hIsoGammaM02Background,"Background");

    //M20
    hIsoGammaM20Signal = (TH1F*)dIsoGammaQA->Get("hIsoGammaM20Signal")->Clone("hIsoGammaM20Signal");
    hIsoGammaM20Background = (TH1F*)dIsoGammaQA->Get("hIsoGammaM20Background")->Clone("hIsoGammaM20Background");
    float integralM20Signal = hIsoGammaM20Signal->Integral(1, hIsoGammaM20Signal->GetNbinsX());
    float integralM20BG = hIsoGammaM20Background->Integral(1, hIsoGammaM20Background->GetNbinsX());
    hIsoGammaM20Signal->Scale(1./integralM20Signal);
    hIsoGammaM20Background->Scale(1./integralM20BG);
    AvoidLogYconflict(hIsoGammaM20Signal);
    AvoidLogYconflict(hIsoGammaM20Background);
    PIsoGammaM20.New(hIsoGammaM20Signal,"Signal");
    PIsoGammaM20.New(hIsoGammaM20Background,"Background");

    //Px
    hIsoGammaPxSignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaPxSignal")->Clone("hIsoGammaPxSignal");
    hIsoGammaPxBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaPxBackground")->Clone("hIsoGammaPxBackground");
    float integralPxSignal = hIsoGammaPxSignal->Integral(1, hIsoGammaPxSignal->GetNbinsX());
    float integralPxBG = hIsoGammaPxBackground->Integral(1, hIsoGammaPxBackground->GetNbinsX());
    hIsoGammaPxSignal->Scale(1./integralPxSignal);
    hIsoGammaPxBackground->Scale(1./integralPxBG);
    AvoidLogYconflict(hIsoGammaPxSignal);
    AvoidLogYconflict(hIsoGammaPxBackground);
    PIsoGammaPx.New(hIsoGammaPxSignal,"Signal");
    PIsoGammaPx.New(hIsoGammaPxBackground,"Background");

    //Py
    hIsoGammaPySignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaPySignal")->Clone("hIsoGammaPySignal");
    hIsoGammaPyBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaPyBackground")->Clone("hIsoGammaPyBackground");
    float integralPySignal = hIsoGammaPySignal->Integral(1, hIsoGammaPySignal->GetNbinsX());
    float integralPyBG = hIsoGammaPyBackground->Integral(1, hIsoGammaPyBackground->GetNbinsX());
    hIsoGammaPySignal->Scale(1./integralPySignal);
    hIsoGammaPyBackground->Scale(1./integralPyBG);
    AvoidLogYconflict(hIsoGammaPySignal);
    AvoidLogYconflict(hIsoGammaPyBackground);
    PIsoGammaPy.New(hIsoGammaPySignal,"Signal");
    PIsoGammaPy.New(hIsoGammaPyBackground,"Background");

    //Pz
    hIsoGammaPzSignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaPzSignal")->Clone("hIsoGammaPzSignal");
    hIsoGammaPzBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaPzBackground")->Clone("hIsoGammaPzBackground");
    float integralPzSignal = hIsoGammaPzSignal->Integral(1, hIsoGammaPzSignal->GetNbinsX());
    float integralPzBG = hIsoGammaPzBackground->Integral(1, hIsoGammaPzBackground->GetNbinsX());
    hIsoGammaPzSignal->Scale(1./integralPzSignal);
    hIsoGammaPzBackground->Scale(1./integralPzBG);
    AvoidLogYconflict(hIsoGammaPzSignal);
    AvoidLogYconflict(hIsoGammaPzBackground);
    PIsoGammaPz.New(hIsoGammaPzSignal,"Signal");
    PIsoGammaPz.New(hIsoGammaPzBackground,"Background");

    //IsoCharged
    hIsoGammaIsoChargedSignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaIsoChargedSignal")->Clone("hIsoGammaIsoChargedSignal");
    hIsoGammaIsoChargedBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaIsoChargedBackground")->Clone("hIsoGammaIsoChargedBackground");
    float integralIsolationSignal = hIsoGammaIsoChargedSignal->Integral(1, hIsoGammaIsoChargedSignal->GetNbinsX());
    float integralIsolationBG = hIsoGammaIsoChargedBackground->Integral(1, hIsoGammaIsoChargedBackground->GetNbinsX());
    hIsoGammaIsoChargedSignal->Scale(1./integralIsolationSignal);
    hIsoGammaIsoChargedBackground->Scale(1./integralIsolationBG);
    AvoidLogYconflict(hIsoGammaIsoChargedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedBackground);
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedSignal,"Signal");
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedBackground,"Background");

    //IsoChargedCorrected
    hIsoGammaIsoChargedCorrectedSignal = (TH1F*)dIsoGammaQA->Get("hIsoGammaIsoChargedCorrectedSignal")->Clone("hIsoGammaIsoChargedCorrectedSignal");
    hIsoGammaIsoChargedCorrectedBackground = (TH1F*)dIsoGammaQA->Get("hIsoGammaIsoChargedCorrectedBackground")->Clone("hIsoGammaIsoChargedCorrectedBackground");
    float integralIsolationCorrectedSignal = hIsoGammaIsoChargedCorrectedSignal->Integral(1, hIsoGammaIsoChargedCorrectedSignal->GetNbinsX());
    float integralIsolationCorrectedBG = hIsoGammaIsoChargedCorrectedBackground->Integral(1, hIsoGammaIsoChargedCorrectedBackground->GetNbinsX());
    hIsoGammaIsoChargedCorrectedSignal->Scale(1./integralIsolationSignal);
    hIsoGammaIsoChargedCorrectedBackground->Scale(1./integralIsolationBG);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedBackground);
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedSignal,"Signal");
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedBackground,"Background");

    //M02vsPt (Slicing 2D hist pr. specified pT-bins)
    TH2F* h2M02vspT;
    TH2F* h2M02vspTSignal;
    TH2F* h2M02vspTBackground;
    TH2F* h2M02vspTAllGamma;
    TH2F* h2M02vspTPionDecayGamma;
    TH2F* h2M02vspTEtaDecayGamma;
    TH2F* h2M02vspTMergedPionGamma;

    h2M02vspT = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pT")->Clone("hIsoGammaM02pT");
    h2M02vspTSignal = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTSignal")->Clone("hIsoGammaM02pTSignal");
    h2M02vspTBackground = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTBackground")->Clone("hIsoGammaM02pTBackground");
    h2M02vspTAllGamma = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTAllGamma")->Clone("hIsoGammaM02pTAllGamma");
    h2M02vspTPionDecayGamma = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTPionDecayGamma")->Clone("hIsoGammaM02pTPionDecayGamma");
    h2M02vspTEtaDecayGamma = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTEtaDecayGamma")->Clone("hIsoGammaM02pTEtaDecayGamma");
    h2M02vspTMergedPionGamma = (TH2F*)dIsoGammaQA->Get("hIsoGammaM02pTMergedPionGamma")->Clone("hIsoGammaM02pTMergedPionGamma");

    int NSlicesPtBins=SlicesPt.size();

    //Initiate arrays for storing Plotting1D's for loop
    Plotting1D PM02vspTSignalBackground[NSlicesPtBins];
    Plotting1D PM02vspTSignalContributions[NSlicesPtBins];

    TH1F* hM02vspTSlice[NSlicesPtBins];
    TH1F* hM02vspTSliceSignal[NSlicesPtBins];
    TH1F* hM02vspTSliceBackground[NSlicesPtBins];
    TH1F* hM02vspTSliceAllGamma[NSlicesPtBins];
    TH1F* hM02vspTSlicePionDecayGamma[NSlicesPtBins];
    TH1F* hM02vspTSliceEtaDecayGamma[NSlicesPtBins];
    TH1F* hM02vspTSliceMergedPionGamma[NSlicesPtBins];

    PlottingGrid PM02Grid;
    PM02Grid.SetAxisLabel("#bf{#it{M}_{02}}", "#bf{#it{N}}");

    for(int i=0;i<NSlicesPtBins-1;i++){
      //Plotting1D PM02vspTSignalBackground;
      //Plotting1D PM02vspTSignalContributions;
      //For saving slices
      //TH1F* hM02vspTSignalBackground=new TH1F(Form("Bkg%i",i));// = 0x0;//Dummy for plotting in the end.
      //TH1F* hM02vspTSliceSignalContributions=new TH1F(Form("Cont%i",i));;// = 0x0;
      
      //cout<<hM02vspTSignalBackground<<

      //Project TH2F's
      string titleM02vspTsplice = "M02:"+to_string(SlicesPt.at(i)).substr(0,4)+"GeV/c<p_{T}<"+to_string(SlicesPt.at(i+1)).substr(0,4)+"GeV/c";
      int pTminbin=h2M02vspT->GetYaxis()->FindBin(SlicesPt.at(i));
      int pTmaxbin=h2M02vspT->GetYaxis()->FindBin(SlicesPt.at(i+1));
      
      //hM02vspTSlices.push_back(h2M02vspT->ProjectionX(titleM02vspTsplice.c_str(),pTminbin,pTmaxbin));
      hM02vspTSlice[i] = (TH1F*)h2M02vspT->ProjectionX(titleM02vspTsplice.c_str(),pTminbin,pTmaxbin);
      hM02vspTSlice[i]->SetTitle(titleM02vspTsplice.c_str());
      hM02vspTSliceSignal[i] = (TH1F*)h2M02vspTSignal->ProjectionX(Form("a%i",i),pTminbin,pTmaxbin);
      hM02vspTSliceBackground[i] = (TH1F*)h2M02vspTBackground->ProjectionX(Form("b%i",i),pTminbin,pTmaxbin);
      hM02vspTSliceAllGamma[i] = (TH1F*)h2M02vspTAllGamma->ProjectionX(Form("c%i",i),pTminbin,pTmaxbin);
      hM02vspTSlicePionDecayGamma[i] = (TH1F*)h2M02vspTPionDecayGamma->ProjectionX(Form("d%i",i),pTminbin,pTmaxbin);
      hM02vspTSliceEtaDecayGamma[i] = (TH1F*)h2M02vspTEtaDecayGamma->ProjectionX(Form("e%i",i),pTminbin,pTmaxbin);
      hM02vspTSliceMergedPionGamma[i] = (TH1F*)h2M02vspTMergedPionGamma->ProjectionX(Form("f%i",i),pTminbin,pTmaxbin);
//Unique hist name.
//Try array.
      AvoidLogYconflict(hM02vspTSlice[i]);
      AvoidLogYconflict(hM02vspTSliceSignal[i]);
      AvoidLogYconflict(hM02vspTSliceBackground[i]);
      AvoidLogYconflict(hM02vspTSliceAllGamma[i]);
      AvoidLogYconflict(hM02vspTSlicePionDecayGamma[i]);
      AvoidLogYconflict(hM02vspTSliceMergedPionGamma[i]);

      //Add to figure and plot:
      //All, Signal and Background
      PM02vspTSignalBackground[i].New(hM02vspTSlice[i],"All");
      PM02vspTSignalBackground[i].New(hM02vspTSliceSignal[i],"Signal");
      PM02vspTSignalBackground[i].New(hM02vspTSliceBackground[i],"Background");
      PM02vspTSignalBackground[i].SetAxisLabel(("M02"));
      PM02vspTSignalBackground[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalBackgroundSig.%s", outputDir.Data(),SlicesPt.at(i),SlicesPt.at(i+1), suffix),kFALSE,kTRUE);
//
      ////All, Signal and Contributions
      PM02vspTSignalContributions[i].New(hM02vspTSlice[i],"All");
      PM02vspTSignalContributions[i].New(hM02vspTSliceSignal[i],"Signal");
      PM02vspTSignalContributions[i].New(hM02vspTSlicePionDecayGamma[i],"#pi^{0}-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceEtaDecayGamma[i],"#eta-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceMergedPionGamma[i],"#pi^{0}");
      PM02vspTSignalContributions[i].SetAxisLabel(("M02"));
      PM02vspTSignalContributions[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalContributions.%s", outputDir.Data(),SlicesPt.at(i),SlicesPt.at(i+1), suffix),kFALSE,kTRUE);

      PM02Grid.New(hM02vspTSlice[i],"All");
      PM02Grid.New(hM02vspTSliceSignal[i],"Signal");
      PM02Grid.New(hM02vspTSliceBackground[i],"Background");
      PM02Grid.NextPad(Form("%.1f < #it{p}_{T} (GeV/c) < %.1f", SlicesPt.at(i),SlicesPt.at(i+1)));
    }
    PM02Grid.Plot(Form("%s/M02_SignalBackgroundSig_Grid.%s", outputDir.Data(),suffix),kFALSE,kFALSE);
    


    //
  } else {
    hMinMassDiffToPi0 = (TH1F*)dIsoGammaQA->Get("hMinMassDiffToPi0");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0);
  }
  PMinMassDiffToPi0.SetAxisLabel("#bf{#Delta#it{m}_{inv} (GeV/#it{c}^{2})}");
  PMinMassDiffToPi0.Plot(Form("%s/MinMassDiffToPi0.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaE.SetAxisLabel("#bf{#it{E}[GeV]}");
  PIsoGammaE.Plot(Form("%s/IsoGammaE.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaM02.SetAxisLabel("#bf{#it{M02}}");
  PIsoGammaM02.Plot(Form("%s/IsoGammaM02.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaM20.SetAxisLabel("#bf{#it{M20}}");
  PIsoGammaM20.Plot(Form("%s/IsoGammaM20.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaPx.SetAxisLabel("#bf{#it{Px}[GeV/#it{c}]}");
  PIsoGammaPx.Plot(Form("%s/IsoGammaPx.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaPy.SetAxisLabel("#bf{#it{Py}[GeV/#it{c}]}");
  PIsoGammaPy.Plot(Form("%s/IsoGammaPy.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaIsoChargedCorrected.SetAxisLabel("#bf{#it{p^{iso}_{T}}[GeV/#it{c}]}");
  PIsoGammaIsoChargedCorrected.Plot(Form("%s/IsoGammaIsoChargedCorrected.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  PIsoGammaPz.SetAxisLabel("#bf{#it{Pz}[GeV/#it{c}]}");
  PIsoGammaPz.Plot(Form("%s/IsoGammaPz.%s", outputDir.Data(), suffix));

  PIsoGammaIsoCharged.SetAxisLabel("#bf{#it{p^{iso}_{T}}[GeV/#it{c}]}");
  PIsoGammaIsoCharged.Plot(Form("%s/IsoGammaIsoCharged.%s", outputDir.Data(), suffix),kFALSE,kTRUE);

  

  EXIT
}

void plotHistosFromTree(TString AnalysisDirectory, bool isDebugRun = false)
{
  GlobalOptions optns(AnalysisDirectory, isDebugRun);

  TString inputFilePath = Form("%s/HistosFromTree%s.root", AnalysisDirectory.Data(), isDebugRun ? "_0" : "");

  if (!std::filesystem::exists(inputFilePath.Data())) {
    inputFilePath = Form("%s/HistosFromTree_0.root", AnalysisDirectory.Data());
    INFO(Form("Did not find HistosFromTree.root, looking for %s now.", inputFilePath.Data()))
  }

  TFile* fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  TDirectory* dJets = (TDirectory*)fIn->Get("Jets");
  if (!dJets)
    FATAL("Dir Jets not found")

  plotJets(dJets, optns);

  TDirectory* dJetQA = (TDirectory*)fIn->Get("JetQA");
  if (!dJetQA)
    FATAL("Dir JetQA not found")

  plotJetQA(dJetQA, optns);

  TDirectory* dIsoGammaQA = (TDirectory*)fIn->Get("IsoGammaQA");
  if (!dIsoGammaQA)
    FATAL("Dir IsoGammaQA not found")

  plotIsoGammaQA(dIsoGammaQA, optns);

  fIn->Close();
}