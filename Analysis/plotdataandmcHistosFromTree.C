#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
//#include "plotHistosfromTree.h"

const char *suffix = "png";
std::vector<float> SlicesPt={3,4,5,6,7,8,9,10,11,12,14,16,18,20};

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

void plotdataandmcHistosFromTree(string outputDirectory="PlotMCandDataTest/100_Charged/Standard",string filepathData="LHC18_pass2/100_Charged/Standard/",string filepathMC="MergedGJ18DandJJ18D/100_Charged/Standard/"){
    //Load data:
    TFile *fileData = TFile::Open(Form("%sHistosFromTree.root",filepathData.c_str()),"READ");
    TFile *fileMC = TFile::Open(Form("%sHistosFromTree.root",filepathMC.c_str()),"READ");

    //Get directories for IsoGammaQA:
    TDirectory *dIsoGammaQAData = (TDirectory *)fileData->Get("IsoGammaQA");
    TDirectory *dIsoGammaQAMC = (TDirectory *)fileMC->Get("IsoGammaQA");

    //Trigger string (to be put in plot):
    TString TriggerString;
    if(outputDirectory.find("100_Charged")!=string::npos){
      TriggerString="Trigger: kINT7";
    }else if(outputDirectory.find("101_Charged")!=string::npos){
      TriggerString="Trigger: EG2";
    }else if(outputDirectory.find("102_Charged")!=string::npos){
      TriggerString="Trigger: EG1";
    }else{
      WARN("No Trigger registered in folder name.")
      TriggerString="?";
    }


    //Create directory for outputfile:
    TString outputDir = Form(outputDirectory.c_str());

    createDirectory(outputDirectory.c_str());

    // Plot response matrix
    Plotting1D PMinMassDiffToPi0;
    TH1F *hMinMassDiffToPi0Signal;
    TH1F *hMinMassDiffToPi0Background;
    TH1F *hMinMassDiffToPi0Data;
    TH1F *hMinMassDiffToPi0;

    Plotting1D PIsoGammaE;
    TH1F *hIsoGammaESignal;
    TH1F *hIsoGammaEBackground;
    TH1F *hIsoGammaEData;
    TH1F *hIsoGammaE;

    Plotting1D PIsoGammaM02;
    TH1F *hIsoGammaM02Signal;
    TH1F *hIsoGammaM02Background;
    TH1F *hIsoGammaM02Data;
    TH1F *hIsoGammaM02;

    Plotting1D PIsoGammaM20;
    TH1F *hIsoGammaM20Signal;
    TH1F *hIsoGammaM20Background;
    TH1F *hIsoGammaM20Data;
    TH1F *hIsoGammaM20;

    Plotting1D PIsoGammaPx;
    TH1F *hIsoGammaPxSignal;
    TH1F *hIsoGammaPxBackground;
    TH1F *hIsoGammaPxData;
    TH1F *hIsoGammaPx;

    Plotting1D PIsoGammaPy;
    TH1F *hIsoGammaPySignal;
    TH1F *hIsoGammaPyBackground;
    TH1F *hIsoGammaPyData;
    TH1F *hIsoGammaPy;

    Plotting1D PIsoGammaPz;
    TH1F *hIsoGammaPzSignal;
    TH1F *hIsoGammaPzBackground;
    TH1F *hIsoGammaPzData;
    TH1F *hIsoGammaPz;

    Plotting1D PIsoGammaIsoCharged;
    TH1F *hIsoGammaIsoChargedSignal;
    TH1F *hIsoGammaIsoChargedBackground;
    TH1F *hIsoGammaIsoChargedData;
    TH1F *hIsoGammaIsoCharged;

    Plotting1D PIsoGammaIsoChargedCorrected;
    TH1F *hIsoGammaIsoChargedCorrectedSignal;
    TH1F *hIsoGammaIsoChargedCorrectedBackground;
    TH1F *hIsoGammaIsoChargedCorrectedData;
    TH1F *hIsoGammaIsoChargedCorrected;


    // Diff to pi0 mass
    hMinMassDiffToPi0Signal = (TH1F *)dIsoGammaQAMC->Get("hMinMassDiffToPi0Signal")->Clone("hMinMassDiffToPi0Signal");
    hMinMassDiffToPi0Background = (TH1F *)dIsoGammaQAMC->Get("hMinMassDiffToPi0Background")->Clone("hMinMassDiffToPi0Background");
    hMinMassDiffToPi0Data = (TH1F *)dIsoGammaQAData->Get("hMinMassDiffToPi0")->Clone("hMinMassDiffToPi0");
    float integralSignal = hMinMassDiffToPi0Signal->Integral(1, hMinMassDiffToPi0Signal->GetNbinsX());
    float integralBG = hMinMassDiffToPi0Background->Integral(1, hMinMassDiffToPi0Background->GetNbinsX());
    float integralData = hMinMassDiffToPi0Data->Integral(1, hMinMassDiffToPi0Data->GetNbinsX());
    hMinMassDiffToPi0Signal->Scale(1. / integralSignal);
    hMinMassDiffToPi0Background->Scale(1. / integralBG);
    hMinMassDiffToPi0Data->Scale(1. / integralData);
    AvoidLogYconflict(hMinMassDiffToPi0Signal);
    AvoidLogYconflict(hMinMassDiffToPi0Background);
    AvoidLogYconflict(hMinMassDiffToPi0Data);
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Signal, "Signal");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Background, "Background");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0Data, "Data");
    PMinMassDiffToPi0.SetLegend(0.4,0.75);
    PMinMassDiffToPi0.NewLatex(0.75,0.8,TriggerString);


    // E
    hIsoGammaESignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaESignal")->Clone("hIsoGammaESignal");
    hIsoGammaEBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaEBackground")->Clone("hIsoGammaEBackground");
    hIsoGammaEData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaE")->Clone("hIsoGammaE");
    float integralESignal = hIsoGammaESignal->Integral(1, hIsoGammaESignal->GetNbinsX());
    float integralEBG = hIsoGammaEBackground->Integral(1, hIsoGammaEBackground->GetNbinsX());
    float integralEData = hIsoGammaEData->Integral(1,hIsoGammaEData->GetNbinsX());
    hIsoGammaESignal->Scale(1. / integralESignal);
    hIsoGammaEBackground->Scale(1. / integralEBG);
    hIsoGammaEData->Scale(1. / integralData);
    AvoidLogYconflict(hIsoGammaESignal);
    AvoidLogYconflict(hIsoGammaEBackground);
    AvoidLogYconflict(hIsoGammaEData);
    PIsoGammaE.New(hIsoGammaESignal, "Signal");
    PIsoGammaE.New(hIsoGammaEBackground, "Background");
    PIsoGammaE.New(hIsoGammaEData,"Data");
    PIsoGammaE.NewLatex(0.75,0.8,TriggerString);

    // M02
    hIsoGammaM02Signal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaM02Signal")->Clone("hIsoGammaM02Signal");
    hIsoGammaM02Background = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaM02Background")->Clone("hIsoGammaM02Background");
    hIsoGammaM02Data = (TH1F *)dIsoGammaQAData->Get("hIsoGammaM02")->Clone("hIsoGammaM02");
    float integralM02Signal = hIsoGammaM02Signal->Integral(1, hIsoGammaM02Signal->GetNbinsX());
    float integralM02BG = hIsoGammaM02Background->Integral(1, hIsoGammaM02Background->GetNbinsX());
    float integralM02Data = hIsoGammaM02Data->Integral(1, hIsoGammaM02Data->GetNbinsX());
    hIsoGammaM02Signal->Scale(1. / integralM02Signal);
    hIsoGammaM02Background->Scale(1. / integralM02BG);
    hIsoGammaM02Data->Scale(1. / integralM02Data);
    AvoidLogYconflict(hIsoGammaM02Signal);
    AvoidLogYconflict(hIsoGammaM02Background);
    AvoidLogYconflict(hIsoGammaM02Data);
    PIsoGammaM02.New(hIsoGammaM02Signal, "Signal");
    PIsoGammaM02.New(hIsoGammaM02Background, "Background");
    PIsoGammaM02.New(hIsoGammaM02Data, "Data");
    PIsoGammaM02.NewLatex(0.75,0.8,TriggerString);

    // M20
    hIsoGammaM20Signal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaM20Signal")->Clone("hIsoGammaM20Signal");
    hIsoGammaM20Background = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaM20Background")->Clone("hIsoGammaM20Background");
    hIsoGammaM20Data = (TH1F *)dIsoGammaQAData->Get("hIsoGammaM20")->Clone("hIsoGammaM20");
    float integralM20Signal = hIsoGammaM20Signal->Integral(1, hIsoGammaM20Signal->GetNbinsX());
    float integralM20BG = hIsoGammaM20Background->Integral(1, hIsoGammaM20Background->GetNbinsX());
    float integralM20Data = hIsoGammaM20Data->Integral(1, hIsoGammaM20Data->GetNbinsX());
    hIsoGammaM20Signal->Scale(1. / integralM20Signal);
    hIsoGammaM20Background->Scale(1. / integralM20BG);
    hIsoGammaM20Data->Scale(1. / integralM20Data);
    AvoidLogYconflict(hIsoGammaM20Signal);
    AvoidLogYconflict(hIsoGammaM20Background);
    AvoidLogYconflict(hIsoGammaM20Data);
    PIsoGammaM20.New(hIsoGammaM20Signal, "Signal");
    PIsoGammaM20.New(hIsoGammaM20Background, "Background");
    PIsoGammaM20.New(hIsoGammaM20Data, "Data");
    PIsoGammaM20.NewLatex(0.75,0.8,TriggerString);

    // Px
    hIsoGammaPxSignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPxSignal")->Clone("hIsoGammaPxSignal");
    hIsoGammaPxBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPxBackground")->Clone("hIsoGammaPxBackground");
    hIsoGammaPxData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaPx")->Clone("hIsoGammaPx");
    float integralPxSignal = hIsoGammaPxSignal->Integral(1, hIsoGammaPxSignal->GetNbinsX());
    float integralPxBG = hIsoGammaPxBackground->Integral(1, hIsoGammaPxBackground->GetNbinsX());
    float integralPxData = hIsoGammaPxData->Integral(1, hIsoGammaPxData->GetNbinsX());
    hIsoGammaPxSignal->Scale(1. / integralPxSignal);
    hIsoGammaPxBackground->Scale(1. / integralPxBG);
    hIsoGammaPxData->Scale(1. / integralPxData);
    AvoidLogYconflict(hIsoGammaPxSignal);
    AvoidLogYconflict(hIsoGammaPxBackground);
    AvoidLogYconflict(hIsoGammaPxData);
    PIsoGammaPx.New(hIsoGammaPxSignal, "Signal");
    PIsoGammaPx.New(hIsoGammaPxBackground, "Background");
    PIsoGammaPx.New(hIsoGammaPxData, "Data");
    PIsoGammaPx.NewLatex(0.75,0.8,TriggerString);

    // Py
    hIsoGammaPySignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPySignal")->Clone("hIsoGammaPySignal");
    hIsoGammaPyBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPyBackground")->Clone("hIsoGammaPyBackground");
    hIsoGammaPyData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaPy")->Clone("hIsoGammaPy");
    float integralPySignal = hIsoGammaPySignal->Integral(1, hIsoGammaPySignal->GetNbinsX());
    float integralPyBG = hIsoGammaPyBackground->Integral(1, hIsoGammaPyBackground->GetNbinsX());
    float integralPyData = hIsoGammaPyData->Integral(1, hIsoGammaPyData->GetNbinsX());
    hIsoGammaPySignal->Scale(1. / integralPySignal);
    hIsoGammaPyBackground->Scale(1. / integralPyBG);
    hIsoGammaPyData->Scale(1. / integralPyData);
    AvoidLogYconflict(hIsoGammaPySignal);
    AvoidLogYconflict(hIsoGammaPyBackground);
    AvoidLogYconflict(hIsoGammaPyData);
    PIsoGammaPy.New(hIsoGammaPySignal, "Signal");
    PIsoGammaPy.New(hIsoGammaPyBackground, "Background");
    PIsoGammaPy.New(hIsoGammaPyData, "Data");
    PIsoGammaPy.NewLatex(0.75,0.8,TriggerString);

    // Pz
    hIsoGammaPzSignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPzSignal")->Clone("hIsoGammaPzSignal");
    hIsoGammaPzBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaPzBackground")->Clone("hIsoGammaPzBackground");
    hIsoGammaPzData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaPz")->Clone("hIsoGammaPz");
    float integralPzSignal = hIsoGammaPzSignal->Integral(1, hIsoGammaPzSignal->GetNbinsX());
    float integralPzBG = hIsoGammaPzBackground->Integral(1, hIsoGammaPzBackground->GetNbinsX());
    float integralPzData = hIsoGammaPzData->Integral(1, hIsoGammaPzData->GetNbinsX());
    hIsoGammaPzSignal->Scale(1. / integralPzSignal);
    hIsoGammaPzBackground->Scale(1. / integralPzBG);
    hIsoGammaPzData->Scale(1. / integralPzData);
    AvoidLogYconflict(hIsoGammaPzSignal);
    AvoidLogYconflict(hIsoGammaPzBackground);
    AvoidLogYconflict(hIsoGammaPzData);
    PIsoGammaPz.New(hIsoGammaPzSignal, "Signal");
    PIsoGammaPz.New(hIsoGammaPzBackground, "Background");
    PIsoGammaPz.New(hIsoGammaPzData, "Data");
    PIsoGammaPz.NewLatex(0.75,0.8,TriggerString);

    // IsoCharged
    hIsoGammaIsoChargedSignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaIsoChargedSignal")->Clone("hIsoGammaIsoChargedSignal");
    hIsoGammaIsoChargedBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaIsoChargedBackground")->Clone("hIsoGammaIsoChargedBackground");
    hIsoGammaIsoChargedData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaIsoCharged")->Clone("hIsoGammaIsoCharged");
    float integralIsolationSignal = hIsoGammaIsoChargedSignal->Integral(1, hIsoGammaIsoChargedSignal->GetNbinsX());
    float integralIsolationBG = hIsoGammaIsoChargedBackground->Integral(1, hIsoGammaIsoChargedBackground->GetNbinsX());
    float integralIsolationData = hIsoGammaIsoChargedData->Integral(1, hIsoGammaIsoChargedData->GetNbinsX());
    hIsoGammaIsoChargedSignal->Scale(1. / integralIsolationSignal);
    hIsoGammaIsoChargedBackground->Scale(1. / integralIsolationBG);
    hIsoGammaIsoChargedData->Scale(1. / integralIsolationData);
    AvoidLogYconflict(hIsoGammaIsoChargedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedBackground);
    AvoidLogYconflict(hIsoGammaIsoChargedData);
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedSignal, "Signal");
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedBackground, "Background");
    PIsoGammaIsoCharged.New(hIsoGammaIsoChargedData, "Data");
    PIsoGammaIsoCharged.NewLatex(0.75,0.8,TriggerString);

    // IsoChargedCorrected
    hIsoGammaIsoChargedCorrectedSignal = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaIsoChargedCorrectedSignal")->Clone("hIsoGammaIsoChargedCorrectedSignal");
    hIsoGammaIsoChargedCorrectedBackground = (TH1F *)dIsoGammaQAMC->Get("hIsoGammaIsoChargedCorrectedBackground")->Clone("hIsoGammaIsoChargedCorrectedBackground");
    hIsoGammaIsoChargedCorrectedData = (TH1F *)dIsoGammaQAData->Get("hIsoGammaIsoChargedCorrected")->Clone("hIsoGammaIsoChargedCorrected");
    float integralIsolationCorrectedSignal = hIsoGammaIsoChargedCorrectedSignal->Integral(1, hIsoGammaIsoChargedCorrectedSignal->GetNbinsX());
    float integralIsolationCorrectedBG = hIsoGammaIsoChargedCorrectedBackground->Integral(1, hIsoGammaIsoChargedCorrectedBackground->GetNbinsX());
    float integralIsolationCorrectedData = hIsoGammaIsoChargedCorrectedData->Integral(1, hIsoGammaIsoChargedCorrectedData->GetNbinsX());
    hIsoGammaIsoChargedCorrectedSignal->Scale(1. / integralIsolationSignal);
    hIsoGammaIsoChargedCorrectedBackground->Scale(1. / integralIsolationBG);
    hIsoGammaIsoChargedCorrectedData->Scale(1. / integralIsolationData);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedSignal);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedBackground);
    AvoidLogYconflict(hIsoGammaIsoChargedCorrectedData);
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedSignal, "Signal");
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedBackground, "Background");
    PIsoGammaIsoChargedCorrected.New(hIsoGammaIsoChargedCorrectedData, "Data");
    PIsoGammaIsoChargedCorrected.NewLatex(0.75,0.8,TriggerString);

    // M02vsPt (Slicing 2D hist pr. specified pT-bins)
    TH2F *h2M02vspT;
    TH2F *h2M02vspTData;
    TH2F *h2M02vspTSignal;
    TH2F *h2M02vspTBackground;
    TH2F *h2M02vspTAllGamma;
    TH2F *h2M02vspTPionDecayGamma;
    TH2F *h2M02vspTEtaDecayGamma;
    TH2F *h2M02vspTMergedPionGamma;

    h2M02vspT = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pT")->Clone("hIsoGammaM02pT");
    h2M02vspTData = (TH2F *)dIsoGammaQAData->Get("hIsoGammaM02pT")->Clone("hIsoGammaM02pT");
    h2M02vspTSignal = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTSignal")->Clone("hIsoGammaM02pTSignal");
    h2M02vspTBackground = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTBackground")->Clone("hIsoGammaM02pTBackground");
    h2M02vspTAllGamma = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTAllGamma")->Clone("hIsoGammaM02pTAllGamma");
    h2M02vspTPionDecayGamma = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTPionDecayGamma")->Clone("hIsoGammaM02pTPionDecayGamma");
    h2M02vspTEtaDecayGamma = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTEtaDecayGamma")->Clone("hIsoGammaM02pTEtaDecayGamma");
    h2M02vspTMergedPionGamma = (TH2F *)dIsoGammaQAMC->Get("hIsoGammaM02pTMergedPionGamma")->Clone("hIsoGammaM02pTMergedPionGamma");

    int NSlicesPtBins = SlicesPt.size();

    // Initiate arrays for storing Plotting1D's for loop
    Plotting1D PM02vspTSignalBackground[NSlicesPtBins];
    Plotting1D PM02vspTSignalContributions[NSlicesPtBins];

    TH1F *hM02vspTSlice[NSlicesPtBins];
    TH1F *hM02vspTSliceData[NSlicesPtBins];
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

      // hM02vspTSlices.push_back(h2M02vspT->ProjectionX(titleM02vspTsplice.c_str(),pTminbin,pTmaxbin));
      hM02vspTSlice[i] = (TH1F *)h2M02vspT->ProjectionX(titleM02vspTsplice.c_str(), pTminbin, pTmaxbin);
      hM02vspTSlice[i]->SetTitle(titleM02vspTsplice.c_str());
      hM02vspTSliceData[i] = (TH1F *)h2M02vspTData->ProjectionX(Form("Data%i",i), pTminbin, pTmaxbin);
      hM02vspTSliceSignal[i] = (TH1F *)h2M02vspTSignal->ProjectionX(Form("a%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceBackground[i] = (TH1F *)h2M02vspTBackground->ProjectionX(Form("b%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceAllGamma[i] = (TH1F *)h2M02vspTAllGamma->ProjectionX(Form("c%i", i), pTminbin, pTmaxbin);
      hM02vspTSlicePionDecayGamma[i] = (TH1F *)h2M02vspTPionDecayGamma->ProjectionX(Form("d%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceEtaDecayGamma[i] = (TH1F *)h2M02vspTEtaDecayGamma->ProjectionX(Form("e%i", i), pTminbin, pTmaxbin);
      hM02vspTSliceMergedPionGamma[i] = (TH1F *)h2M02vspTMergedPionGamma->ProjectionX(Form("f%i", i), pTminbin, pTmaxbin);
      // Unique hist name.
      // Try array.
      AvoidLogYconflict(hM02vspTSlice[i]);
      AvoidLogYconflict(hM02vspTSliceData[i]);
      AvoidLogYconflict(hM02vspTSliceSignal[i]);
      AvoidLogYconflict(hM02vspTSliceBackground[i]);
      AvoidLogYconflict(hM02vspTSliceAllGamma[i]);
      AvoidLogYconflict(hM02vspTSlicePionDecayGamma[i]);
      AvoidLogYconflict(hM02vspTSliceMergedPionGamma[i]);

      // Add to figure and plot:
      // All, Signal and Background
      PM02vspTSignalBackground[i].New(hM02vspTSlice[i], "All");
      PM02vspTSignalBackground[i].New(hM02vspTSliceData[i], "Data");
      PM02vspTSignalBackground[i].New(hM02vspTSliceSignal[i], "Signal");
      PM02vspTSignalBackground[i].New(hM02vspTSliceBackground[i], "Background");
      PM02vspTSignalBackground[i].SetAxisLabel(("M02"));
      PM02vspTSignalBackground[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalBackgroundSig.%s", outputDir.Data(), SlicesPt.at(i), SlicesPt.at(i + 1), suffix), kFALSE, kTRUE);
      //
      ////All, Signal and Contributions
      PM02vspTSignalContributions[i].New(hM02vspTSlice[i], "All");
      PM02vspTSignalContributions[i].New(hM02vspTSliceData[i], "Data");
      PM02vspTSignalContributions[i].New(hM02vspTSliceSignal[i], "Signal");
      PM02vspTSignalContributions[i].New(hM02vspTSlicePionDecayGamma[i], "#pi^{0}-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceEtaDecayGamma[i], "#eta-decay");
      PM02vspTSignalContributions[i].New(hM02vspTSliceMergedPionGamma[i], "#pi^{0}");
      PM02vspTSignalContributions[i].SetAxisLabel(("M02"));
      PM02vspTSignalContributions[i].NewLatex(0.75,0.8,TriggerString);
      PM02vspTSignalContributions[i].Plot(Form("%s/M02_pTmin%.1f_pTmax%.1f_SignalContributions.%s", outputDir.Data(), SlicesPt.at(i), SlicesPt.at(i + 1), suffix), kFALSE, kTRUE);

      PM02Grid.New(hM02vspTSlice[i], "All");
      PM02Grid.New(hM02vspTSliceData[i], "Data");
      PM02Grid.New(hM02vspTSliceSignal[i], "Signal");
      PM02Grid.New(hM02vspTSliceBackground[i], "Background");
      PM02Grid.NewLatex(0.7,0.8,TriggerString);
      PM02Grid.NextPad(Form("%.1f < #it{p}_{T} (GeV/c) < %.1f", SlicesPt.at(i), SlicesPt.at(i + 1)));
    }
    PM02Grid.Plot(Form("%s/M02_SignalBackgroundSig_Grid.%s", outputDir.Data(), suffix), kFALSE, kFALSE);


    hMinMassDiffToPi0 = (TH1F *)dIsoGammaQAData->Get("hMinMassDiffToPi0");
    PMinMassDiffToPi0.New(hMinMassDiffToPi0);

    PMinMassDiffToPi0.SetAxisLabel("#bf{#Delta#it{m}_{inv} (GeV/#it{c}^{2})}");
    PMinMassDiffToPi0.Plot(Form("%s/MinMassDiffToPi0.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

    PIsoGammaE.SetAxisLabel("#bf{#it{E}[GeV]}");
    PIsoGammaE.Plot(Form("%s/IsoGammaE.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

    PIsoGammaM02.SetAxisLabel("#bf{#it{M02}}");
    PIsoGammaM02.Plot(Form("%s/IsoGammaM02.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

    PIsoGammaM20.SetAxisLabel("#bf{#it{M20}}");
    PIsoGammaM20.Plot(Form("%s/IsoGammaM20.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

    PIsoGammaPx.SetAxisLabel("#bf{#it{Px}[GeV/#it{c}]}");
    PIsoGammaPx.Plot(Form("%s/IsoGammaPx.%s", outputDir.Data(), suffix));

    PIsoGammaPy.SetAxisLabel("#bf{#it{Py}[GeV/#it{c}]}");
    PIsoGammaPy.Plot(Form("%s/IsoGammaPy.%s", outputDir.Data(), suffix));

    PIsoGammaIsoChargedCorrected.SetAxisLabel("#bf{#it{p^{iso}_{T}}[GeV/#it{c}]}");
    PIsoGammaIsoChargedCorrected.Plot(Form("%s/IsoGammaIsoChargedCorrected.%s", outputDir.Data(), suffix), kFALSE, kTRUE);

    PIsoGammaPz.SetAxisLabel("#bf{#it{Pz}[GeV/#it{c}]}");
    PIsoGammaPz.Plot(Form("%s/IsoGammaPz.%s", outputDir.Data(), suffix));

    PIsoGammaIsoCharged.SetAxisLabel("#bf{#it{p^{iso}_{T}}[GeV/#it{c}]}");
    PIsoGammaIsoCharged.Plot(Form("%s/IsoGammaIsoCharged.%s", outputDir.Data(), suffix), kFALSE, kTRUE);


}