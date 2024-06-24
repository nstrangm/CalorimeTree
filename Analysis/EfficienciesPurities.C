#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <filesystem>
#include "Utilities.h"
#include "Cuts.h"

const char* suffix = "png";

int nPtBinEdges=5;
float pTbinEdges[5] = {10,15,20,40,80};
float pTbinsToplot[4] = {12.5,17.5,30,60};

/*
This macro calculates the acceptance, efficiency and purity of prompt photons in MC.
It will also calculate the data-driven purity with the ABCD-method.

It takes as its arguments two paths to folders with root files produced by 
"makeHistosFromTree.C" and a path to an output directory. 

Mandatory arguments:
-One input file has to be from MC and one from actual data.
-The output directory is created if it does not already exist.
Optional arguments:
-Optional flags are supplied in the macro call, to allow for calculations to be skipped.

NOTES:
This script assumes that the same binning is used for all histograms of the same type.
I.e. histogram binning in MC and Data must be the same.
*/

void EfficienciesPurities(string AnalysisDirectoryMC,string AnalysisDirectoryData, const char* outputDir, bool doABCD = true, bool doPurity = true, bool doAcceptance = true, bool doEfficiency = true){
  //GlobalOptions optns(AnalysisDirectory, jobId);

  TString inputFilePathMC = Form("%s/HistosFromTree.root", AnalysisDirectoryMC.c_str());
  TString inputFilePathData = Form("%s/HistosFromTree.root", AnalysisDirectoryData.c_str());

  TString TriggerString;
  if(AnalysisDirectoryMC.find("100_Charged")!=string::npos){
    TriggerString="Trigger: kINT7";
  }else if(AnalysisDirectoryMC.find("101_Charged")!=string::npos){
    TriggerString="Trigger: EG2";
  }else if(AnalysisDirectoryMC.find("102_Charged")!=string::npos){
    TriggerString="Trigger: EG1";
  }else{
    WARN("No Trigger registered in folder name.")
    TriggerString="?";
  }

  TFile* fInMC = new TFile(inputFilePathMC, "READ");
  TFile* fInData = new TFile(inputFilePathData, "READ");
  TDirectory* dIsoGammaMC = (TDirectory*)fInMC->Get("IsoGammas");
  TDirectory* dIsoGammaData = (TDirectory*)fInData->Get("IsoGammas");

  if(doPurity || doEfficiency || doAcceptance || doABCD){
    createDirectory(outputDir);
    //Purity and Efficiency from MC:
    if(doPurity || doEfficiency || doAcceptance){

      TH1F* hIsoGammaPt;
      TH1F* hGammaGenPt;
      hIsoGammaPt = (TH1F*)dIsoGammaMC->Get("hIsoGammaPt")->Clone("hIsoGammaPt");
      hIsoGammaPt->Rebin(4);
      hGammaGenPt = (TH1F*)dIsoGammaMC->Get("hGammaGenPt")->Clone("hGammaGenPt");
      hGammaGenPt->Rebin(4 );
      

      if(doPurity){
        Plotting1D PIsoGammaPurity;
        TH1F* hIsoGammaPurity;
        hIsoGammaPurity = (TH1F*)dIsoGammaMC->Get("hIsoGammaPtSignal")->Clone("hIsoGammaPtSignal");
        hIsoGammaPurity->Rebin(4);
        hIsoGammaPurity->Divide(hIsoGammaPt);
        hIsoGammaPurity->SetTitle("Purity");
        PIsoGammaPurity.New(hIsoGammaPurity,TriggerString);
        PIsoGammaPurity.SetLegend(0.4,0.75,0.2,0.4);
        PIsoGammaPurity.SetAxisLabel("#bf{#it{p^{iso}_{T}}[GeV/#it{c}]}");
        PIsoGammaPurity.Plot(Form("%s/Purity.%s", outputDir, suffix));
      }
      if(doEfficiency){
        Plotting1D PIsoGammaEfficiency;
        TH1F* hIsoGammaEfficiency;
        hIsoGammaEfficiency = (TH1F*)dIsoGammaMC->Get("hIsoGammaPtSignal")->Clone("hIsoGammaPtSignal");
        hIsoGammaEfficiency->Rebin(4);
        hIsoGammaEfficiency->Divide(hGammaGenPt);
        hIsoGammaEfficiency->SetTitle("Efficiency");
        PIsoGammaEfficiency.New(hIsoGammaEfficiency,TriggerString);
        PIsoGammaEfficiency.SetLegend(0.2,0.4,0.6,0.7);
        PIsoGammaEfficiency.SetAxisLabel("#bf{#it{p_{T}}[GeV/#it{c}]}");
        PIsoGammaEfficiency.Plot(Form("%s/Efficiency.%s", outputDir, suffix));
      }
      if(doAcceptance){
        Plotting1D PIsoGammaAcceptance;
        TH1F* hIsoGammaAcceptance;
        hIsoGammaAcceptance = (TH1F*)dIsoGammaMC->Get("hIsoGammaPtSignal")->Clone("hIsoGammaPtSignal");
        hIsoGammaAcceptance->Rebin(4);
        hIsoGammaAcceptance->Divide(hGammaGenPt);
        hIsoGammaAcceptance->SetTitle("Acceptance");
        PIsoGammaAcceptance.New(hIsoGammaAcceptance,TriggerString);
        PIsoGammaAcceptance.SetLegend(0.4,0.75,0.2,0.4);
        PIsoGammaAcceptance.SetAxisLabel("#bf{#it{p_{T}}[GeV/#it{c}]}");
        PIsoGammaAcceptance.Plot(Form("%s/Acceptance.%s", outputDir, suffix));
      }
    }
    if(doABCD){
      //Load THnSparse:
      //Load MC data:
      THnSparse *PtIsovsM02vsPtMC = (THnSparse*)dIsoGammaMC->Get("hIsoGammaIsovsM02vsPtnoCuts");
      THnSparse *PtIsovsM02vsPtMCBackground = (THnSparse*)dIsoGammaMC->Get("hIsoGammaIsovsM02vsPtnoCutsBackground");
      //Load actual data:
      THnSparse *PtIsovsM02vsPtData = (THnSparse*)dIsoGammaData->Get("hIsoGammaIsovsM02vsPtnoCuts");
    
      TAxis *axis = PtIsovsM02vsPtMC->GetAxis(2);
      int nbins = nPtBinEdges-1;//axis->GetNbins()/2000;

      TAxis *axisPtIso = PtIsovsM02vsPtMC->GetAxis(0);
      int nbinsIso = axisPtIso->GetNbins();

      //Arrays for holding values from the various regions of the PtIso vs. M02 distribution in MC and data.
      float AMC[nbins];
      float BMC[nbins];
      float CMC[nbins];
      float DMC[nbins];
      float AData[nbins];
      float BData[nbins];
      float CData[nbins];
      float DData[nbins];
      //Arrays for holding alfa and raw PABCD
      float alfa[nbins];
      float PABCDraw[nbins];
      float PABCD[nbins];

      //Loop over pT bins and calculate A,B,C and D:
      for(int bin = 1; bin<=nbins; ++bin){
      //Calculating the correction factor alfa
        //For projecting THnsparse
        TH2F *h2M02vsPtprojectionMC;
        TH2F *h2M02vsPtprojectionMCBackground;

        //TH2F *h2M02vsPtprojection= (TH2F*)Project2D(PtIsovsM02vsPt, binlowedge,binupedge);
        //Restrict pt-range
        PtIsovsM02vsPtMC->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        PtIsovsM02vsPtMCBackground->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        cout<<(PtIsovsM02vsPtMC->GetAxis(2)->GetName())<<"\n";
        h2M02vsPtprojectionMC = (TH2F*)PtIsovsM02vsPtMC->Projection(0,1);
        h2M02vsPtprojectionMCBackground = (TH2F*)PtIsovsM02vsPtMCBackground->Projection(0,1);
        //find bins
        int pTIsoMinBin = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(1.5);
        int pTIsoMaxBin = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(4.);
        int M02MinBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(0.1);
        int M02MinBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(0.3);
        int M02MaxBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(0.4);
        int M02MaxBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(2.0);

        AMC[bin-1]=h2M02vsPtprojectionMCBackground->Integral(1,pTIsoMinBin,M02MinBin1,M02MinBin2);
        BMC[bin-1]=h2M02vsPtprojectionMC->Integral(1,pTIsoMinBin,M02MaxBin1,M02MaxBin2);
        CMC[bin-1]=h2M02vsPtprojectionMC->Integral(pTIsoMaxBin,nbinsIso,M02MinBin1,M02MinBin2);
        DMC[bin-1]=h2M02vsPtprojectionMC->Integral(pTIsoMaxBin,nbinsIso,M02MaxBin1,M02MaxBin2);

        //Freeing memory:
        delete h2M02vsPtprojectionMC;
        delete h2M02vsPtprojectionMCBackground;

        cout<<"alfa:"<<"\n";
        cout<<((AMC[bin-1]/CMC[bin-1])/(BMC[bin-1]/DMC[bin-1]))<<"\n";

      }
      for(int bin = 1; bin<=nbins; ++bin){
        //For claculating the PABCD-raw:
        TH2F *h2M02vsPtprojectionData;
        //Identifying binning for pt-range (hardcoded granularity)
        float binlowedge = axis->GetBinLowEdge(bin);
        float binupedge = axis->GetBinUpEdge(2000);
        //Restrict Pt-range
        PtIsovsM02vsPtData->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        //Project onto pt-range:
        h2M02vsPtprojectionData = (TH2F*)PtIsovsM02vsPtData->Projection(0,1);

        //find bins
        int pTIsoMinBin = h2M02vsPtprojectionData->GetXaxis()->FindFixBin(1.5);
        int pTIsoMaxBin = h2M02vsPtprojectionData->GetXaxis()->FindFixBin(4.);
        int M02MinBin1 = h2M02vsPtprojectionData->GetYaxis()->FindFixBin(0.1);
        int M02MinBin2 = h2M02vsPtprojectionData->GetYaxis()->FindFixBin(0.3);
        int M02MaxBin1 = h2M02vsPtprojectionData->GetYaxis()->FindFixBin(0.4);
        int M02MaxBin2 = h2M02vsPtprojectionData->GetYaxis()->FindFixBin(2.0);

        //Integrate
        AData[bin-1]=h2M02vsPtprojectionData->Integral(1,pTIsoMinBin,M02MinBin1,M02MinBin2);
        BData[bin-1]=h2M02vsPtprojectionData->Integral(1,pTIsoMinBin,M02MaxBin1,M02MaxBin2);
        CData[bin-1]=h2M02vsPtprojectionData->Integral(pTIsoMaxBin,nbinsIso,M02MinBin1,M02MinBin2);
        DData[bin-1]=h2M02vsPtprojectionData->Integral(pTIsoMaxBin,nbinsIso,M02MaxBin1,M02MaxBin2);
        //Freeing memory:
        delete h2M02vsPtprojectionData;

        cout<<"PABCDraw:"<<"\n";
        cout<<((CData[bin-1]/AData[bin-1])/(DData[bin-1]/BData[bin-1]))<<"\n";
        PABCDraw[bin-1]=1-((CData[bin-1]/AData[bin-1])/(DData[bin-1]/BData[bin-1]));

        //Save:
        PABCD[bin-1]=1-((CData[bin-1]/AData[bin-1])/(DData[bin-1]/BData[bin-1]))*((AMC[bin-1]/CMC[bin-1])/(BMC[bin-1]/DMC[bin-1]));
        cout<<"PABCD"<<"\n";
        cout<<"Corrected:"<<PABCD[bin-1]<<" raw:"<<PABCDraw[bin-1]<<"\n";
      }
      auto c1 = new TCanvas("c1","PABCD corrected",200,10,700,500);
      TGraphErrors* PABCDplot = new TGraphErrors(nbins,pTbinsToplot,PABCD,0,0);
      PABCDplot->SetMarkerStyle(21);
      PABCDplot->SetMarkerSize(2);
      PABCDplot->SetMarkerColor(8);
      PABCDplot->SetTitle("P_{ABCD} vs p_{T};p_{T};P_{ABCD}");
      PABCDplot->Draw("AP");

      auto c2 = new TCanvas("c2","PABCDraw corrected",200,10,700,500);
      TGraphErrors* PABCDrawplot = new TGraphErrors(nbins,pTbinsToplot,PABCDraw,0,0);
      PABCDrawplot->SetMarkerStyle(21);
      PABCDrawplot->SetMarkerSize(2);
      PABCDrawplot->SetMarkerColor(8);
      PABCDrawplot->SetTitle("P_{ABCD,raw} vs p_{T};p_{T};P_{ABCD,raw}");
      PABCDrawplot->Draw("AP");

    }
  } 
}