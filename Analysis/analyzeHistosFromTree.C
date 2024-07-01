#include "Utilities.h"
#include "Binnings.h"
#include "analyseHistosFromTree.h"

void analyzeJets(TDirectory *dJets, GlobalOptions optns)
{

  if (optns.isMC)
  {
    TDirectory *dJESRSlices = dJets->mkdir("JESRSlices", "JESRSlices", kTRUE);
    dJESRSlices->cd();
    TH2F *hDLSubPLDivPLJetPtVsPLJetPt = (TH2F *)dJets->Get("hDLSubPLDivPLJetPtVsPLJetPt");
    TH1F *hJES = new TH1F("hJES", "hJES;pT,JES", nBins_JESR, binning_JESR);
    TH1F *hJER = new TH1F("hJER", "hJER;pT,JER", nBins_JESR, binning_JESR);
    for (int iJESRBin = 0; iJESRBin < nBins_JESR; iJESRBin++)
    {
      int pTBinFrom = hDLSubPLDivPLJetPtVsPLJetPt->GetYaxis()->FindBin(binning_JESR[iJESRBin] + 1E-5);
      int pTBinTo = hDLSubPLDivPLJetPtVsPLJetPt->GetYaxis()->FindBin(binning_JESR[iJESRBin + 1] - 1E-5);
      INFO(Form("%d: %d - %d", iJESRBin, pTBinFrom, pTBinTo))
      TH1F *hDLSubPLDivPLJetPtVsPLJetPtProjection = (TH1F *)hDLSubPLDivPLJetPtVsPLJetPt->ProjectionX(Form("hDLSubPLDivPLJetPtVsPLJetPtProjection_%.1f-%.1f", binning_JESR[iJESRBin], binning_JESR[iJESRBin + 1]), pTBinFrom, pTBinTo);
      hDLSubPLDivPLJetPtVsPLJetPtProjection->Write(Form("hDLSubPLDivPLJetPtVsPLJetPtProjection_%.1f-%.1f", binning_JESR[iJESRBin], binning_JESR[iJESRBin + 1]), TObject::kWriteDelete);
      hJES->SetBinContent(iJESRBin+1, hDLSubPLDivPLJetPtVsPLJetPtProjection->GetMean());
      hJES->SetBinError(iJESRBin+1, hDLSubPLDivPLJetPtVsPLJetPtProjection->GetMeanError());
      hJER->SetBinContent(iJESRBin+1, hDLSubPLDivPLJetPtVsPLJetPtProjection->GetStdDev());
      hJER->SetBinError(iJESRBin+1, hDLSubPLDivPLJetPtVsPLJetPtProjection->GetStdDevError());
    }
    dJets->cd();
    hJES->Write("hJES", TObject::kWriteDelete);
    hJER->Write("hJER", TObject::kWriteDelete);
  }
}

void CalculateEffPurABCD(string AnalysisDirectoryMC,string AnalysisDirectoryData, const char* outputDir, bool doABCD = true, bool doPurity = true, bool doAcceptance = true, bool doEfficiency = true)
{
  GlobalOptions optns(AnalysisDirectoryData, 1);

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

  TDirectory* dGammaMC = (TDirectory*)fInMC->Get("Gammas");
  TDirectory* dGammaData = (TDirectory*)fInData->Get("Gammas");

  if(doPurity || doEfficiency || doAcceptance || doABCD){
    //Outputdirectory:
    createDirectory(outputDir);
    //Output file:
    TFile *OutFile=(TFile*)TFile::Open(Form("%s/AnlysedHistosFromTree.root",outputDir),"RECREATE");
    //Purity and Efficiency from MC:
    if(doPurity || doEfficiency || doAcceptance){
      TH1F* hIsoGammaPt;
      TH1F* hGammaGenPt;
      TH1F* hIsoGammaPtSignal;
      //Read data:
      hIsoGammaPt = (TH1F*)dIsoGammaMC->Get("hIsoGammaPt")->Clone("hIsoGammaPt");
      hGammaGenPt = (TH1F*)dIsoGammaMC->Get("hGammaGenPt")->Clone("hGammaGenPt");
      hIsoGammaPtSignal = (TH1F*)dIsoGammaMC->Get("hIsoGammaPtSignal")->Clone("hIsoGammaPtSignal");

      if(doPurity){
        TH1F* hIsoGammaPurity;
        TH1F* hIsoGammaPtRebinned;
        //Rebin:
        hIsoGammaPtRebinned = (TH1F*)hIsoGammaPt->Rebin(npTbinEdges2-1, "hIsoGammaPtRebinned",pTbinEdges2);
        hIsoGammaPurity = (TH1F*)hIsoGammaPtSignal->Rebin(npTbinEdges2-1, "hPurity",pTbinEdges2);
        //Calculate Purity
        hIsoGammaPurity->Divide(hIsoGammaPtRebinned);
        hIsoGammaPurity->SetTitle("Purity");
        hIsoGammaPurity->SetName("Purity");
        //hIsoGammaPurity->Write();

        delete hIsoGammaPtRebinned;
      }
      if(doEfficiency){
        TH1F* hIsoGammaEfficiency;
        TH1F* hGammaGenPtRebinned;
        //Rebin:
        hGammaGenPtRebinned = (TH1F*)hGammaGenPt->Rebin(npTbinEdges2-1, "hGammGenPtRebinned",pTbinEdges2);
        hIsoGammaEfficiency = (TH1F*)hIsoGammaPt->Rebin(npTbinEdges2-1, "hEfficiency", pTbinEdges2);
        //Calculate efficiency
        hIsoGammaEfficiency->Divide(hGammaGenPtRebinned);
        hIsoGammaEfficiency->SetTitle("Efficiency");
        hIsoGammaEfficiency->SetName("Efficiency");
        //hIsoGammaEfficiency->Write();

        delete hGammaGenPtRebinned;
      }
      if(doAcceptance){
        //Initialise histograms
        TH1F* hIsoGammaAcceptance;
        TH1F* hGammaGenPtAcceptanceCut;
        TH1F* hGammaGenPtRebinned;
        //Read data:
        hGammaGenPtAcceptanceCut = (TH1F*)dIsoGammaMC->Get("hGammaGenPtAcceptanceCut")->Clone("hGammaGenPtAcceptanceCut");
        //Rebin:
        hIsoGammaAcceptance = (TH1F*)hGammaGenPtAcceptanceCut->Rebin(npTbinEdges2-1,"hAcceptance",pTbinEdges2);
        hGammaGenPtRebinned = (TH1F*)hGammaGenPt->Rebin(npTbinEdges2-1,"hGammaGenPtRebinned",pTbinEdges2);
        //Calculate acceptance:
        hIsoGammaAcceptance->Divide(hGammaGenPtRebinned);
        hIsoGammaAcceptance->SetTitle("Acceptance");
        hIsoGammaAcceptance->SetName("Acceptance");
        //hIsoGammaAcceptance->Write();

        delete hGammaGenPtAcceptanceCut;
        delete hGammaGenPtRebinned;
      }
    }
    if(doABCD){
      std::vector<ABCD> ABCDs;
      TDirectory *dABCD = DefineABCDHistos(OutFile,npTbinEdges-1, pTbinEdges);
      CalculateABCDData(ABCDs,dGammaData,dGammaMC,dABCD,npTbinEdges-1,pTbinEdges,pTIsoABCDCuts,M02ABCDCuts);
      FillABCDHistos(ABCDs,dABCD);
    }
    OutFile->Write();
    OutFile->Close();
    fInData->Close();
    fInMC->Close();
  }
} 

void analyzeHistosFromTree(TString AnalysisDirectory, string AnalysisDirectoryMC,string AnalysisDirectoryData, bool doABCD = true, bool doPurity = true, bool doAcceptance = true, bool doEfficiency = true)
{
  GlobalOptions optns(AnalysisDirectoryData, 1);

  if (optns.doJets)
  {
    TString inputFilePath = Form("%s/HistosFromTree.root", AnalysisDirectoryData.c_str());

    TFile *fIn = new TFile(inputFilePath, "UPDATE");
    if (!fIn)
      FATAL(Form("File %s not found", inputFilePath.Data()))

    TDirectory *dJets = (TDirectory *)fIn->Get("Jets");
    if (!dJets)
      FATAL("Dir Jets not found")

    analyzeJets(dJets, optns);
  }
  if (optns.doIsoGamma)
  {
    CalculateEffPurABCD(AnalysisDirectoryMC,AnalysisDirectoryData, AnalysisDirectory.Data(), doABCD, doPurity, doAcceptance, doEfficiency);
  }
}
