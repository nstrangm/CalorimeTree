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
  
  TDirectory* dClusterMC = (TDirectory*)fInMC->Get("Clusters");

  TDirectory* dGammaMC = (TDirectory*)fInMC->Get("Gammas");
  TDirectory* dGammaData = (TDirectory*)fInData->Get("Gammas");

  if(doPurity || doEfficiency || doAcceptance || doABCD){
    //Outputdirectory:
    createDirectory(outputDir);
    //Output file:
    TFile *OutFile=(TFile*)TFile::Open(Form("%s/AnlysedHistosFromTree.root",outputDir),"RECREATE");
    //Purity and Efficiency from MC:
    if(doPurity || doEfficiency || doAcceptance){
      //Initiate relevant histograms:
      //Generator level:
      TH1F* NGen;
      TH1F* NGenAcceptanceCut;
      TH1F* NGenSignal;
      TH1F* rawNGen;
      TH1F* rawNGenAcceptanceCut;
      TH1F* rawNGenSignal;
      //MC reconstructed:
      TH1F* NSignalClustCut;
      TH1F* NClustCut;
      TH1F* NSignalClustIsoCut;
      TH1F* NClustIsoCut;
      TH1F* rawNSignalClustCut;
      TH1F* rawNClustCut;
      TH1F* rawNSignalClustIsoCut;
      TH1F* rawNClustIsoCut;

      //Read relevant histograms:
      //Generator level:
      rawNGen = (TH1F*)dGammaMC->Get("hGammaGenPt");
      NGen = (TH1F*)rawNGen->Rebin(npTbinEdges2-1, "NGen",pTbinEdges2);
      rawNGenSignal = (TH1F*)dIsoGammaMC->Get("hGammaGenPtSignal");
      NGenSignal = (TH1F*)rawNGenSignal->Rebin(npTbinEdges2-1, "NGenSignal",pTbinEdges2);
      rawNGenAcceptanceCut = (TH1F*)dGammaMC->Get("hGammaGenPtAcceptanceCut");
      NGenAcceptanceCut = (TH1F*)rawNGenAcceptanceCut->Rebin(npTbinEdges2-1, "NGenAcceptanceCut",pTbinEdges2);
      //MC reconstructed:
      rawNSignalClustCut = (TH1F*)dGammaMC->Get("hIsoGammaPtSignal");
      NSignalClustCut= (TH1F*)rawNSignalClustCut->Rebin(npTbinEdges2-1, "NSignalClustCut",pTbinEdges2);
      rawNClustCut = (TH1F*)dGammaMC->Get("hIsoGammaPt");
      NClustCut = (TH1F*)rawNClustCut->Rebin(npTbinEdges2-1, "NClustCut",pTbinEdges2);
      rawNSignalClustIsoCut = (TH1F*)dIsoGammaMC->Get("hIsoGammaPtSignal");
      NSignalClustIsoCut = (TH1F*)rawNSignalClustIsoCut->Rebin(npTbinEdges2-1, "NSignalClustIsoCut", pTbinEdges2);
      rawNClustIsoCut = (TH1F*)dIsoGammaMC->Get("hIsoGammaPt");
      NClustIsoCut = (TH1F*)rawNClustIsoCut->Rebin(npTbinEdges2-1, "NClustIsoCut",pTbinEdges2);

      if(doPurity){
        TH1F* Purity;
        TH1F* PurityClusterCuts;

        Purity = (TH1F*)NSignalClustIsoCut->Clone("Purity");
        Purity->Divide(NClustIsoCut);
      //  Purity->Write();

        PurityClusterCuts = (TH1F*)NSignalClustCut->Clone("PurityClusterCuts");
        PurityClusterCuts->Divide(NClustCut);
        for(int i=0;i<npTbinEdges2-1;i++){
          cout<<"pT-bin:"<<pTbinEdges2[i]<<";"<<pTbinEdges2[i+1]<<"\n";
          cout<<"Integral:"<<"\n";
          cout<<(NSignalClustIsoCut->Integral(i,i+1)/NClustIsoCut->Integral(i,i+1));
        }
      //  PurityClusterCuts->Write();
      }
      if(doAcceptance){
        TH1F* Acceptance;
        Acceptance = (TH1F*)NGenAcceptanceCut->Clone("Acceptance");
        Acceptance->Divide(NGen);
      //  Acceptance->Write();
      }
      if(doEfficiency){
        TH1F* FullEfficiency;
        TH1F* EfficiencyClusterCuts;
        TH1F* EfficiencyClusterIsolationCuts;

        FullEfficiency = (TH1F*)NSignalClustIsoCut->Clone("FullEfficiency");
        FullEfficiency->Divide(NGenSignal);
      //  FullEfficiency->Write();

        EfficiencyClusterCuts = (TH1F*)NSignalClustCut->Clone("EfficiencyClusterCuts");
        EfficiencyClusterCuts->Divide(NGenSignal);
      //  EfficiencyClusterCuts->Write();

        EfficiencyClusterIsolationCuts = (TH1F*)NSignalClustIsoCut->Clone("EfficiencyClusterIsolationCuts");
        EfficiencyClusterIsolationCuts->Divide(NSignalClustCut);
       // EfficiencyClusterIsolationCuts->Write();

        

      }

      delete NGen;
      delete NGenAcceptanceCut;
      delete NGenSignal;
      delete NSignalClustCut;
      delete NClustCut;
      delete NSignalClustIsoCut;
      delete NClustIsoCut;
      delete rawNGen;
      delete rawNGenAcceptanceCut;
      delete rawNGenSignal;
      delete rawNSignalClustCut;
      delete rawNClustCut;
      delete rawNSignalClustIsoCut;
      delete rawNClustIsoCut;

      
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