#include "Utilities.h"
#include "Binnings.h"

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

void analyzeHistosFromTree(TString AnalysisDirectory)
{

  GlobalOptions optns(AnalysisDirectory, 1);

  TString inputFilePath = Form("%s/HistosFromTree.root", AnalysisDirectory.Data());

  TFile *fIn = new TFile(inputFilePath, "UPDATE");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  if (optns.doJets)
  {
    TDirectory *dJets = (TDirectory *)fIn->Get("Jets");
    if (!dJets)
      FATAL("Dir Jets not found")

    analyzeJets(dJets, optns);
  }
}