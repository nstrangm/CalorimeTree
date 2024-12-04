#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "analyseExclGammaJet.h"


void determineNTriggers(TDirectory *dSignalPhoton, TDirectory *dReferencePhoton, TDirectory *dSignalMerged, TDirectory *dReferenceMerged)
{ 
    // get hIsoGammaPt
    TH1F *hIsoGammaPtPhotonSignal = (TH1F*)dSignalPhoton->Get("hIsoGammaPt");
    hIsoGammaPtPhotonSignal->SetName("hIsoGammaPtPhotonSignal");
    TH1F *hIsoGammaPtPhotonReference = (TH1F*)dReferencePhoton->Get("hIsoGammaPt");
    hIsoGammaPtPhotonReference->SetName("hIsoGammaPtPhotonReference");

    // get number of entries
    nTriggersSignalPhoton = hIsoGammaPtPhotonSignal->Integral();
    nTriggersReferencePhoton = hIsoGammaPtPhotonReference->Integral();

    TH1F* hIsoGammaPtMergedSignal = (TH1F*)dSignalMerged->Get("hIsoGammaPt");
    hIsoGammaPtMergedSignal->SetName("hIsoGammaPtMergedSignal");
    TH1F* hIsoGammaPtMergedReference = (TH1F*)dReferenceMerged->Get("hIsoGammaPt");
    hIsoGammaPtMergedReference->SetName("hIsoGammaPtMergedReference");

    nTriggersSignalMerged = hIsoGammaPtMergedSignal->Integral();
    nTriggersReferenceMerged = hIsoGammaPtMergedReference->Integral();
}

// returns a histogram shifted by x
template<class T>
TH1F* shiftTH1ByX(T *h, double x, TString name)
{
   // create a new histogram that has same amount of bins but min and max are shifted by x
   double min = h->GetXaxis()->GetXmin();
   double max = h->GetXaxis()->GetXmax();
   int nbins = h->GetNbinsX();

   // we need to keep asme bin ranges otherwise divide will break later
   TH1F *hShifted = new TH1F(name, name, nbins, min, max);
   
   // loop over all bins and copy the content of h - x to hShifted
   for(int i = 1; i <= nbins; i++) {
     double binContent = h->GetBinContent(i);
     double binError = h->GetBinError(i);
     double binCenter = h->GetBinCenter(i);
     double shiftedBinCenter = binCenter + x;
     
     // Find bin in shifted histogram corresponding to shiftedBinCenter
     int shiftedBin = hShifted->FindBin(shiftedBinCenter);
     
     // Only fill if bin is valid and shifted value is >= 0
     if(shiftedBin > 0 && shiftedBin <= nbins && shiftedBinCenter >= 0) {
       hShifted->SetBinContent(shiftedBin, binContent);
       hShifted->SetBinError(shiftedBin, binError);
     }
   }

   hShifted->SetTitle(Form("Rho distribution shifted by %.1f GeV", x));

   return hShifted;
}

void getRhoDistributions(TDirectory *dSignal, TDirectory *dReference, TString label)
{
  // get rho vs pt
  // create new directory in fOut
  double nTriggersSignal = 0;
  double nTriggersReference = 0;
  if(label == "Photons")
  {
    nTriggersSignal = nTriggersSignalPhoton;
    nTriggersReference = nTriggersReferencePhoton;
  } 
  else if (label == "MergedPi0s")
  {
    nTriggersSignal = nTriggersSignalMerged;
    nTriggersReference = nTriggersReferenceMerged;
  }
  TDirectory *dRho = fOut->mkdir(Form("Rho_%s", label.Data()));

  TH2F *hRhoSignalPt = (TH2F*)dSignal->Get("hIsoGammaRhoPt");
  hRhoSignalPt->SetName("hIsoGammaRhoPtSignal");
  TH2F *hRhoReferencePt = (TH2F*)dReference->Get("hIsoGammaRhoPt");
  hRhoReferencePt->SetName("hIsoGammaRhoPtReference");
  
  hRhoSignalPt->RebinX(4);
  hRhoReferencePt->RebinX(4);
  double dPt = hRhoSignalPt->GetYaxis()->GetBinWidth(1);
  hRhoSignalPt->Scale(1./(nTriggersSignal*dPt));
  hRhoReferencePt->Scale(1./(nTriggersReference*dPt));

  // project them onto the x axis
  TH1F *hRhoSignal = (TH1F*)hRhoSignalPt->ProjectionX();
  hRhoSignal->SetName("hRhoSignal");
  TH1F *hRhoReference = (TH1F*)hRhoReferencePt->ProjectionX(); 
  hRhoReference->SetName("hRhoReference");

  // create many shifted histograms
  // create vector with possible shifts from -1 to 3 GeV in steps of 0.1 GeV
  std::vector<double> ptShifts;
  for(double x = -1; x <= 5; x+=0.1)
  {
    ptShifts.push_back(x);
  }
  std::vector<TH1F*> hRhoReferenceShifted;
  for(double x : ptShifts)
  {
    hRhoReferenceShifted.push_back((TH1F*)shiftTH1ByX(hRhoReference, x, Form("hRhoReferenceShifted_%.1f", x)));
  }

  // determine ratio of hRhoSignal and hRhoReferenceShifted
  std::vector<TH1F*> hRhoRatios;
  std::vector<TF1*> fRhoRatiosFit;
  for(int i = 0; i < hRhoReferenceShifted.size(); i++)
  {
    TH1F* hShifted = hRhoReferenceShifted[i];
    TH1F* hRatio = (TH1F*)hRhoSignal->Clone(Form("hRhoRatio_%.1f", ptShifts[i]));
    hRatio->Divide(hShifted,hRhoSignal,1,1);
    hRhoRatios.push_back(hRatio);
    
    TF1* fit = new TF1(Form("fRhoRatiosFit_%f", ptShifts[i]), "pol1", hRatio->GetXaxis()->GetXmin(), hRatio->GetXaxis()->GetXmax());
    fRhoRatiosFit.push_back(fit);
  }
  double smallestSlope = 1000000;
  int bestShiftIndex = 0;
  for(int i=0; i<hRhoRatios.size(); i++)
  {
    hRhoRatios[i]->Fit(fRhoRatiosFit[i], "R0");
    double slope = fRhoRatiosFit[i]->GetParameter(1);
    if(TMath::Abs(slope) < smallestSlope)
    {
      smallestSlope = TMath::Abs(slope);
      bestShiftIndex = i;
    }
  }

  dRho->cd();

  // write all histograms of interest to dRho
  hRhoSignal->Write("hRhoSignal");
  hRhoReference->Write("hRhoReference");

  // write all hRhoRatios and fRhoRatiosFit to dRho
  for(int i=0; i<hRhoRatios.size(); i++)
  {
    hRhoRatios[i]->Write(Form("hRhoRatio_%.1f", ptShifts[i]));
    fRhoRatiosFit[i]->Write(Form("fRhoRatiosFit_%.1f", ptShifts[i]));
  }
  
  // write hRhoReferenceShifted
  for(int i=0; i<hRhoReferenceShifted.size(); i++)
  {
    hRhoReferenceShifted[i]->Write(Form("hRhoReferenceShifted_%.1f", ptShifts[i]));
  }
  INFO(Form("Best shift index: %d with shift by %.1f GeV", bestShiftIndex, ptShifts[bestShiftIndex]));
  hRhoRatios[bestShiftIndex]->Write("hRhoRatioBestShift");
  fRhoRatiosFit[bestShiftIndex]->Write("fRhoRatiosFitBestShift");


}

void getRecoilJetPtDistributions(TDirectory *dSignal, TDirectory *dReference, TString label)
{

    // choose nTriggersSignal and nTriggersReference depending on label
    double nTriggersSignal = 0;
    double nTriggersReference = 0;
    if(label == "Photons")
    {
      nTriggersSignal = nTriggersSignalPhoton;
      nTriggersReference = nTriggersReferencePhoton;
    }
    else if(label == "MergedPi0s")
    {
      nTriggersSignal = nTriggersSignalMerged;
      nTriggersReference = nTriggersReferenceMerged;
    }
    // create outputfolder for recoil jet
    TDirectory *dRecoilJet = fOut->mkdir(Form("RecoilJet_%s", label.Data()));

    // get THnSparseF for signal and reference
    THnSparseF *hIsoGammaJetDeltaPhiJetPtGammaPt = (THnSparseF*)dSignal->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
    hIsoGammaJetDeltaPhiJetPtGammaPt->SetName("hIsoGammaJetDeltaPhiJetPtGammaPtSignal");
    THnSparseF *hIsoGammaJetDeltaPhiJetPtGammaPtReference = (THnSparseF*)dReference->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
    hIsoGammaJetDeltaPhiJetPtGammaPtReference->SetName("hIsoGammaJetDeltaPhiJetPtGammaPtReference");

    // do projection onto Jet and DeltaPhi axis
    TH2F *hIsoGammaJetDeltaPhiJetPtSignal = (TH2F*)hIsoGammaJetDeltaPhiJetPtGammaPt->Projection(1,0);
    hIsoGammaJetDeltaPhiJetPtSignal->SetName("hIsoGammaJetDeltaPhiJetPtSignal");
    TH2F *hIsoGammaJetDeltaPhiJetPtReference = (TH2F*)hIsoGammaJetDeltaPhiJetPtGammaPtReference->Projection(1,0);
    hIsoGammaJetDeltaPhiJetPtReference->SetName("hIsoGammaJetDeltaPhiJetPtReference");

    // normalize by number of triggers 
    // rebin in pt by factor 4 in y axis
    hIsoGammaJetDeltaPhiJetPtSignal->RebinY(4);
    hIsoGammaJetDeltaPhiJetPtReference->RebinY(4);

    double dPt = hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->GetBinWidth(1);
    double dPhi = hIsoGammaJetDeltaPhiJetPtSignal->GetXaxis()->GetBinWidth(1);
    hIsoGammaJetDeltaPhiJetPtSignal->Scale(1./(nTriggersSignal*dPt*dPhi));
    hIsoGammaJetDeltaPhiJetPtReference->Scale(1./(nTriggersReference*dPt*dPhi));
    
    hIsoGammaJetDeltaPhiJetPtSignal->GetXaxis()->SetTitle("#Delta#phi");
    hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->SetTitle("#it{p}_{T, ch jet}^{reco} (GeV/#it{c})");
    hIsoGammaJetDeltaPhiJetPtSignal->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi}");


    // write all histograms of interest to dRecoilJet
    dRecoilJet->cd();
    hIsoGammaJetDeltaPhiJetPtSignal->Write("hIsoGammaJetDeltaPhiJetPtSignal");
    hIsoGammaJetDeltaPhiJetPtReference->Write("hIsoGammaJetDeltaPhiJetPtReference");

    // project the two dimensional distributions ontot the jet pt axis
    // do projections in phi bins
    for(auto &phiBin : phiBins)
    {
      hIsoGammaJetDeltaPhiJetPtSignal->GetXaxis()->SetRangeUser(phiBin.first, phiBin.second);
      hIsoGammaJetDeltaPhiJetPtReference->GetXaxis()->SetRangeUser(phiBin.first, phiBin.second);
      TH1F* hRecoilJetPtSignal = (TH1F*) hIsoGammaJetDeltaPhiJetPtSignal->ProjectionY();
      TH1F* hRecoilJetPtReference = (TH1F*)hIsoGammaJetDeltaPhiJetPtReference->ProjectionY();
      hRecoilJetPtSignal->SetName(Form("hRecoilJetPtSignal_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtReference->SetName(Form("hRecoilJetPtReference_%.2f_%.2f", phiBin.first, phiBin.second));

      // get ratio of signal over reference
      TH1F* hRecoilJetPtRatio = (TH1F*)hRecoilJetPtSignal->Clone(Form("hRecoilJetPtRatio_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtRatio->Divide(hRecoilJetPtReference);
      hRecoilJetPtRatio->SetName(Form("hRecoilJetPtSignalReferenceRatio_%.2f_%.2f", phiBin.first, phiBin.second));

      // determine scaling factor for reference in range cRefFitRange
      TF1* fCRefFit = new TF1(Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second), "pol0", cRefFitRange.first, cRefFitRange.second);
      hRecoilJetPtRatio->Fit(fCRefFit, "R0");
      double cRef = fCRefFit->GetParameter(0);
      INFO(Form("Determined scaling factor for reference in range %.2f_%.2f: %.3f", cRefFitRange.first, cRefFitRange.second, cRef));
      
      TH1F* hRecoilJetReferenceScaled = (TH1F*)hRecoilJetPtReference->Clone(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetReferenceScaled->Scale(1./cRef);
      hRecoilJetReferenceScaled->SetName(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      
      TH1F* hRecoilJetPtSignal_Subtracted = (TH1F*)hRecoilJetPtSignal->Clone(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted->Add(hRecoilJetReferenceScaled, -1);
      hRecoilJetPtSignal_Subtracted->SetName(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));
      
      // write all histograms of interest to dRecoilJet
      dRecoilJet->cd();
      hRecoilJetPtSignal->Write(Form("hRecoilJetPtSignal_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtReference->Write(Form("hRecoilJetPtReference_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtRatio->Write(Form("hRecoilJetPtSignalReferenceRatio_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetReferenceScaled->Write(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted->Write(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));
      fCRefFit->Write(Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second));
    }

}

std::pair<double, double> calculateABCDPurity(THnSparseF* hM02vsIsoPt, double ptMin, double ptMax)
{
    hM02vsIsoPt->Sumw2();

    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();
    // show the axis ranges
    INFO(Form("Axis 0: %.2f_%.2f", hM02vsIsoPt->GetAxis(0)->GetXmin(), hM02vsIsoPt->GetAxis(0)->GetXmax()));
    INFO(Form("Axis 1: %.2f_%.2f", hM02vsIsoPt->GetAxis(1)->GetXmin(), hM02vsIsoPt->GetAxis(1)->GetXmax()));
    INFO(Form("Axis 2: %.2f_%.2f", hM02vsIsoPt->GetAxis(2)->GetXmin(), hM02vsIsoPt->GetAxis(2)->GetXmax()));
    
    // find bins for region A
    INFO(Form("Region A: %.2f_%.2f, %.2f_%.2f", regionA.minM02, regionA.maxM02, regionA.minIso, regionA.maxIso));
    std::tuple<int, int, int, int> binsA = {hM02vsIsoPt->GetAxis(1)->FindBin(regionA.minM02), hM02vsIsoPt->GetAxis(1)->FindBin(regionA.maxM02), hM02vsIsoPt->GetAxis(0)->FindBin(regionA.minIso), hM02vsIsoPt->GetAxis(0)->FindBin(regionA.maxIso)};
  
    // Integrate over region B
    std::tuple<int, int, int, int> binsB = {hM02vsIsoPt->GetAxis(1)->FindBin(regionB.minM02), hM02vsIsoPt->GetAxis(1)->FindBin(regionB.maxM02), hM02vsIsoPt->GetAxis(0)->FindBin(regionB.minIso), hM02vsIsoPt->GetAxis(0)->FindBin(regionB.maxIso)};
    INFO(Form("Region B: %.2f_%.2f, %.2f_%.2f", regionB.minM02, regionB.maxM02, regionB.minIso, regionB.maxIso));
    // Integrate over region C
    std::tuple<int, int, int, int> binsC = {hM02vsIsoPt->GetAxis(1)->FindBin(regionC.minM02), hM02vsIsoPt->GetAxis(1)->FindBin(regionC.maxM02), hM02vsIsoPt->GetAxis(0)->FindBin(regionC.minIso), hM02vsIsoPt->GetAxis(0)->FindBin(regionC.maxIso)};
    INFO(Form("Region C: %.2f_%.2f, %.2f_%.2f", regionC.minM02, regionC.maxM02, regionC.minIso, regionC.maxIso));
    // Integrate over region D
    std::tuple<int, int, int, int> binsD = {hM02vsIsoPt->GetAxis(1)->FindBin(regionD.minM02), hM02vsIsoPt->GetAxis(1)->FindBin(regionD.maxM02), hM02vsIsoPt->GetAxis(0)->FindBin(regionD.minIso), hM02vsIsoPt->GetAxis(0)->FindBin(regionD.maxIso)};
    INFO(Form("Region D: %.2f_%.2f, %.2f_%.2f", regionD.minM02, regionD.maxM02, regionD.minIso, regionD.maxIso));

    TH1D* hA;
    TH1D* hB;
    TH1D* hC;
    TH1D* hD;
    
    hM02vsIsoPt->GetAxis(1)->SetRange(std::get<0>(binsA), std::get<1>(binsA));
    hM02vsIsoPt->GetAxis(0)->SetRange(std::get<2>(binsA), std::get<3>(binsA));
    hA = (TH1D*)hM02vsIsoPt->Projection(2,"E");
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->SetRange(std::get<0>(binsB), std::get<1>(binsB));
    hM02vsIsoPt->GetAxis(0)->SetRange(std::get<2>(binsB), std::get<3>(binsB));
    hB = (TH1D*)hM02vsIsoPt->Projection(2,"E");
    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->SetRange(std::get<0>(binsC), std::get<1>(binsC));
    hM02vsIsoPt->GetAxis(0)->SetRange(std::get<2>(binsC), std::get<3>(binsC));
    hC = (TH1D*)hM02vsIsoPt->Projection(2,"E");
    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->SetRange(std::get<0>(binsD), std::get<1>(binsD));
    hM02vsIsoPt->GetAxis(0)->SetRange(std::get<2>(binsD), std::get<3>(binsD));
    hD = (TH1D*)hM02vsIsoPt->Projection(2,"E");
    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();

    double nAErr;
    double nA = hA->IntegralAndError(hA->GetXaxis()->FindBin(ptMin), hA->GetXaxis()->FindBin(ptMax), nAErr);
    double nBErr;
    double nB = hB->IntegralAndError(hB->GetXaxis()->FindBin(ptMin), hB->GetXaxis()->FindBin(ptMax), nBErr);
    double nCErr;
    double nC = hC->IntegralAndError(hC->GetXaxis()->FindBin(ptMin), hC->GetXaxis()->FindBin(ptMax), nCErr);
    double nDErr;
    double nD = hD->IntegralAndError(hD->GetXaxis()->FindBin(ptMin), hD->GetXaxis()->FindBin(ptMax), nDErr);



    // calculate purity
    double purity = 0;
    double purityErr = 0;
    if((nA!=0) && (nB!=0) && (nD!=0)){
      purity = 1 - ((nC/nA)/(nD/nB));
      purityErr = TMath::Sqrt(((nA * nA * nB * nB * nC * nC * nDErr * nDErr) + (nA * nA * nB * nB * nD * nD * nCErr * nCErr) + (nA * nA * nC * nC * nD * nD * nBErr * nBErr) + (nAErr * nAErr * nB * nB * nC * nC * nD * nD)) / (pow(nA, 4) * pow(nD, 4)));
    }


    INFO(Form("Purity: %.3f #pm %.3f", purity, purityErr));
    INFO(Form("\t nA: %.3f #pm %.3f", nA, nAErr));
    INFO(Form("\t nB: %.3f #pm %.3f", nB, nBErr));
    INFO(Form("\t nC: %.3f #pm %.3f", nC, nCErr));
    INFO(Form("\t nD: %.3f #pm %.3f", nD, nDErr));
    return std::make_pair(purity, purityErr);

}



void processPurity(TDirectory *dSignal,TString label, GlobalOptions optns)
{
    ExclusiveTriggerParticleSelection exclTrigSelection(optns);
    TDirectory *dPurity = fOut->mkdir(Form("Purity_%s", label.Data()));
    THnSparseF *hM02vsIsoPt = (THnSparseF*)dSignal->Get("hIsoGammaIsovsM02vsPt");
    hM02vsIsoPt->SetName("hM02vsIsoPt");

    // Set range of pt axis to be used for purity calculation
    double ptMin = exclTrigSelection.getPtMinSignal();
    double ptMax = exclTrigSelection.getPtMaxSignal();
    hM02vsIsoPt->GetAxis(2)->SetRangeUser(ptMin, ptMax);
    TH2F* hM02vsIsoSignal = (TH2F*)hM02vsIsoPt->Projection(0, 1);
    hM02vsIsoSignal->SetName("hM02vsIsoSignal");
    
    ptMin = exclTrigSelection.getPtMinReference();  
    ptMax = exclTrigSelection.getPtMaxReference();
    hM02vsIsoPt->GetAxis(2)->SetRangeUser(ptMin, ptMax);
    TH2F* hM02vsIsoReference = (TH2F*)hM02vsIsoPt->Projection(0, 1);
    hM02vsIsoReference->SetName("hM02vsIsoReference");

    hM02vsIsoPt->GetAxis(0)->UnZoom();
    hM02vsIsoPt->GetAxis(1)->UnZoom();
    hM02vsIsoPt->GetAxis(2)->UnZoom();

    // calculate purity
    std::pair<double, double> puritySignal = calculateABCDPurity(hM02vsIsoPt,exclTrigSelection.getPtMinSignal(), exclTrigSelection.getPtMaxSignal());
    std::pair<double, double> purityReference = calculateABCDPurity(hM02vsIsoPt,exclTrigSelection.getPtMinReference(), exclTrigSelection.getPtMaxReference());

    // write output to dPurity
    dPurity->cd();
    hM02vsIsoSignal->Write("hM02vsIsoSignal");
    hM02vsIsoReference->Write("hM02vsIsoReference");
    TF1* fPuritySignal = new TF1("fPuritySignal", "pol0", 0, 1);
    fPuritySignal->SetParameter(0, puritySignal.first);
    fPuritySignal->SetParError(0, puritySignal.second);
    fPuritySignal->Write("fPuritySignal");
    TF1* fPurityReference = new TF1("fPurityReference", "pol0", 0, 1);
    fPurityReference->SetParameter(0, purityReference.first);
    fPurityReference->SetParError(0, purityReference.second);
    fPurityReference->Write("fPurityReference");
}




void analyseExclGammaJet(TString AnalysisDirectory, bool isDebugRun = false)
{
  GlobalOptions optns(AnalysisDirectory, isDebugRun);



  TString inputFilePath = Form("%s/HistosFromTree.root", AnalysisDirectory.Data());
  fOut = new TFile(Form("%s/ExclGammaJet.root", AnalysisDirectory.Data()), "RECREATE");

  if (!std::filesystem::exists(inputFilePath.Data()))
  {
    inputFilePath = Form("%s/HistosFromTree_0.root", AnalysisDirectory.Data());
    INFO(Form("Did not find HistosFromTree.root, looking for %s now.", inputFilePath.Data()))
  }

  TFile *fIn = new TFile(inputFilePath, "READ");
  if (!fIn)
    FATAL(Form("File %s not found", inputFilePath.Data()))

  // get all required input folders

  // dGammas
  TDirectory *dGammas = (TDirectory *)fIn->Get("Gammas");
  // trigger photons
  TDirectory *dTriggerPhotonsSignal = (TDirectory *)fIn->Get("TriggerPhotonsSignal");
  TDirectory *dTriggerPhotonsReference = (TDirectory *)fIn->Get("TriggerPhotonsReference");
  TDirectory *dTriggerPhotonsSignalQA = nullptr;
  TDirectory *dTriggerPhotonsReferenceQA = nullptr;
  if(optns.doQA)
  {
    dTriggerPhotonsSignalQA = (TDirectory *)fIn->Get("TriggerPhotonSignalQA");
    dTriggerPhotonsReferenceQA = (TDirectory *)fIn->Get("TriggerPhotonReferenceQA");
  }

  // trigger merged pi0s
  TDirectory *dTriggerMergedPi0sSignal = (TDirectory *)fIn->Get("TriggerMergedPi0sSignal");
  TDirectory *dTriggerMergedPi0sReference = (TDirectory *)fIn->Get("TriggerMergedPi0sReference");
  TDirectory *dTriggerMergedPi0sSignalQA = nullptr;
  TDirectory *dTriggerMergedPi0sReferenceQA = nullptr;
  if(optns.doQA)
  {
    dTriggerMergedPi0sSignalQA = (TDirectory *)fIn->Get("TriggerMergedPi0sSignalQA");
    dTriggerMergedPi0sReferenceQA = (TDirectory *)fIn->Get("TriggerMergedPi0sReferenceQA");
  }

  // trigger photon - jet correlations
  TDirectory *dTriggerGammaJetCorrelationsSignal = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationsSignal");
  TDirectory *dTriggerGammaJetCorrelationsReference = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationsReference");
  TDirectory *dTriggerGammaJetCorrelationsSignalQA = nullptr;
  TDirectory *dTriggerGammaJetCorrelationsReferenceQA = nullptr;
  if(optns.doQA)
  {
    dTriggerGammaJetCorrelationsSignalQA = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationSignalQA");
    dTriggerGammaJetCorrelationsReferenceQA = (TDirectory *)fIn->Get("TriggerGammaJetCorrelationReferenceQA");
  }

  // trigger merged pi0 - jet correlations
  TDirectory *dTriggerMergedPi0JetCorrelationsSignal = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationsSignal");
  TDirectory *dTriggerMergedPi0JetCorrelationsReference = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationsReference");
  TDirectory *dTriggerMergedPi0JetCorrelationsSignalQA = nullptr;
  TDirectory *dTriggerMergedPi0JetCorrelationsReferenceQA = nullptr;
  if(optns.doQA)
  {
    dTriggerMergedPi0JetCorrelationsSignalQA = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationSignalQA");
    dTriggerMergedPi0JetCorrelationsReferenceQA = (TDirectory *)fIn->Get("TriggerMergedPi0JetCorrelationReferenceQA");
  }

  // check which directories exist and print those that do not
  if(!dTriggerPhotonsSignal)
    INFO("TriggerPhotonsSignal not found");
  if(!dTriggerPhotonsReference)
    INFO("TriggerPhotonsReference not found");
  if(!dTriggerMergedPi0sSignal)
    INFO("TriggerMergedPi0sSignal not found");
  if(!dTriggerMergedPi0sReference)
    INFO("TriggerMergedPi0sReference not found");
  if(!dTriggerGammaJetCorrelationsSignal)
    INFO("TriggerGammaJetCorrelationsSignal not found");
  if(!dTriggerGammaJetCorrelationsReference)
    INFO("TriggerGammaJetCorrelationsReference not found");
  if(!dTriggerPhotonsSignalQA)
    INFO("TriggerPhotonsSignalQA not found");
  if(!dTriggerPhotonsReferenceQA)
    INFO("TriggerPhotonsReferenceQA not found");
  if(!dTriggerMergedPi0sSignalQA)
    INFO("TriggerMergedPi0sSignalQA not found");
  if(!dTriggerMergedPi0sReferenceQA)
    INFO("TriggerMergedPi0sReferenceQA not found");
  if(!dTriggerGammaJetCorrelationsSignalQA)
    INFO("TriggerGammaJetCorrelationsSignalQA not found");
  if(!dTriggerGammaJetCorrelationsReferenceQA)
    INFO("TriggerGammaJetCorrelationsReferenceQA not found");
  if(!dTriggerMergedPi0JetCorrelationsSignalQA)
    INFO("TriggerMergedPi0JetCorrelationsSignalQA not found");
  if(!dTriggerMergedPi0JetCorrelationsReferenceQA)
    INFO("TriggerMergedPi0JetCorrelationsReferenceQA not found");


  // determine number of triggers
  determineNTriggers(dTriggerPhotonsSignal, dTriggerPhotonsReference, dTriggerMergedPi0sSignal, dTriggerMergedPi0sReference);
  INFO(Form("Determined number of photon triggers: Signal = %f, Reference = %f", nTriggersSignalPhoton, nTriggersReferencePhoton));
  INFO(Form("Determined number of merged pi0 triggers: Signal = %f, Reference = %f", nTriggersSignalMerged, nTriggersReferenceMerged));
  getRhoDistributions(dTriggerPhotonsSignalQA, dTriggerPhotonsReferenceQA, "Photons");
  getRhoDistributions(dTriggerMergedPi0sSignalQA, dTriggerMergedPi0sReferenceQA, "MergedPi0s");

  getRecoilJetPtDistributions(dTriggerGammaJetCorrelationsSignal, dTriggerGammaJetCorrelationsReference, "Photons");
  getRecoilJetPtDistributions(dTriggerMergedPi0JetCorrelationsSignal, dTriggerMergedPi0JetCorrelationsReference, "MergedPi0s");

  // purity calculation needs to be called with none isolated photons
  processPurity(dGammas,"Photons", optns);


}