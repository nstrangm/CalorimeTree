#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "analyseExclGammaJet.h"
#include "TVirtualFitter.h"
#include "TMath.h"

// takes as input a ratio histogram and returns the optimal fit function and used range
// it begins using the fitting range specified in cRefFitRange
// going from the upper edge it goes in steps of 1 GeV and performs a pol1 fit until it finds a range where the fit is stable
std::pair<TF1*, std::pair<double, double>> findOptimalCRefFit(TH1F* ratio, TString name){
    double ptMax = cRefFitRange.second;
    double ptMin = cRefFitRange.first;
    if(ptMax<ptMin){
        FATAL("ptMax is less than ptMin");
    }
    while(ptMax - ptMin > minFitInterval){
        TF1* fCRefFitTmp = new TF1(Form("fCRefFitTmp_%s", name.Data()), "pol1", ptMin, ptMax);
        ratio->Fit(fCRefFitTmp, "R0");
        double cRef = fCRefFitTmp->GetParameter(0);
        double cRefError = fCRefFitTmp->GetParError(0);
        double cRefSlope = fCRefFitTmp->GetParameter(1);
        double cRefSlopeError = fCRefFitTmp->GetParError(1);
        // check if slope is consistent with 0 within errors
        bool isConsistent = (TMath::Abs(cRefSlope) < cRefSlopeError);
        delete fCRefFitTmp; // Delete the temporary fit object to prevent memory leak
        if(isConsistent){
            break;
        }
        ptMax -= 1;
    }

    // do fit one more time with final fit range
    TF1* fCRefFit = new TF1(Form("%s", name.Data()), "pol1", ptMin, ptMax);
    ratio->Fit(fCRefFit, "R0");
    // check one more time if slope is consistent with 0, otherwise throw warning
    double cRefSlope = fCRefFit->GetParameter(1);
    double cRefSlopeError = fCRefFit->GetParError(1);
    if(TMath::Abs(cRefSlope) > cRefSlopeError){
        INFO(Form("Slope is not consistent with 0 within errors for %s", name.Data()));
    }
    return std::make_pair(fCRefFit, std::make_pair(ptMin, ptMax));
}


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
   
   // loop over shifted histogram
   for(int i = 1; i <= nbins; i++) {
     double binCenter = hShifted->GetBinCenter(i);
     double shiftedBinCenter = binCenter - x;
     double content = 0;
     double error = 0;
     int shiftedBin = -1;
     if(shiftedBinCenter >= min && shiftedBinCenter <= max){
      shiftedBin = h->FindBin(shiftedBinCenter);
      content = h->GetBinContent(shiftedBin);
      error = h->GetBinError(shiftedBin);
     }
     hShifted->SetBinContent(i, content);
     hShifted->SetBinError(i, error);
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
  
  // hRhoSignalPt->RebinX(4);
  // hRhoReferencePt->RebinX(4);
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
    hRatio->Divide(hRhoSignal,hShifted,1,1);
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

  TH1F* hBestShift = new TH1F("hBestShift", "best shift value", 1, -0.5, 0.5);
  hBestShift->SetBinContent(1, ptShifts[bestShiftIndex]);
  hBestShift->SetBinError(1, 0);
  hBestShift->Write("hBestShift");


}

void getRecoilJetPtDistributions(TDirectory *dSignal, TDirectory *dReference, TString label, GlobalOptions optns)
{
    DLJetCuts jetCuts(optns);
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

    // get THnSparseF for signal and reference (0 to 2pi)
    THnSparseF *hIsoGammaJetDeltaPhi2piJetPtGammaPtSignal = (THnSparseF*)dSignal->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
    hIsoGammaJetDeltaPhi2piJetPtGammaPtSignal->SetName("hIsoGammaJetDeltaPhi2piJetPtGammaPtSignal");
    THnSparseF *hIsoGammaJetDeltaPhi2piJetPtGammaPtReference = (THnSparseF*)dReference->Get("hIsoGammaJetDeltaPhi2piJetPtGammaPt");
    hIsoGammaJetDeltaPhi2piJetPtGammaPtReference->SetName("hIsoGammaJetDeltaPhi2piJetPtGammaPtReference");

    // get THnSparseF for signal and reference (0 to pi)
    THnSparseF *hIsoGammaJetDeltaPhiJetPtGammaPtSignal = (THnSparseF*)dSignal->Get("hIsoGammaJetDeltaPhiJetPtGammaPt");
    hIsoGammaJetDeltaPhiJetPtGammaPtSignal->SetName("hIsoGammaJetDeltaPhiJetPtGammaPtSignal");
    THnSparseF *hIsoGammaJetDeltaPhiJetPtGammaPtReference = (THnSparseF*)dReference->Get("hIsoGammaJetDeltaPhiJetPtGammaPt");
    hIsoGammaJetDeltaPhiJetPtGammaPtReference->SetName("hIsoGammaJetDeltaPhiJetPtGammaPtReference");

    // do projection onto Jet and DeltaPhi axis (0, pi)
    TH2F *hIsoGammaJetDeltaPhiJetPtSignal = (TH2F*)hIsoGammaJetDeltaPhiJetPtGammaPtSignal->Projection(1,0);
    hIsoGammaJetDeltaPhiJetPtSignal->SetName("hIsoGammaJetDeltaPhiJetPtSignal");
    TH2F *hIsoGammaJetDeltaPhiJetPtReference = (TH2F*)hIsoGammaJetDeltaPhiJetPtGammaPtReference->Projection(1,0);
    hIsoGammaJetDeltaPhiJetPtReference->SetName("hIsoGammaJetDeltaPhiJetPtReference");

    // do projection onto Jet and DeltaPhi axis (0, 2pi)
    TH2F *hIsoGammaJetDeltaPhi2piJetPtSignal = (TH2F*)hIsoGammaJetDeltaPhi2piJetPtGammaPtSignal->Projection(1,0);
    hIsoGammaJetDeltaPhi2piJetPtSignal->SetName("hIsoGammaJetDeltaPhi2piJetPtSignal");
    TH2F *hIsoGammaJetDeltaPhi2piJetPtReference = (TH2F*)hIsoGammaJetDeltaPhi2piJetPtGammaPtReference->Projection(1,0);
    hIsoGammaJetDeltaPhi2piJetPtReference->SetName("hIsoGammaJetDeltaPhi2piJetPtReference");

    // get mirrored histograms (0, 2pi)
    TH2F *hIsoGammaJetDeltaPhiJetPtSignalMirrored = mirrorTH2F(hIsoGammaJetDeltaPhiJetPtSignal);
    TH2F *hIsoGammaJetDeltaPhiJetPtReferenceMirrored = mirrorTH2F(hIsoGammaJetDeltaPhiJetPtReference);
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->SetName("hIsoGammaJetDeltaPhiJetPtSignalMirrored");
    hIsoGammaJetDeltaPhiJetPtReferenceMirrored->SetName("hIsoGammaJetDeltaPhiJetPtReferenceMirrored");




    // normalize by number of triggers 
    // rebin in pt by factor 4 in y axis
    hIsoGammaJetDeltaPhiJetPtSignal->RebinY(2); // 0 to pi
    hIsoGammaJetDeltaPhiJetPtReference->RebinY(2);
    hIsoGammaJetDeltaPhi2piJetPtSignal->RebinY(2); // 0 to 2pi
    hIsoGammaJetDeltaPhi2piJetPtReference->RebinY(2);
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->RebinY(2); // 0 to 2pi mirrored from 0 to pi
    hIsoGammaJetDeltaPhiJetPtReferenceMirrored->RebinY(2);

    double dPt = hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->GetBinWidth(1);
    double dPhi = hIsoGammaJetDeltaPhiJetPtSignal->GetXaxis()->GetBinWidth(1);
    hIsoGammaJetDeltaPhiJetPtSignal->Scale(1./(nTriggersSignal*dPt*dPhi));
    hIsoGammaJetDeltaPhiJetPtReference->Scale(1./(nTriggersReference*dPt*dPhi));
    dPt = hIsoGammaJetDeltaPhiJetPtSignalMirrored->GetYaxis()->GetBinWidth(1);
    dPhi = hIsoGammaJetDeltaPhiJetPtSignalMirrored->GetXaxis()->GetBinWidth(1);
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->Scale(1./(nTriggersSignal*dPt*dPhi));
    hIsoGammaJetDeltaPhiJetPtReferenceMirrored->Scale(1./(nTriggersReference*dPt*dPhi));
    dPt = hIsoGammaJetDeltaPhi2piJetPtSignal->GetYaxis()->GetBinWidth(1);
    dPhi = hIsoGammaJetDeltaPhi2piJetPtSignal->GetXaxis()->GetBinWidth(1);
    hIsoGammaJetDeltaPhi2piJetPtSignal->Scale(1./(nTriggersSignal*dPt*dPhi));
    hIsoGammaJetDeltaPhi2piJetPtReference->Scale(1./(nTriggersReference*dPt*dPhi));
    
    hIsoGammaJetDeltaPhiJetPtSignal->GetXaxis()->SetTitle("#Delta#phi");
    hIsoGammaJetDeltaPhiJetPtSignal->GetYaxis()->SetTitle("#it{p}_{T, ch jet}^{reco} (GeV/#it{c})");
    hIsoGammaJetDeltaPhiJetPtSignal->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi}");
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->GetXaxis()->SetTitle("#Delta#phi");
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->GetYaxis()->SetTitle("#it{p}_{T, ch jet}^{reco} (GeV/#it{c})");
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi}");
    hIsoGammaJetDeltaPhi2piJetPtSignal->GetXaxis()->SetTitle("#Delta#phi");
    hIsoGammaJetDeltaPhi2piJetPtSignal->GetYaxis()->SetTitle("#it{p}_{T, ch jet}^{reco} (GeV/#it{c})");
    hIsoGammaJetDeltaPhi2piJetPtSignal->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{3}N}{d#it{p}_{T, ch jet}^{reco} d#Delta#phi}");

        // get histograms shifted by -pi/2
    TH2F* hIsoGammaJetDeltaPhiJetPtSignalShifted = shiftDeltaPhiRange(hIsoGammaJetDeltaPhi2piJetPtSignal, -TMath::Pi()/2);
    hIsoGammaJetDeltaPhiJetPtSignalShifted->SetName("hIsoGammaJetDeltaPhiJetPtSignalShifted");
    TH2F* hIsoGammaJetDeltaPhiJetPtReferenceShifted = shiftDeltaPhiRange(hIsoGammaJetDeltaPhi2piJetPtReference, -TMath::Pi()/2);
    hIsoGammaJetDeltaPhiJetPtReferenceShifted->SetName("hIsoGammaJetDeltaPhiJetPtReferenceShifted");

    // write all histograms of interest to dRecoilJet
    dRecoilJet->cd();
    hIsoGammaJetDeltaPhiJetPtSignal->Write("hIsoGammaJetDeltaPhiJetPtSignal");
    hIsoGammaJetDeltaPhiJetPtReference->Write("hIsoGammaJetDeltaPhiJetPtReference");
    hIsoGammaJetDeltaPhi2piJetPtSignal->Write("hIsoGammaJetDeltaPhi2piJetPtSignal");
    hIsoGammaJetDeltaPhi2piJetPtReference->Write("hIsoGammaJetDeltaPhi2piJetPtReference");
    hIsoGammaJetDeltaPhiJetPtSignalMirrored->Write("hIsoGammaJetDeltaPhiJetPtSignalMirrored");
    hIsoGammaJetDeltaPhiJetPtReferenceMirrored->Write("hIsoGammaJetDeltaPhiJetPtReferenceMirrored");
    hIsoGammaJetDeltaPhiJetPtSignalShifted->Write("hIsoGammaJetDeltaPhiJetPtSignalShifted");
    hIsoGammaJetDeltaPhiJetPtReferenceShifted->Write("hIsoGammaJetDeltaPhiJetPtReferenceShifted");

    // arrays for integration in phi ranges
    double phiBinsMid[phiBins.size() - 1 ];
    double phiBinsUpError[phiBins.size() - 1 ];
    double phiBinsDownError[phiBins.size() - 1 ];
    for(int i=1; i<phiBins.size(); i++)
    {
      phiBinsMid[i-1] = phiBins[i].first + (phiBins[i].second - phiBins[i].first)/2;
      phiBinsUpError[i-1] = TMath::Abs(phiBins[i].first - phiBins[i].second)/2;
      phiBinsDownError[i-1] = TMath::Abs(phiBins[i].second - phiBins[i].first)/2;
    }
    double ptIntegralPhi[recoilJetPtBinsIntegration.size()][phiBins.size() - 1 ];
    double ptIntegralErrorPhi[recoilJetPtBinsIntegration.size()][phiBins.size() - 1 ];

    double cRefValues[phiBins.size() - 1];
    double cRefValuesError[phiBins.size() - 1];
    double cRefSlopeValues[phiBins.size() - 1];
    double cRefSlopeValuesError[phiBins.size() - 1];

    // do projections in phi bins
    int iPhiBin = 0;
    for(auto &phiBin : phiBins)
    {
      // now use histograms from 0 to pi
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
      auto fitResult = findOptimalCRefFit(hRecoilJetPtRatio, Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second));
      TF1* fCRefFit = new TF1(Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second), "pol0", fitResult.second.first, fitResult.second.second);
      hRecoilJetPtRatio->Fit(fCRefFit, "R0");
      TH1F* hCRefFitConfidence = new TH1F(Form("hCRefFitConfidence_%.2f_%.2f", phiBin.first, phiBin.second), "Confidence interval of fCRefFit", 200, fitResult.second.first, fitResult.second.second);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCRefFitConfidence);
      // determing fCRefFitPol1
      TF1* fCRefFitPol1 = new TF1(Form("fCRefFitPol1_%.2f_%.2f", phiBin.first, phiBin.second), "[0]*x + [1]", fitResult.second.first, fitResult.second.second);
      hRecoilJetPtRatio->Fit(fCRefFitPol1, "R0");
      TH1F* hCRefFitPol1Confidence = new TH1F(Form("hCRefFitPol1Confidence_%.2f_%.2f", phiBin.first, phiBin.second), "Confidence interval of fCRefFitPol1", 200, fitResult.second.first, fitResult.second.second);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hCRefFitPol1Confidence);
      double cRef = fCRefFit->GetParameter(0);
      double cRefSlope = fCRefFitPol1->GetParameter(0);
      double cRefSlopeError = fCRefFitPol1->GetParError(0);
      INFO(Form("Determined scaling factor for reference in range %.2f_%.2f: %.3f", fitResult.second.first, fitResult.second.second, cRef));
      
      TH1F* hRecoilJetReferenceScaled = (TH1F*)hRecoilJetPtReference->Clone(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetReferenceScaled->Scale(1./cRef);
      hRecoilJetReferenceScaled->SetName(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      
      TH1F* hRecoilJetPtSignal_Subtracted = (TH1F*)hRecoilJetPtSignal->Clone(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted->Add(hRecoilJetReferenceScaled, -1);
      hRecoilJetPtSignal_Subtracted->SetName(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));

      // do a rebinned distribution
      TH1F* hRecoilJetReferenceScaled_Rebinned = (TH1F*) hRecoilJetReferenceScaled->Clone(Form("hRecoilJetReferenceScaled_Rebinned_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetReferenceScaled_Rebinned->Rebin(4);
      hRecoilJetReferenceScaled_Rebinned->SetName(Form("hRecoilJetReferenceScaled_Rebinned_%.2f_%.2f", phiBin.first, phiBin.second));
      TH1F* hRecoilJetPtSignal_Subtracted_Rebinned = (TH1F*)hRecoilJetPtSignal->Clone(Form("hRecoilJetPtSignal_Subtracted_Rebinned_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted_Rebinned->Rebin(4);
      hRecoilJetPtSignal_Subtracted_Rebinned->Add(hRecoilJetReferenceScaled_Rebinned, -1);

      for (int i=0; i<recoilJetPtBinsIntegration.size(); i++){
        // check if phi min and phi max are the first of phiBins. Then we should continue because this is always the back-to-back integrated bin
        if(iPhiBin == 0) continue;

        double ptIntegral, ptIntegralError;
        int binMin = hRecoilJetPtSignal_Subtracted->FindBin(recoilJetPtBinsIntegration[i].first);
        int binMax = hRecoilJetPtSignal_Subtracted->FindBin(recoilJetPtBinsIntegration[i].second);
        ptIntegral = hRecoilJetPtSignal_Subtracted->IntegralAndError(binMin, binMax, ptIntegralError);
        // phi min and max
        ptIntegralPhi[i][iPhiBin-1] = ptIntegral;
        ptIntegralErrorPhi[i][iPhiBin-1] = ptIntegralError;
      }

      // store cRef values also in vector
      if(iPhiBin > 0){
      cRefValues[iPhiBin-1] = cRef;
      cRefValuesError[iPhiBin-1] = fCRefFit->GetParError(0);
      }

      // write all histograms of interest to dRecoilJet
      dRecoilJet->cd();
      hRecoilJetPtSignal->Write(Form("hRecoilJetPtSignal_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtReference->Write(Form("hRecoilJetPtReference_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtRatio->Write(Form("hRecoilJetPtSignalReferenceRatio_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetReferenceScaled->Write(Form("hRecoilJetReferenceScaled_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted->Write(Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second));
      hRecoilJetPtSignal_Subtracted_Rebinned->Write(Form("hRecoilJetPtSignal_Subtracted_Rebinned_%.2f_%.2f", phiBin.first, phiBin.second));
      fCRefFit->Write(Form("fCRefFit_%.2f_%.2f", phiBin.first, phiBin.second));
      fCRefFitPol1->Write(Form("fCRefFitPol1_%.2f_%.2f", phiBin.first, phiBin.second));
      hCRefFitConfidence->Write(Form("hCRefFitConfidence_%.2f_%.2f", phiBin.first, phiBin.second));
      hCRefFitPol1Confidence->Write(Form("hCRefFitPol1Confidence_%.2f_%.2f", phiBin.first, phiBin.second));


      iPhiBin++;
    }

    // write recoil jet signal subtracted phi graphs to dRecoilJet
    dRecoilJet->cd();
    
    // create TGraphAsymmErrors for each recoil jetpt integrations
    TGraphAsymmErrors* gRecoilJetPtIntegralPhi[recoilJetPtBinsIntegration.size()];
    for(int i=0; i<recoilJetPtBinsIntegration.size(); i++){
      gRecoilJetPtIntegralPhi[i] = new TGraphAsymmErrors(phiBins.size() - 1, phiBinsMid, ptIntegralPhi[i], phiBinsDownError, phiBinsUpError, ptIntegralErrorPhi[i], ptIntegralErrorPhi[i]);
      gRecoilJetPtIntegralPhi[i]->SetName(Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBinsIntegration[i].first, recoilJetPtBinsIntegration[i].second));
      gRecoilJetPtIntegralPhi[i]->Write(Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBinsIntegration[i].first, recoilJetPtBinsIntegration[i].second));
    }

    // create TGraphAsymmErrors for each cRef values
    TGraphAsymmErrors* gCRefValues = new TGraphAsymmErrors(phiBins.size() - 1, phiBinsMid, cRefValues, phiBinsDownError, phiBinsUpError, cRefValuesError, cRefValuesError);
    gCRefValues->SetName("gCRefValuesVsPhi");
    gCRefValues->Write("gCRefValuesVsPhi");

    // create TGraphAsymmErrors for each cRefSlope values
    TGraphAsymmErrors* gCRefSlopeValues = new TGraphAsymmErrors(phiBins.size() - 1, phiBinsMid, cRefSlopeValues, phiBinsDownError, phiBinsUpError, cRefSlopeValuesError, cRefSlopeValuesError);
    gCRefSlopeValues->SetName("gCRefSlopeValuesVsPhi");
    gCRefSlopeValues->Write("gCRefSlopeValuesVsPhi");

    // Study of correlations with rapidity gaps
    THnSparseF *hIsoGammaJetDeltaPhi2piDeltaEtaJetPt = (THnSparseF*)dSignal->Get("hIsoGammaJetDeltaPhi2piDeltaEtaJetPt");
    // axis 0: DeltaPhi axis 1: DeltaEta axis 2: jet pt
    dPt = hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(2)->GetBinWidth(1);
    dPhi = hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(0)->GetBinWidth(1);
    double dEta = hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(1)->GetBinWidth(1);

    // select events that are back to back in phi (which is phibin 0)
    hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(0)->SetRangeUser(phiBins[0].first, phiBins[0].second);

    // project onto dEta and jet pt
    TH2F* hIsoGammaJetDeltaEtaJetPt = (TH2F*)hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->Projection(2,1);
    hIsoGammaJetDeltaEtaJetPt->SetName("hIsoGammaJetDeltaEtaJetPt");
    hIsoGammaJetDeltaEtaJetPt->GetXaxis()->SetTitle("#Delta#eta");
    hIsoGammaJetDeltaEtaJetPt->GetYaxis()->SetTitle("#it{p}_{T}^{jet} (GeV/#it{c})");
    hIsoGammaJetDeltaEtaJetPt->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{3}N}{d#Delta#eta d#it{p}_{T}^{jet} d#Delta#phi}");
    hIsoGammaJetDeltaEtaJetPt->Scale(1./(nTriggersSignal*dPt*dPhi*dEta));

    // loop over dEta bins and do projections
    dRecoilJet->cd();
    hIsoGammaJetDeltaEtaJetPt->Write("hIsoGammaJetDeltaEtaJetPt");

    TH1F* hRecoilJetPtDEta[dEtaBins.size()-1];
    TH1F* hRecoilJetPtDEtaRatio[dEtaBins.size()-1];
    TH1F* hRecoilJetPtDEtaScaled[dEtaBins.size()-1];
    TH1F* hRecoilJetPtDEtaSubtracted[dEtaBins.size()-1];
    TF1* fRecoilJetPtDEtaRatioFit[dEtaBins.size()-1];
    int i = 0;
     for (auto dEtaBin : dEtaBins) {
      // select events in this dEta bin
      float dEtaMax = jetCuts.JetEtaPhiMinMax[0][1] + 0.67;
      if(dEtaBin.second > dEtaMax)
      {
        dEtaBin.second = dEtaMax;
      }

      hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(1)->SetRangeUser(dEtaBin.first, dEtaBin.second);
      double dEtaWidth = dEtaBin.second - dEtaBin.first;
      // project onto jet pt
      hRecoilJetPtDEta[i] = (TH1F*)hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->Projection(2);
      hRecoilJetPtDEta[i]->SetName(Form("hRecoilJetPtDEta_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      hRecoilJetPtDEta[i]->GetXaxis()->SetTitle("#it{p}_{T}^{jet} (GeV/#it{c})");
      hRecoilJetPtDEta[i]->GetYaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{2}N}{d#it{p}_{T}^{jet} d#Delta#phi d#Delta#eta}");
      hRecoilJetPtDEta[i]->Rebin(4);
      dPt = hRecoilJetPtDEta[i]->GetXaxis()->GetBinWidth(1);
      hRecoilJetPtDEta[i]->Scale(1./(nTriggersSignal*dPt*dPhi*dEtaWidth));

      // write histogram to directory
      dRecoilJet->cd();
      hRecoilJetPtDEta[i]->Write(Form("hRecoilJetPtDEta_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      i++;
     }

    for (int i = 0; i < dEtaBins.size()-1; i++) {
      auto dEtaBin = dEtaBins[i];
      // select events in this dEta bin
      float dEtaMax = jetCuts.JetEtaPhiMinMax[0][1] + 0.67;
      if(dEtaBin.second > dEtaMax)
      {
        dEtaBin.second = dEtaMax;
      }

      // reset axis range
      hIsoGammaJetDeltaPhi2piDeltaEtaJetPt->GetAxis(1)->UnZoom();

      hRecoilJetPtDEtaRatio[i] = (TH1F*)hRecoilJetPtDEta[i]->Clone(Form("hRecoilJetPtDEtaRatio_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      hRecoilJetPtDEtaRatio[i]->Divide(hRecoilJetPtDEta[i], hRecoilJetPtDEta[dEtaBins.size()-1]);
      hRecoilJetPtDEtaRatio[i]->Write(Form("hRecoilJetPtDEtaRatio_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
       
      // fit the ratio with a pol0
      fRecoilJetPtDEtaRatioFit[i] = new TF1(Form("fRecoilJetPtDEtaRatioFit_%.2f_%.2f", dEtaBin.first, dEtaBin.second), "pol0", 0, 100);
      hRecoilJetPtDEtaRatio[i]->Fit(fRecoilJetPtDEtaRatioFit[i], "R0");
      fRecoilJetPtDEtaRatioFit[i]->Write(Form("fRecoilJetPtDEtaRatioFit_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      
      // scale the distributions with 1/fRecoilJetPtDEtaRatioFit->GetParameter(0)
      hRecoilJetPtDEtaScaled[i] = (TH1F*)hRecoilJetPtDEta[i]->Clone(Form("hRecoilJetPtDEtaScaled_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      hRecoilJetPtDEtaScaled[i]->Scale(1./fRecoilJetPtDEtaRatioFit[i]->GetParameter(0));
      hRecoilJetPtDEtaScaled[i]->Write(Form("hRecoilJetPtDEtaScaled_%.2f_%.2f", dEtaBin.first, dEtaBin.second));

      hRecoilJetPtDEtaSubtracted[i] = (TH1F*)hRecoilJetPtDEtaScaled[i]->Clone(Form("hRecoilJetPtDEtaSubtracted_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
      hRecoilJetPtDEtaSubtracted[i]->Add(hRecoilJetPtDEta[dEtaBins.size()-1], -1);
      hRecoilJetPtDEtaSubtracted[i]->Write(Form("hRecoilJetPtDEtaSubtracted_%.2f_%.2f", dEtaBin.first, dEtaBin.second));
    }



}

std::pair<double, double> calculateABCDPurity(TH2F* hM02vsIsoPt, double ptMin, double ptMax)
{
    hM02vsIsoPt->Sumw2();


    // find bins for region A
    INFO(Form("Region A: %.2f_%.2f, %.2f_%.2f", regionA.minM02, regionA.maxM02, regionA.minIso, regionA.maxIso));
    std::tuple<int, int, int, int> binsA = {hM02vsIsoPt->GetXaxis()->FindBin(regionA.minM02), hM02vsIsoPt->GetXaxis()->FindBin(regionA.maxM02), hM02vsIsoPt->GetYaxis()->FindBin(regionA.minIso), hM02vsIsoPt->GetYaxis()->FindBin(regionA.maxIso)};
  
    // Integrate over region B
    std::tuple<int, int, int, int> binsB = {hM02vsIsoPt->GetXaxis()->FindBin(regionB.minM02), hM02vsIsoPt->GetXaxis()->FindBin(regionB.maxM02), hM02vsIsoPt->GetYaxis()->FindBin(regionB.minIso), hM02vsIsoPt->GetYaxis()->FindBin(regionB.maxIso)};
    INFO(Form("Region B: %.2f_%.2f, %.2f_%.2f", regionB.minM02, regionB.maxM02, regionB.minIso, regionB.maxIso));
    // Integrate over region C
    std::tuple<int, int, int, int> binsC = {hM02vsIsoPt->GetXaxis()->FindBin(regionC.minM02), hM02vsIsoPt->GetXaxis()->FindBin(regionC.maxM02), hM02vsIsoPt->GetYaxis()->FindBin(regionC.minIso), hM02vsIsoPt->GetYaxis()->FindBin(regionC.maxIso)};
    INFO(Form("Region C: %.2f_%.2f, %.2f_%.2f", regionC.minM02, regionC.maxM02, regionC.minIso, regionC.maxIso));
    // Integrate over region D
    std::tuple<int, int, int, int> binsD = {hM02vsIsoPt->GetXaxis()->FindBin(regionD.minM02), hM02vsIsoPt->GetXaxis()->FindBin(regionD.maxM02), hM02vsIsoPt->GetYaxis()->FindBin(regionD.minIso), hM02vsIsoPt->GetYaxis()->FindBin(regionD.maxIso)};
    INFO(Form("Region D: %.2f_%.2f, %.2f_%.2f", regionD.minM02, regionD.maxM02, regionD.minIso, regionD.maxIso));


    double nAErr;
    double nA = hM02vsIsoPt->IntegralAndError(std::get<0>(binsA), std::get<1>(binsA), std::get<2>(binsA), std::get<3>(binsA), nAErr);
    double nBErr;
    double nB = hM02vsIsoPt->IntegralAndError(std::get<0>(binsB), std::get<1>(binsB), std::get<2>(binsB), std::get<3>(binsB), nBErr);
    double nCErr;
    double nC = hM02vsIsoPt->IntegralAndError(std::get<0>(binsC), std::get<1>(binsC), std::get<2>(binsC), std::get<3>(binsC), nCErr);
    double nDErr;
    double nD = hM02vsIsoPt->IntegralAndError(std::get<0>(binsD), std::get<1>(binsD), std::get<2>(binsD), std::get<3>(binsD), nDErr);

    // calculate purity
    double purity = 0;
    double purityErr = 0;
    if((nA!=0) && (nB!=0) && (nD!=0)){
      purity = 1 - ((nC/nA)/(nD/nB));
      purityErr = TMath::Sqrt(((nA * nA * nB * nB * nC * nC * nDErr * nDErr) + (nA * nA * nB * nB * nD * nD * nCErr * nCErr) + (nA * nA * nC * nC * nD * nD * nBErr * nBErr) + (nAErr * nAErr * nB * nB * nC * nC * nD * nD)) / (pow(nA, 4) * pow(nD, 4)));
    }


    INFO(Form("Purity: %.3f #pm %.3f", purity, purityErr));
    INFO(Form("Entries of region histograms:"));
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
    std::pair<double, double> puritySignal = calculateABCDPurity(hM02vsIsoSignal,exclTrigSelection.getPtMinSignal(), exclTrigSelection.getPtMaxSignal());
    std::pair<double, double> purityReference = calculateABCDPurity(hM02vsIsoReference,exclTrigSelection.getPtMinReference(), exclTrigSelection.getPtMaxReference());

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

  getRecoilJetPtDistributions(dTriggerGammaJetCorrelationsSignal, dTriggerGammaJetCorrelationsReference, "Photons", optns);
  getRecoilJetPtDistributions(dTriggerMergedPi0JetCorrelationsSignal, dTriggerMergedPi0JetCorrelationsReference, "MergedPi0s", optns);

  // purity calculation needs to be called with none isolated photons
  processPurity(dGammas,"Photons", optns);


}