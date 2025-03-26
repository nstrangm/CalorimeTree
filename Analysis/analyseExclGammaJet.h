#ifndef _ANALYSEEXCLGAMMAJET_H_
#define _ANALYSEEXCLGAMMAJET_H_

#include "THnSparse.h"

double nTriggersSignalPhoton = 0;
double nTriggersReferencePhoton = 0;
double nTriggersSignalMerged = 0;
double nTriggersReferenceMerged = 0;
std::vector<float> SlicesPt={3,4,5,6,7,9,11};

// to be used fo cluster QA ( or whatever)
TFile* fOut = nullptr;

// phi bins for analysis. First bin always has to be back to back
// std::vector<std::pair<double, double>> phiBins = {{TMath::Pi()-0.6, TMath::Pi()+0.6}, {1.57,2.07}, {2.61,2.76},{3.02,3.14}};
// finer phi bins
std::vector<std::pair<double, double>> phiBins = {{TMath::Pi()-0.6, TMath::Pi()+0.6}, {1.57,2.07},{2.08,2.22},{2.23,2.37},{2.38,2.52},{2.53,2.67},{2.68,2.82},{2.83,2.97},{2.98,3.12}};

std::vector<std::pair<double, double>> dEtaBins = {{0.0, 0.5},{0.5,1}, {1, 1.3}};
std::pair<double, double> cRefFitRange = {-50, 0};

std::vector<std::pair<double, double>> recoilJetPtBinsIntegration = {{20,30},{30,40},{40,50}};

// cuts for ABCD purity calculation
struct Region
{
  double minM02;
  double maxM02;
  double minIso;
  double maxIso;
};

Region regionA = {0.1, 0.3, -1000, 1.5};
Region regionB = {0.4, 2.0, -1000, 1.5};
Region regionC = {0.4, 2.0, 4, 1000};
Region regionD = {0.1, 0.3, 4, 1000};

// Function to mirror a TH2 histogram around x = pi to extend range from to 0, 2pi. The input histogram is expected to range from 0 to pi.
TH2F* mirrorTH2F(TH2F* h)
{
  // Get the original histogram dimensions
  int nBinsX = h->GetNbinsX();
  int nBinsY = h->GetNbinsY();
  double xMin = h->GetXaxis()->GetXmin();
  double xMax = h->GetXaxis()->GetXmax();
  double yMin = h->GetYaxis()->GetXmin();
  double yMax = h->GetYaxis()->GetXmax();
  
  // Create a new histogram with range 0 to 2pi on x-axis
  TH2F* hMirrored = new TH2F(Form("%s_mirrored", h->GetName()),Form("%s_mirrored", h->GetTitle()),nBinsX*2, 0, 2*TMath::Pi(),nBinsY, yMin, yMax);
  
  // Copy axis titles
  hMirrored->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hMirrored->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  
  // Fill the new histogram
  for (int ix = 1; ix <= nBinsX; ix++) {
    for (int iy = 1; iy <= nBinsY; iy++) {
      double content = h->GetBinContent(ix, iy);
      double error = h->GetBinError(ix, iy);
      
      // Fill the original part (0 to pi)
      hMirrored->SetBinContent(ix, iy, content);
      hMirrored->SetBinError(ix, iy, error);
      
      // Fill the mirrored part (pi to 2pi)
      int mirroredBin = 2*nBinsX + 1 - ix;
      hMirrored->SetBinContent(mirroredBin, iy, content);
      hMirrored->SetBinError(mirroredBin, iy, error);
    }
  }
  return hMirrored;
}



#endif //_ANALYSEEXCLGAMMAJET_H_
