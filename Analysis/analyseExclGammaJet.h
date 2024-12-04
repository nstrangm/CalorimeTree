#ifndef _ANALYSEEXCLGAMMAJET_H_
#define _ANALYSEEXCLGAMMAJET_H_


double nTriggersSignalPhoton = 0;
double nTriggersReferencePhoton = 0;
double nTriggersSignalMerged = 0;
double nTriggersReferenceMerged = 0;
std::vector<float> SlicesPt={3,4,5,6,7,9,11};

// to be used fo cluster QA ( or whatever)
TFile* fOut = nullptr;

// phi bins for analysis. First bin always has to be back to back
std::vector<std::pair<double, double>> phiBins = {{TMath::Pi()-0.6, TMath::Pi()+0.6}, {1.57,2.07}, {2.61,2.76},{3.02,3.14}};

std::pair<double, double> cRefFitRange = {-50, 0};

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


#endif //_ANALYSEEXCLGAMMAJET_H_
