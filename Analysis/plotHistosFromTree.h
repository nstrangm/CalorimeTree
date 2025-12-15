#ifndef _PLOTHISTOSFROMTREE_H_
#define _PLOTHISTOSFROMTREE_H_

std::vector<float> SlicesPt={3,4,5,6,7,9,11};

// to be used fo cluster QA ( or whatever)
  const Int_t nPtBins = 9;
  Double_t ptBins[nPtBins+1] = {10, 15, 20, 25 ,30, 40, 50, 60,70,80};

  const Int_t nPtBinsCoarse = 5;
  Double_t ptBinsCoarse[nPtBinsCoarse+1] = {10, 20, 30, 40, 50, 80};

// used for jets
  int const nJetPtBins = 9;
  Double_t jetPtBins[nJetPtBins+1] = {0,5,10,15,20,30,40,50,150,200};

  const Int_t nBinsOccupancy = 4;
  Double_t binsOccupancy[nBinsOccupancy+1] = {0,1000,2000,5000,10000};

  Int_t fancyColors[8];

  std::pair<double, double> jetExamplePtRange = {20, 30};

  const int nJetConstituentBins = 8;
  Double_t jetConstituentBins[nJetConstituentBins+1] = {0,5,10,15,20,30,40,50,100};
  TString legendHeader = "ALICE work-in-progress";
  TString legendHeaderJets = "ALICE work-in-progress";
#endif //_PLOTHISTOSFROMTREE_H_
