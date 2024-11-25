#ifndef _PLOTHISTOSFROMTREE_H_
#define _PLOTHISTOSFROMTREE_H_

std::vector<float> SlicesPt={3,4,5,6,7,9,11};

// to be used fo cluster QA ( or whatever)
  const Int_t nPtBins = 9;
  Double_t ptBins[nPtBins+1] = {10, 15, 20, 25 ,30, 40, 50, 60,70,80};

  const Int_t nPtBinsCoarse = 5;
  Double_t ptBinsCoarse[nPtBinsCoarse+1] = {10, 20, 30, 40, 50, 80};

#endif //_PLOTHISTOSFROMTREE_H_