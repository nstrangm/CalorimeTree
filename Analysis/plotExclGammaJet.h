#ifndef _PLOTEXCLGAMMAJET_H_
#define _PLOTEXCLGAMMAJET_H_
#include "THnSparse.h"


// information about the run configuration
std::pair<double, double> centralityRange;
double jetRadius = 0.4;
std::pair<double, double> signalTriggerPt;
std::pair<double, double> referenceTriggerPt;

TString outputDir;
TString suffix = "pdf";

const Int_t nPtBins = 9;
Double_t ptBins[nPtBins+1] = {10, 15, 20, 25 ,30, 40, 50, 60,70,80};

const Int_t nPtBinsCoarse = 5;
Double_t ptBinsCoarse[nPtBinsCoarse+1] = {10, 20, 30, 40, 50, 80};

// used for jets
int const nJetPtBins = 9;
Double_t jetPtBins[nJetPtBins+1] = {0,5,10,15,20,30,40,50,150,200};

const Int_t nBinsOccupancy = 4;
Double_t binsOccupancy[nBinsOccupancy+1] = {0,1000,2000,5000,10000};



#endif //_PLOTEXCLGAMMAJET_H_