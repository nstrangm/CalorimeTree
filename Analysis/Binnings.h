#ifndef _BINNINGS_H_
#define _BINNINGS_H_
//For Jets
const int nBins_JESR = 26;
const float binning_JESR[nBins_JESR + 1] = {0., 5., 6., 7., 8., 9., 10., 11., 12., 13.,
                                            14., 15., 17., 20., 25., 30., 35., 40., 45., 50.,
                                            60., 70., 85., 100., 120., 150., 200.};

//For IsoGammas
//For ABCD
int npTbinEdges=9;
float pTbinEdges[9] = {10,12.5,15,17.5,20,30,40,60,80};
float pTbinsToplot[8] = {11.25,13.75,16.25,18.75,25,35,50,70};
//For efficiencies, purity and acceptance
int npTbinEdges2=10;
Double_t pTbinEdges2[10] = {0,10,12.5,15,17.5,20,30,40,60,80};
Double_t pTbinsToplot2[8] = {11.25,13.75,16.25,18.75,25,35,50,70};
//Set pTIsoCuts [min isolated, max isolated, min antiisolated, max antiisolated]
float pTIsoABCDCuts[4] = {-50,1.5,4.0,100};
float M02ABCDCuts[4] = {0.1,0.3,0.4,2.0};

#endif //_BINNINGS_H_