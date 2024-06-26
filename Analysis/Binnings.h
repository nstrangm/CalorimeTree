#ifndef _BINNINGS_H_
#define _BINNINGS_H_
//For Jets
const int nBins_JESR = 26;
const float binning_JESR[nBins_JESR + 1] = {0., 5., 6., 7., 8., 9., 10., 11., 12., 13.,
                                            14., 15., 17., 20., 25., 30., 35., 40., 45., 50.,
                                            60., 70., 85., 100., 120., 150., 200.};

//For IsoGammas
int npTbinEdges=5;
float pTbinEdges[5] = {10,15,20,40,80};
float pTbinsToplot[4] = {12.5,17.5,30,60};
//Set pTIsoCuts [min isolated, max isolated, min antiisolated, max antiisolated]
float pTIsoABCDCuts[4] = {0,1.5,4.0,30};
float M02ABCDCuts[4] = {0.1,0.3,0.4,2.0};

#endif //_BINNINGS_H_