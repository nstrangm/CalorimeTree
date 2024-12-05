#ifndef _PLOTEXCLGAMMAJET_H_
#define _PLOTEXCLGAMMAJET_H_


// information about the run configuration
std::pair<double, double> centralityRange;
double jetRadius = 0.4;
std::pair<double, double> signalTriggerPt;
std::pair<double, double> referenceTriggerPt;

TString outputDir;
TString suffix = "pdf";



#endif //_PLOTEXCLGAMMAJET_H_