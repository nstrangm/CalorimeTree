#include "TMath.h"
#include "Cuts.h"
#include "Utilities.h"

enum NLFunction : int { kTestBeamShaper,
                        kTestBeamShaperWoScale,
                        kTestBeamFinalMC };

enum FinetuningFunction : int { kTestBeamDefaultMCRun1,
                                kTestBeamFinalMCRun1,
                                kTestBeamFinalMCRun2,
                                kTestBeamFinalMCRun2Sys1,
                                kTestBeamFinalMCRun2Sys2,
                                kTestBeamFinalMCRun213TeV,
                                kTestBeamFinalMCRun2Variation1,
                                kRun1Variation1,
                                kRun1Variation2 };

double applyNonLin(double clusE, int NLFunction)
{
  double fNonLinearityParams[5];
  double energy = clusE;
  switch (NLFunction) {
    case kTestBeamShaper: {
      fNonLinearityParams[0] = 1.91897;
      fNonLinearityParams[1] = 0.0264988;
      fNonLinearityParams[2] = 0.965663;
      fNonLinearityParams[3] = -187.501;
      fNonLinearityParams[4] = 2762.51;
      energy /= (1.0505 * (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Log(energy)) / (1 + (fNonLinearityParams[2] * TMath::Exp((energy - fNonLinearityParams[3]) / fNonLinearityParams[4]))));
      break;
    }
    case kTestBeamShaperWoScale: {
      fNonLinearityParams[0] = 1.91897;
      fNonLinearityParams[1] = 0.0264988;
      fNonLinearityParams[2] = 0.965663;
      fNonLinearityParams[3] = -187.501;
      fNonLinearityParams[4] = 2762.51;
      energy /= (1.0 * (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Log(energy)) / (1 + (fNonLinearityParams[2] * TMath::Exp((energy - fNonLinearityParams[3]) / fNonLinearityParams[4]))));
      break;
    }
    case kTestBeamFinalMC: {
      fNonLinearityParams[0] = 1.09357;
      fNonLinearityParams[1] = 0.0192266;
      fNonLinearityParams[2] = 0.291993;
      fNonLinearityParams[3] = 370.927;
      fNonLinearityParams[4] = 694.656;
      energy /= (1.00 * (fNonLinearityParams[0] + fNonLinearityParams[1] * TMath::Log(energy)) / (1 + (fNonLinearityParams[2] * TMath::Exp((energy - fNonLinearityParams[3]) / fNonLinearityParams[4]))));
      break;
    }
  }
  return energy;
}

// ###################################################################################
// ############################ Fine tuning ##########################################
// ###################################################################################

double applyFinetuning(double clusE, int NLFunction)
{
  double newEnergy = clusE;
  float fNLAfterburnerPara[11];
  Double_t p0, p1, p2, p3, p4, p5;
  switch (NLFunction) {
    case kTestBeamDefaultMCRun1: {
      // cout << "Applying kTestBeamDefaultMCRun1" << endl;
      fNLAfterburnerPara[0] = 1.00483;
      fNLAfterburnerPara[1] = -3.88706;
      fNLAfterburnerPara[2] = -0.141753;
      // Iteration-2 parameters
      fNLAfterburnerPara[3] = 1.;
      fNLAfterburnerPara[4] = 0.99; // factor to be applied on 1 cell clusters only
      fNLAfterburnerPara[5] = 0;
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      newEnergy /= fNLAfterburnerPara[3];

      break;
    }
    case kTestBeamFinalMCRun1: {
      // cout << "Applying kTestBeamFinalMCRun1" << endl;
      fNLAfterburnerPara[0] = 0.987912;
      fNLAfterburnerPara[1] = -2.94105;
      fNLAfterburnerPara[2] = -0.273207;
      // Iteration-2 parameters
      fNLAfterburnerPara[3] = 1.0125; // Run1 additional correction factor of 1.25%
      fNLAfterburnerPara[4] = 0.99;   // factor to be applied on 1 cell clusters only
      fNLAfterburnerPara[5] = 0;
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      newEnergy /= fNLAfterburnerPara[3];

      break;
    }
    case kTestBeamFinalMCRun2: {
      // There are no extracted parameters yet
      //   Iteration-1 paramters
      fNLAfterburnerPara[0] = 0.979235;
      fNLAfterburnerPara[1] = -3.17131;
      fNLAfterburnerPara[2] = -0.464198;

      //  Iteration-2 paramters
      fNLAfterburnerPara[3] = 1.0363369;
      fNLAfterburnerPara[4] = 0.5659247074;
      fNLAfterburnerPara[5] = -2.7818482972;
      fNLAfterburnerPara[6] = 1.0437012864;
      fNLAfterburnerPara[7] = 0.3620283273;
      fNLAfterburnerPara[8] = -2.8321172480;

      // interpolation between PCM-EMC and EMC-EMC
      fNLAfterburnerPara[9] = 1.0025; // pcm-emc would be 0

      // additional fine tuning for 1 cell clusters
      fNLAfterburnerPara[10] = 0.99;

      // iteration 1 FT (direct fit with exponential on ratio)
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      // iteration 2 FT (indirect fit with exponentials on mass pos.)
      if ((fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7] * newEnergy + fNLAfterburnerPara[8])) != 0) {
        newEnergy /= ((fNLAfterburnerPara[3] - TMath::Exp(-fNLAfterburnerPara[4] * newEnergy + fNLAfterburnerPara[5])) / (fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7] * newEnergy + fNLAfterburnerPara[8])));
      }
      // interpolation between PCMEMC and EMC (0.25%)
      newEnergy /= fNLAfterburnerPara[9];
      break;
    }
    case kTestBeamFinalMCRun2Sys1: {
      // There are no extracted parameters yet
      //   Iteration-1 paramters
      fNLAfterburnerPara[0] = 0.979235;
      fNLAfterburnerPara[1] = -3.17131;
      fNLAfterburnerPara[2] = -0.464198;

      //  Iteration-2 paramters
      fNLAfterburnerPara[3] = 1.0363369;
      fNLAfterburnerPara[4] = 0.5659247074;
      fNLAfterburnerPara[5] = -2.7818482972;
      fNLAfterburnerPara[6] = 1.0437012864;
      fNLAfterburnerPara[7] = 0.3620283273;
      fNLAfterburnerPara[8] = -2.8321172480;

      // for PCM-EMC
      fNLAfterburnerPara[9] = 1.; // pcm-emc would be 0

      // additional fine tuning for 1 cell clusters
      fNLAfterburnerPara[10] = 0.99;

      // iteration 1 FT (direct fit with exponential on ratio)
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      // iteration 2 FT (indirect fit with exponentials on mass pos.)
      //  if((fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7]*newEnergy+fNLAfterburnerPara[8]) ) != 0){
      //    newEnergy /= ( (fNLAfterburnerPara[3] - TMath::Exp(-fNLAfterburnerPara[4]*newEnergy+fNLAfterburnerPara[5]) )/(fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7]*newEnergy+fNLAfterburnerPara[8]) ) );
      //  }
      // interpolation between PCMEMC and EMC (0.25%)
      newEnergy /= fNLAfterburnerPara[9];
      break;
    }
    // energy /= FunctionNL_kSDM(energy, 0.984314, -3.30941, -0.399441);
    //         if(cluster->GetNCells() == 1){ // different fine tuning for 1 cell clusters
    //           energy /= 0.99;
    //         }
    case kTestBeamFinalMCRun2Sys2: {
      // There are no extracted parameters yet
      //   Iteration-1 paramters
      fNLAfterburnerPara[0] = 0.984314;
      fNLAfterburnerPara[1] = -3.30941;
      fNLAfterburnerPara[2] = -0.399441;

      //  Iteration-2 paramters
      fNLAfterburnerPara[3] = 1.0363369;
      fNLAfterburnerPara[4] = 0.5659247074;
      fNLAfterburnerPara[5] = -2.7818482972;
      fNLAfterburnerPara[6] = 1.0437012864;
      fNLAfterburnerPara[7] = 0.3620283273;
      fNLAfterburnerPara[8] = -2.8321172480;

      // for PCM-EMC
      fNLAfterburnerPara[9] = 1.; // pcm-emc would be 0

      // additional fine tuning for 1 cell clusters
      fNLAfterburnerPara[10] = 0.99;

      // iteration 1 FT (direct fit with exponential on ratio)
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      // iteration 2 FT (indirect fit with exponentials on mass pos.)
      //  if((fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7]*newEnergy+fNLAfterburnerPara[8]) ) != 0){
      //    newEnergy /= ( (fNLAfterburnerPara[3] - TMath::Exp(-fNLAfterburnerPara[4]*newEnergy+fNLAfterburnerPara[5]) )/(fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7]*newEnergy+fNLAfterburnerPara[8]) ) );
      //  }
      // interpolation between PCMEMC and EMC (0.25%)
      newEnergy /= fNLAfterburnerPara[9];
      break;
    }
    case kTestBeamFinalMCRun2Variation1: {
      // There are no extracted parameters yet
      //   Iteration-1 paramters
      fNLAfterburnerPara[0] = 0.99759;
      fNLAfterburnerPara[1] = -3.21271;
      fNLAfterburnerPara[2] = -0.363656;

      //  Iteration-2 paramters
      fNLAfterburnerPara[3] = 1.0363369;
      fNLAfterburnerPara[4] = 0.5659247074;
      fNLAfterburnerPara[5] = -2.7818482972;
      fNLAfterburnerPara[6] = 1.0437012864;
      fNLAfterburnerPara[7] = 0.3620283273;
      fNLAfterburnerPara[8] = -2.8321172480;

      // interpolation between PCM-EMC and EMC-EMC
      fNLAfterburnerPara[9] = 1.0025;

      // additional fine tuning for 1 cell clusters
      fNLAfterburnerPara[10] = 0.99;

      // iteration 1 FT (direct fit with exponential on ratio)
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      // iteration 2 FT (indirect fit with exponentials on mass pos.)
      if ((fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7] * newEnergy + fNLAfterburnerPara[8])) != 0) {
        newEnergy /= ((fNLAfterburnerPara[3] - TMath::Exp(-fNLAfterburnerPara[4] * newEnergy + fNLAfterburnerPara[5])) / (fNLAfterburnerPara[6] - TMath::Exp(-fNLAfterburnerPara[7] * newEnergy + fNLAfterburnerPara[8])));
      }
      // interpolation between PCMEMC and EMC (0.25%)
      newEnergy /= fNLAfterburnerPara[9];
      break;
    }
    case kTestBeamFinalMCRun213TeV: {
      // There are no extracted parameters yet
      // Iteration-1 paramters
      fNLAfterburnerPara[0] = 0.987912;
      fNLAfterburnerPara[1] = -2.94105;
      fNLAfterburnerPara[2] = -0.273207;
      // Iteration-2 parameters
      fNLAfterburnerPara[3] = 1.00349; // additional factor from average to improve agreement
      fNLAfterburnerPara[4] = 0.99;    // factor to be applied on 1 cell clusters only
      fNLAfterburnerPara[5] = 0;
      newEnergy /= fNLAfterburnerPara[0] + TMath::Exp(fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
      newEnergy /= fNLAfterburnerPara[3];

      break;
    }
    case kRun1Variation1: // has to be applied on top of existing NL for PCM-EMC
    {
      p0 = 0.997021;
      p1 = 0.000763945;
      p2 = -2.76814e-05;
      p3 = 4.87087e-07;
      p4 = -4.12892e-09;
      p5 = 1.3477e-11;
      newEnergy *= (p0 + p1 * newEnergy + p2 * pow(newEnergy, 2) + p3 * pow(newEnergy, 3) + p4 * pow(newEnergy, 4) + p5 * pow(newEnergy, 5));
      // newEnergy*= (0.997021 + 0.000763945 * newEnergy + -2.76814e-05 * pow(newEnergy,2) +4.87087e-07 * pow(newEnergy,3) +-4.12892e-09 * pow(newEnergy,4) +1.3477e-11 * pow(newEnergy,5)) ;
      break;
    }
    case kRun1Variation2: // has to be applied on top of existing NL for PCM-EMC
    {
      p0 = 1.00345;
      p1 = -0.00104611;
      p2 = 3.47909e-05;
      p3 = -5.57732e-07;
      p4 = 4.3178e-09;
      p5 = -1.29549e-11;
      newEnergy *= (p0 + p1 * newEnergy + p2 * pow(newEnergy, 2) + p3 * pow(newEnergy, 3) + p4 * pow(newEnergy, 4) + p5 * pow(newEnergy, 5));
      break;
    }
  }
  return newEnergy;
}

void applyNonLinAndFineTuningCorrection(std::vector<Cluster>& IsoGammas, IsoGammaCuts isoGammaCuts, GlobalOptions optns)
{
  for (int iCluster = 0; iCluster < (int)IsoGammas.size(); iCluster++) {
    float clusterECorrected = IsoGammas.at(iCluster).E;
    if (!optns.isMC) {
      clusterECorrected = applyNonLin(clusterECorrected, kTestBeamShaperWoScale);
    } else {
      clusterECorrected = applyNonLin(clusterECorrected, kTestBeamFinalMC);
      if (isoGammaCuts.NonLinMode == 0) {
        clusterECorrected = applyFinetuning(clusterECorrected, kTestBeamFinalMCRun2);
      } else if (isoGammaCuts.NonLinMode == 2) {
        clusterECorrected = applyFinetuning(clusterECorrected, kTestBeamFinalMCRun2Sys1);
      } else if (isoGammaCuts.NonLinMode == 3) {
        clusterECorrected = applyFinetuning(clusterECorrected, kTestBeamFinalMCRun2Sys2);
      } else {
        clusterECorrected = applyFinetuning(clusterECorrected, kTestBeamFinalMCRun2);
      }
    }
    float clusterECorrection = clusterECorrected / IsoGammas.at(iCluster).E;
    IsoGammas.at(iCluster).px *= clusterECorrection;
    IsoGammas.at(iCluster).py *= clusterECorrection;
    IsoGammas.at(iCluster).pz *= clusterECorrection;
    IsoGammas.at(iCluster).E *= clusterECorrection;
  }
}