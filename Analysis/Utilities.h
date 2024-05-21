#ifndef _UTILITIES_
#define _UTILITIES_

#include "../Plotting/PlottingClass.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <yaml-cpp/yaml.h>

class GlobalOptions
{
public:
  bool isMC = false;
  TString parentDir;
  TString dataSet;
  TString dataSetLabel;
  TString trainConfig;
  TString cutString;
  int TreeFormat = 0;
  bool doQA = true;
  bool doIsoGamma = false;
  bool doJets = false;
  bool doPi0 = false;
  GlobalOptions(TString AnalysisDirectory);
  // GlobalOptions(bool userWantsMC, bool userWantsQA, TString EventCutString, TString IsoGammaCutString, TString JetCutString, TString Pi0CutString);
  ~GlobalOptions(){};
};

// GlobalOptions::GlobalOptions(bool userWantsMC, bool userWantsQA, TString EventCutString, TString IsoGammaCutString, TString JetCutString, TString Pi0CutString)
GlobalOptions::GlobalOptions(TString AnalysisDirectory)
{
  std::stringstream ss((std::string)AnalysisDirectory.Data());
  std::string word;
  std::string words[4];

  for (int i = 0; i < 4; i++)
  {
    std::getline(ss, word, '/');
    words[i] = word;
  }
  parentDir = (TString)words[0].c_str();
  dataSet = (TString)words[1].c_str();
  trainConfig = (TString)words[2].c_str();
  cutString = (TString)words[3].c_str();

  INFO(Form("parentDir = %s, dataSet = %s, trainConfig = %s, cutString = %s", parentDir.Data(), dataSet.Data(), trainConfig.Data(), cutString.Data()))

  const char *YAMLFilePath = Form("%s/RunConfig.yaml", parentDir.Data());
  YAML::Node config = YAML::LoadFile(YAMLFilePath);

  doIsoGamma = config["doIsoGamma"].as<bool>();
  doJets = config["doJets"].as<bool>();
  doPi0 = config["doPi0"].as<bool>();

  if (!config[(std::string)dataSet])
    FATAL(Form("Dataset %s not found in YAML file %s", dataSet.Data(), YAMLFilePath))

  YAML::Node dataSetConfig = config[(std::string)dataSet];
  dataSetLabel = dataSetConfig["label"].as<string>().c_str();
  isMC = dataSetConfig["isMC"].as<bool>();

  std::cout << R"(
  ________/\\\\\\\\\_________________/\\\\\\__________________________________________________________________________/\\\\\\\\\\\\\\\_____________________________________________
   _____/\\\////////_________________\////\\\_________________________________________________________________________\///////\\\/////______________________________________________
    ___/\\\/_____________________________\/\\\_________________________________/\\\__________________________________________\/\\\___________________________________________________
     __/\\\______________/\\\\\\\\\_______\/\\\________/\\\\\_____/\\/\\\\\\\__\///_____/\\\\\__/\\\\\_______/\\\\\\\\________\/\\\________/\\/\\\\\\\______/\\\\\\\\______/\\\\\\\\__
      _\/\\\_____________\////////\\\______\/\\\______/\\\///\\\__\/\\\/////\\\__/\\\__/\\\///\\\\\///\\\___/\\\/////\\\_______\/\\\_______\/\\\/////\\\___/\\\/////\\\___/\\\/////\\\_
       _\//\\\______________/\\\\\\\\\\_____\/\\\_____/\\\__\//\\\_\/\\\___\///__\/\\\_\/\\\_\//\\\__\/\\\__/\\\\\\\\\\\________\/\\\_______\/\\\___\///___/\\\\\\\\\\\___/\\\\\\\\\\\__
        __\///\\\___________/\\\/////\\\_____\/\\\____\//\\\__/\\\__\/\\\_________\/\\\_\/\\\__\/\\\__\/\\\_\//\\///////_________\/\\\_______\/\\\_________\//\\///////___\//\\///////___
         ____\////\\\\\\\\\_\//\\\\\\\\/\\__/\\\\\\\\\__\///\\\\\/___\/\\\_________\/\\\_\/\\\__\/\\\__\/\\\__\//\\\\\\\\\\_______\/\\\_______\/\\\__________\//\\\\\\\\\\__\//\\\\\\\\\\_
          _______\/////////___\////////\//__\/////////_____\/////_____\///__________\///__\///___\///___\///____\//////////________\///________\///____________\//////////____\//////////__
        )" << '\n';

  LOG(Form("Analyzing %s%s%s with%s QA in %s", doIsoGamma ? "Isolated Gammas" : "", doJets ? " + Jets" : "", doPi0 ? " + Pi0" : "", doQA ? "" : "out", isMC ? "MC" : "data"))
}

// Returns the fraction of a cone at Eta = cEta and radius r that lies within the acceptance of a detector covering +-maxEta
float CalculateIsoCorrectionFactor(float cEta, float maxEta, float r)
{
  // Distance of cluster eta to detector limit eta
  float distanceEta = TMath::Abs(maxEta) - TMath::Abs(cEta);

  // If cluster is further away from border than r, no need to carry on. Everything fits
  if (distanceEta > r)
    return 1.;

  // Calculate angle between eta acceptance limit and outermost cone eta
  float angleBeta = 2 * TMath::ACos(distanceEta / r);

  // Calculate sector of circle outside of acceptance
  float areaExcess = 0.5 * r * r * (angleBeta - TMath::Sin(angleBeta));

  // Calculate fraction of area covered
  float totalArea = TMath::Pi() * r * r;
  float coveredArea = totalArea - areaExcess;

  // Return the fraction of the cone within the acceptance
  return coveredArea / totalArea;
}

void setExpArray(Double_t array[], int nPoints, float min, float max)
{
  double x = TMath::Power(max, 1. / (1. - nPoints)) / TMath::Power(min, nPoints / (1. - nPoints));
  double y = TMath::Log10(min / x);
  for (int i = 0; i < nPoints; i++)
    array[i] = x * TMath::Power(10, y * i);
}

#endif