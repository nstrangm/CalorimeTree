#ifndef _UTILITIES_
#define _UTILITIES_

#include "PlottingClass.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <sys/stat.h>
#include <sys/types.h>

class GlobalOptions
{
public:
  bool isDebugRun = false; // Active if jobId = 0 to only run over one file for debugging
  // int jobId = 0; // Slurm jobId
  bool isMC = false;
  TString dataSet;
  TString dataSetLabel;
  TString trainConfig;
  TString cutString;
  TString analysisDirPath;
  int TreeFormat = 0;
  bool doQA = true;
  bool doIsoGamma = false;
  bool doJets = false;
  bool domPi0 = false;
  bool doGGPi0 = false;
  bool doCombineExclGammaJet = false;
  GlobalOptions(TString AnalysisDirectory, int jobId);
  // GlobalOptions(bool userWantsMC, bool userWantsQA, TString EventCutString, TString IsoGammaCutString, TString JetCutString, TString Pi0CutString);
  ~GlobalOptions() {};
};

// GlobalOptions::GlobalOptions(bool userWantsMC, bool userWantsQA, TString EventCutString, TString IsoGammaCutString, TString JetCutString, TString Pi0CutString)
GlobalOptions::GlobalOptions(TString AnalysisDirectory, int jobId)
{
  if (jobId < 0)
    FATAL("Negative jobId")
  if (!jobId)
    LOG("This is a debug run")

  if (jobId == 0)
  {
    isDebugRun = true;
    debugLevel = 3;
  }

  analysisDirPath = AnalysisDirectory;
  std::stringstream ss((std::string)AnalysisDirectory.Data());
  std::string word;
  std::string words[3];

  for (int i = 0; i < 3; i++)
  {
    std::getline(ss, word, '/');
    words[i] = word;
  }
  dataSet = (TString)words[0].c_str();
  trainConfig = (TString)words[1].c_str();
  cutString = (TString)words[2].c_str();

  INFO(Form("dataSet = %s, trainConfig = %s, cutString = %s", dataSet.Data(), trainConfig.Data(), cutString.Data()))

  YAML::Node config = YAML::LoadFile("RunConfig.yaml");

  doIsoGamma = config["doIsoGamma"].as<bool>();
  doJets = config["doJets"].as<bool>();
  domPi0 = config["domPi0"].as<bool>();
  doGGPi0 = config["doGGPi0"].as<bool>();
  doCombineExclGammaJet = config["doCombineExclGammaJet"].as<bool>();

  if (!config[(std::string)dataSet])
    FATAL(Form("Dataset %s not found in YAML file RunConfig.yaml", dataSet.Data()))

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

  INFO(Form("Analyzing %s%s%s%s with%s QA in %s", doIsoGamma ? "Isolated Gammas" : "", doJets ? " + Jets" : "", domPi0 ? " + mPi0" : "", doGGPi0 ? " + GGPi0" : "", doQA ? "" : "out", isMC ? "MC" : "data"))
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

void createDirectory(TString path)
{
  struct stat info;
  if (stat(path.Data(), &info) != 0)
  {
    if (mkdir(path.Data(), 0755) == -1) // Directory does not exist, create it
      std::cerr << "Error creating directory: " << path.Data() << std::endl;
  }
  return;
}

// define centrality enums
enum CentralityEnum{
  k0_10,
  k10_30,
  k30_50,
  k50_90,
  k0_90
};
enum JetRadiusEnum{
  kR02,
  kR03,
  kR04,
  kR05,
  kR06
};

class FilePath{
  public:
    TString filePath;
    FilePath(TString filePath);
    // centralit enum
    CentralityEnum centrality;
    JetRadiusEnum jetRadius;
};
// Extract the centrality and jet radius from the filePath. The string is expected to be in the format of "Run3-0-10/JetRadius_R02/ExclGammaJet.root"
FilePath::FilePath(TString filePath){
  this->filePath = filePath;
  // Extract centrality from the file path
  if (filePath.Contains("Run3-0-10")) {
    centrality = k0_10;
  } else if (filePath.Contains("Run3-10-30")) {
    centrality = k10_30;
  } else if (filePath.Contains("Run3-30-50")) {
    centrality = k30_50;
  } else if (filePath.Contains("Run3-50-90")) {
    centrality = k50_90;
  } else if (filePath.Contains("Run3-0-90")) {
    centrality = k0_90;
  } else {
    std::cerr << "Unknown centrality in file path: " << filePath.Data() << std::endl;
  }
  
  // Extract jet radius from the file path
  if (filePath.Contains("JetRadius_R02") || filePath.Contains("R02_")) {
    jetRadius = kR02;
  } else if (filePath.Contains("JetRadius_R03") || filePath.Contains("R03_")) {
    jetRadius = kR03;
  } else if (filePath.Contains("JetRadius_R04") || filePath.Contains("R04_")) {
    jetRadius = kR04;
  } else if (filePath.Contains("JetRadius_R05") || filePath.Contains("R05_")) {
    jetRadius = kR05;
  } else if (filePath.Contains("JetRadius_R06") || filePath.Contains("R06_")) {
    jetRadius = kR06;
  } else if (filePath.Contains("Standard")) {
    jetRadius = kR04; // Default radius for Standard configuration
  } else {
    std::cerr << "Unknown jet radius in file path: " << filePath.Data() << std::endl;
  }
  
}

class CombineExclGammaJetOptions{
  public:
    TString dataset;
    std::vector<TString> jetRadiusFolders;
    std::vector<GlobalOptions> jetRadiusOptions; // also store global options to retrieve used cuts
    std::vector<TString> centralityFolders;
    TString analysisDirPath; 
    CombineExclGammaJetOptions(TString AnalysisDirectory, TString configYaml);
    std::vector<FilePath> inputFiles;
    TString GetInputFilePath(CentralityEnum centrality, JetRadiusEnum jetRadius){
      for (auto inputFile : inputFiles){
        if (inputFile.centrality == centrality && inputFile.jetRadius == jetRadius){
          return inputFile.filePath;
        }
      }
      return "";
    }
    TFile* GetRootFile(CentralityEnum centrality, JetRadiusEnum jetRadius){
      return TFile::Open(GetInputFilePath(centrality, jetRadius).Data());
    }
    ~CombineExclGammaJetOptions() {};
    void printOptions();
};

CombineExclGammaJetOptions::CombineExclGammaJetOptions(TString AnalysisDirectory, TString configYaml){
  analysisDirPath = AnalysisDirectory;
  YAML::Node config = YAML::LoadFile(configYaml.Data());
  YAML::Node combineExclGammaJetConfig = config["combineExclGammaJet"];
  dataset = combineExclGammaJetConfig["dataset"].as<string>().c_str();
  // Load jetRadiusFolders from comma separated list in YAML
  std::string jetRadiusFoldersStr = combineExclGammaJetConfig["jetRadiusFolders"].as<string>();
  std::stringstream ss(jetRadiusFoldersStr);
  std::string folder;
  
  // Parse comma-separated values into vector
  while (std::getline(ss, folder, ',')) {
    jetRadiusFolders.push_back(TString(folder.c_str()));
  }
  
  // Load centralityFolders from comma separated list in YAML
  std::string centralityFoldersStr = combineExclGammaJetConfig["centralityFolders"].as<string>();
  std::stringstream ss2(centralityFoldersStr);
  
  // Parse comma-separated values into vector
  while (std::getline(ss2, folder, ',')) {
    centralityFolders.push_back(TString(folder.c_str()));
  }
  
  // loop over all jet radius folders and centrality folders and create a FilePath object for each
  for (auto jetRadiusFolder : jetRadiusFolders){
    // ugly workaround but only store for one centrality folder 
    // centrality does not matter for the global options
    // which we will only use in this case to get the trigger cuts
    TString jetRadiusFolderPath =  dataset + "/" + centralityFolders[0] + "/" + jetRadiusFolder;
    INFO(Form("Jet Radius Folder: %s", jetRadiusFolderPath.Data()));
    jetRadiusOptions.emplace_back(GlobalOptions(jetRadiusFolderPath,0));
    for (auto centralityFolder : centralityFolders){
      inputFiles.push_back(FilePath(dataset + "/" + centralityFolder + "/" + jetRadiusFolder + "/ExclGammaJet.root"));
    }
  }
}

void CombineExclGammaJetOptions::printOptions(){
  INFO(Form("Dataset: %s", dataset.Data()));
  INFO(Form("Analysis Directory: %s", analysisDirPath.Data()));
  INFO(Form("Jet Radius Folders:"));
  for (auto folder : jetRadiusFolders){
    INFO(Form("  %s", folder.Data()));
  }
  INFO(Form("Centrality Folders:"));
  for (auto folder : centralityFolders){
    INFO(Form("  %s", folder.Data()));
  }
}

#endif
