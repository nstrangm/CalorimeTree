//Load libraries
//#include "TH1.h"
//#include "TH2.h"
//#include "TF1.h"
//#include "TF2.h"
//#include "TLine.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TGraphAsymmErrors.h"
//#include "TCanvas.h"
//#include "TLegend.h"
//#include "TObjString.h"
//#include "TLatex.h"
//#include "TString.h"
//#include "TStyle.h"
//#include "TColor.h"
//#include "TFile.h"
//#include "TRandom.h"
//#include "TTree.h"
//#include "TEllipse.h"
//#include "TCurlyLine.h"
//#include "TMath.h"
//#include <iostream>
//#include <string>
//#include <vector>
//#include <time.h>
//#include "TDatabasePDG.h"
//#include "../Analysis/Logging.h"
#include "TFile.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <filesystem>
#include "../Analysis/Utilities.h"
#include "../Analysis/PlottingClass.h"
#include <yaml-cpp/yaml.h>

const char* suffix = "png";

void createDirectory(const char* path)
{
  struct stat info;
  if (stat(path, &info) != 0) {
    if (mkdir(path, 0755) == -1) // Directory does not exist, create it
      std::cerr << "Error creating directory: " << path << std::endl;
  }
  return;
}

void CompareCuts(const char* path){
    //Read config file
    YAML::Node config = YAML::LoadFile("CompareCutsSetup.yaml");

    //Prepare Output directory
    createDirectory(path);

    //Get number of files
    int nfiles = config["Nfiles"].as<int>();
    //Get number of plots
    int nplots = config["Nplots"].as<int>();
    //Get number of groups
    int ngroups = config["Ngroups"].as<int>();

    //Get plots
    YAML::Node toplot = config["ToPlot"];
    YAML::Emitter outplot1;
    outplot1 << toplot[0];
    std::string plotstring1 = outplot1.c_str();
    

    //Get filepaths
    YAML::Node filepaths = config["FilePaths"];

    //Get groups
    YAML::Node groups = config["Groups"];

    //Get identifiers
    YAML::Node fileIdent = config["Identifier"];

    //Get x-ranges
    YAML::Node xmax = config["xmax"];
    YAML::Node xmin = config["xmin"];

    //Get directory 
    YAML::Node directory = config["Directory"];
    YAML::Emitter outdir;
    outdir << directory;
    std::string outdirstring = outdir.c_str();

    //For the case where all QA plots are compared:
    std::vector<std::string> objectNames;
    if(plotstring1=="All"){
        nplots=0;
        YAML::Emitter out;
        out << filepaths[0];
        std::string filename = out.c_str();
        //Read current file
        TFile *file = TFile::Open(filename.c_str(),"READ");
        //Navigate to directory
        TDirectory* dh = (TDirectory*)file->Get(outdirstring.c_str());
        //Get entries in tdirectory:
        TList *listofkeys = dh->GetListOfKeys();
        TIter next(listofkeys);
        TKey *key;
        // Loop over the keys and store the object names
        while ((key = (TKey*)next())) {
            // Get the name of the object
            std::string name = key->GetName();
            cout<< key->GetClassName()<<"\n";
            // Add the name to the vector

            if(strcmp(key->GetClassName(),"TH1F")==0){
                nplots+=1;
                objectNames.push_back(name);
                cout<<name<<"\n";
            }
            
        }
        cout<<"\n";
    }
    
    //std::string All[10]={"hIsoGammaIsoCharged","hIsoGammaIsoChargedCorrected","hIsoGammaPx","hIsoGammaPy","hIsoGammaPz","hIsoGammaE","hIsoGammaM02","hIsoGammaM20","hIsoGammaEtaPhi","hIsoGammaEBeforeAfterNL"};
    //if(plotstring1=="All"){
    //    nplots=10;
    //}

    std::string toplotstring = "";

    //Loop through plots to compare:
    for(Int_t j=0; j<nplots; j++){
        //Read name of current QAplot

        if(plotstring1=="All"){
            toplotstring = objectNames.at(j);
        }else{
           YAML::Emitter outplot;
            outplot << toplot[j];
            toplotstring = outplot.c_str(); 
        }
        
        //Initiate plot class:
        Plotting1D Phdum;
        TH1F* hdum;
        //Begin file loop
        for(Int_t i=0;i<nfiles;i++){
            
            //Read current filename
            YAML::Emitter out;
            out << filepaths[i];
            std::string filename = out.c_str();
            cout << filename <<"\n";

            //Read current fileIdent
            YAML::Emitter outid;
            outid << fileIdent[i];
            std::string fileident = outid.c_str();
            cout << fileident <<"\n";

            //Read current file
            TFile *file = TFile::Open(filename.c_str(),"READ");
            //file->ls();
            //Navigate to directory
            TDirectory* dh = (TDirectory*)file->Get(outdirstring.c_str());

            //Loop through groups:
            for(int k=0;k<ngroups;k++){
                YAML::Emitter outgroup;
                outgroup << groups[k];
                std::string groupstring = outgroup.c_str();

                //Initiate new histogram:
                TH1F* h = (TH1F*)dh->Get((toplotstring+groupstring).c_str())->Clone(toplotstring.c_str());
                float integralSignal = h->Integral(1, h->GetNbinsX());
                h->Scale(1. / integralSignal);
                //Adjust bin contents to avoid logy conflict.
                double min_nonzero = std::numeric_limits<double>::max();
                for (int i = 1; i <= h->GetNbinsX(); ++i) {
                    float content = h->GetBinContent(i);
                    if (content > 1e-10 && content < min_nonzero) {
                        min_nonzero = content;
                    }
                    if (h->GetBinContent(i) == 0) {
                        h->SetBinContent(i, 1e-10);
                    }
                }
                h->SetMinimum(min_nonzero);

                //Add to ouput vector:
                Phdum.New(h, fileident+groupstring);
            }

            
        }
        //Plot finalplots
        Phdum.Plot(Form("%s/%s.%s", path,toplotstring.c_str(), suffix),false,true);
    }
    
    


    return;
}