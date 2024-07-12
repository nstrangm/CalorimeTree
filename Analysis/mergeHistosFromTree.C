#include "../Analysis/Utilities.h"
#include "../Analysis/Cuts.h"
#include "../Analysis/PhysicsObjects.h"
#include "../Analysis/TreeUtilities.h"
#include "../Analysis/HistogramLibrary.h"
#include "../Analysis/ClusterECorrections.h"

/*
This script reads two root files containing histograms.
It then reads the contents of both and merges the contents
of the directories that appear in both files.
*/

std::set<std::string> getDirectories(TFile *file)
{
    std::set<std::string> directories;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        if (std::string(key->GetClassName()) == "TDirectoryFile")
        {
            directories.insert(key->GetName());
        }
    }
    return directories;
}

void mergeHistosFromTree(TString outputDir = "MergedGJ18DandJJ18D/100_Charged/Standard/", string filepath1 = "GJ18D/100_Charged/Standard", string filepath2 = "JJ18D/100_Charged/Standard", float w1 = 1, float w2 = 1,bool addtoRunConfig=true)
{
    ENTER
    // Open files
    TFile *file1 = TFile::Open(Form("%s/HistosFromTree.root",filepath1.c_str()), "READ");
    TFile *file2 = TFile::Open(Form("%s/HistosFromTree.root",filepath2.c_str()), "READ");

    GlobalOptions optns1(filepath1.c_str(),1);
    GlobalOptions optns2(filepath2.c_str(),1);

    if(optns1.isMC!=optns2.isMC){
        FATAL("Trying to merge MC data and non-MC data!");
    }

    //Add merge file config to RunConfig
    if(addtoRunConfig){
        //Load configuration file:
        std::stringstream ss((std::string)outputDir.Data());
        std::string word;
        std::string words[3];
        for (int i = 0; i < 3; i++)
        {
            std::getline(ss, word, '/');
            words[i] = word;
        }
        YAML::Node config = YAML::LoadFile("RunConfig.yaml");
        YAML::Node mergeconfig = config[(words[0])];
        //Is MC?
        mergeconfig["isMC"]=optns1.isMC;
        //Paths to raw datasets
        mergeconfig["path1"]=std::string(Form("%s, %s",optns1.analysisDirPath.Data(),optns2.analysisDirPath.Data()));
        //Labels
        mergeconfig["label"]=Form("%s, %s",optns1.dataSet.Data(),optns2.dataSet.Data());
        //trainconfigs
        mergeconfig["TrainConfig1"]=optns1.trainConfig.Data();
        mergeconfig["TrainConfig2"]=optns2.trainConfig.Data();
        std::ofstream fout("RunConfig.yaml");
        fout << config;
    }

    createDirectory(outputDir.Data());

    std::set<std::string> directories1 = getDirectories(file1);
    std::set<std::string> directories2 = getDirectories(file2);

    // Find common directories
    std::set<std::string> commonDirectories;
    for (const auto &dir : directories1)
    {
        if (directories2.find(dir) != directories2.end())
        {
            commonDirectories.insert(dir);
        }else{
            WARN(Form("%s in directory 1, but not in directory 2.",dir.c_str()))
        }
    }

    // Loop through common directories

    TFile *TestOut = TFile::Open(Form("%sHistosFromTree.root", outputDir.Data()), "Recreate");
    for (const auto &dir : commonDirectories)
    {
        TDirectory *dh1 = (TDirectory *)file1->Get(dir.c_str());
        TDirectory *dh2 = (TDirectory *)file2->Get(dir.c_str());

        LOG(Form("Merging histograms from the directory: %s",dir.c_str()))

        // Get list of keys1:
        TList *listofkeys1 = dh1->GetListOfKeys();
        TList *listofkeys2 = dh2->GetListOfKeys();

        TDirectory *outdir = TestOut->mkdir(dir.c_str());

        TIter next(listofkeys1);
        TKey *key;

        std::string toplotstring = "";
        std::vector<std::string> objectNames;
        std::vector<std::string> objectClass;
        int nplots = 0;
        while ((key = (TKey *)next()))
        {
            // Get the name of the object
            std::string name = key->GetName();
            std::string classname = key->GetClassName();
            cout<<name<<"\n";
            cout<<classname<<"\n";
            //  Add the name to the vector
            if ((strcmp(key->GetClassName(), "TH1F") == 0)||(strcmp(key->GetClassName(), "TH2F") == 0)||(strcmp(key->GetClassName(), "THnSparseT<TArrayF>") == 0))
            {
                nplots += 1;
                objectNames.push_back(name);
                objectClass.push_back(classname);
            }
        }
        for (int i = 0; i < nplots; i++)
        {
            if(strcmp(objectClass.at(i).c_str(), "TH1F") == 0){
                TH1F *h1 = (TH1F *)dh1->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                TH1F *h2 = (TH1F *)dh2->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                if(!h1){
                    cout<<"h1 not loaded"<<"\n";
                }
                if(!h2){
                    cout<<"h2 not loaded"<<"\n";
                }
                h1->Scale(w1);
                h1->Add(h2, w2);
                outdir->cd();
                h1->Write();

            }
            if(strcmp(objectClass.at(i).c_str(), "TH2F") == 0){
                TH2F *h1 = (TH2F *)dh1->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                TH2F *h2 = (TH2F *)dh2->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                if(!h1){
                    cout<<"h1 not loaded"<<"\n";
                }
                if(!h2){
                    cout<<"h2 not loaded"<<"\n";
                }
                h1->Scale(w1);
                h1->Add(h2, w2);
                outdir->cd();
                h1->Write();
            }
            if(strcmp(objectClass.at(i).c_str(), "THnSparseT<TArrayF>") == 0){
                THnSparseT<TArrayF> *h1 = (THnSparseT<TArrayF> *)dh1->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                THnSparseT<TArrayF> *h2 = (THnSparseT<TArrayF> *)dh2->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
                if(!h1){
                    cout<<"h1 not loaded"<<"\n";
                }
                if(!h2){
                    cout<<"h2 not loaded"<<"\n";
                }
                h1->Scale(w1);
                h1->Add(h2, w2);
                outdir->cd();
                h1->Write();
            }
        }
    }
    TestOut->Close();
    file1->Close();
    file2->Close();

    EXIT
}