#include "../Analysis/Utilities.h"
#include "../Analysis/Cuts.h"
#include "../Analysis/PhysicsObjects.h"
#include "../Analysis/TreeUtilities.h"
#include "../Analysis/HistogramLibrary.h"
#include "../Analysis/ClusterECorrections.h"

/*
This script reads two root files containing histograms.
It then reads the contents of both and merges the contents
of those files that appear in both files.
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

void mergeHistosFromTree(TString outputDir = "MergedGJ18DandJJ18D/100_Charged/Standard/", string filepath1 = "GJ18D/100_Charged/Standard/HistosFromTree.root", string filepath2 = "JJ18D/100_Charged/Standard/HistosFromTree.root", float w1 = 1, float w2 = 1)
{
    // Open files
    TFile *file1 = TFile::Open(filepath1.c_str(), "READ");
    TFile *file2 = TFile::Open(filepath2.c_str(), "READ");

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
        }
    }

    // Loop through common directories

    TFile *TestOut = TFile::Open(Form("%sHistosFromTree.root", outputDir.Data()), "Recreate");
    for (const auto &dir : commonDirectories)
    {
        TDirectory *dh1 = (TDirectory *)file1->Get(dir.c_str());
        TDirectory *dh2 = (TDirectory *)file2->Get(dir.c_str());

        // Get list of keys1:
        TList *listofkeys1 = dh1->GetListOfKeys();
        TList *listofkeys2 = dh2->GetListOfKeys();

        TDirectory *outdir = TestOut->mkdir(dir.c_str());

        TIter next(listofkeys1);
        TKey *key;

        std::string toplotstring = "";
        std::vector<std::string> objectNames;
        int nplots = 0;
        while ((key = (TKey *)next()))
        {
            // Get the name of the object
            std::string name = key->GetName();
            // cout<< key->GetClassName()<<"\n";
            // cout<< dh2->FindKey(name.c_str())<<"\n";
            //  Add the name to the vector
            if (strcmp(key->GetClassName(), "TH1F") == 0)
            {
                nplots += 1;
                objectNames.push_back(name);
                cout << name << "\n";
            }
        }
        for (int i = 0; i < nplots; i++)
        {
            TH1 *h1 = (TH1 *)dh1->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());
            TH1 *h2 = (TH1 *)dh2->Get(objectNames.at(i).c_str())->Clone(objectNames.at(i).c_str());

            // Apply weight to h1:
            h1->Scale(w1);
            h1->Add(h2, w2);
            // h1->Draw();
            outdir->cd();
            h1->Write();
            outdir->ls();
            // if(i==nplots-1){
            //     h1->Draw();
            //     h2->Draw("same");
        }

        // Loop through objects
    }

    // Get directory of histograms
    // TDirectory* dh1 = (TDirectory*)file1->Get(outdirstring.c_str());
    // TDirectory* dh2 = (TDirectory*)file1->Get(outdirstring.c_str());
    ////Get entries in tdirectory:
    // TList *listofkeys1 = dh1->GetListOfKeys();
    // TList *listofkeys2 = dh2->GetListOfKeys();
    ////Keep only entries that are in both dh1 and dh2.
    // string Directory = "IsoGammaQA";

    // Loop over the keys and store the object names
    // TIter next(listofkeys1);
    // TKey *key;
    // while ((key = (TKey*)next())) {
    //    // Get the name of the object
    //    std::string name = key->GetName();
    //    cout<< key->GetClassName()<<"\n";
    //    cout<< dh2->FindKey(name.Data())<<"\n";
    //    // Add the name to the vector
    //    if(strcmp(key->GetClassName(),"TH1F")==0){
    //        nplots+=1;
    //        objectNames.push_back(name);
    //        cout<<name<<"\n";
    //    }
    //
    //}
}