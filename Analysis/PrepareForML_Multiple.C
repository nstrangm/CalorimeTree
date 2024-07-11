#include "Utilities.h"
#include "Cuts.h"
#include "PhysicsObjects.h"
#include "TreeUtilities.h"
#include "HistogramLibrary.h"
#include "ClusterECorrections.h"

void RestructureTreeForML_TreeUtils(TString AnalysisDirectory,TString inputfilename,GlobalOptions optns);


void PrepareForML_Multiple(TString AnalysisDirectory,int GroupID){
    ENTER
    INFO("started PrepareForML_Multiple")
    if (GroupID < 0){
        FATAL("Negative GroupID")
    }
        

    GlobalOptions optns(AnalysisDirectory, 1);

    //Open InputFile names:
    ifstream inputFileNames(Form("%s/%s/InputFiles/InputFiles_group_%i.txt",optns.dataSet.Data(),optns.trainConfig.Data(),GroupID));//

    if(!inputFileNames.is_open()){
        FATAL(Form("Specified inputFilenamesFile=%s does not exist. Either analysis directory or GroupID is wrong!",Form("%s/%s/InputFiles/InputFiles_group_%i.txt",optns.dataSet.Data(),optns.trainConfig.Data(),GroupID)));
        return;
    }
    //Read Inputfilenamesfile:
    string line;
    cout<<"File content:"<< endl;
    while (getline(inputFileNames,line)){
        //cout<<line<<endl;
        RestructureTreeForML_TreeUtils(AnalysisDirectory,TString(line.c_str()),optns);
    }

    inputFileNames.close();

    INFO(Form("Done with %s",Form("%s/%s/InputFiles/InputFiles_group_%i.txt",optns.dataSet.Data(),optns.trainConfig.Data(),GroupID)))
    EXIT
    //
}

void RestructureTreeForML_TreeUtils(TString AnalysisDirectory,TString inputfilename,GlobalOptions optns){
//Input:
    //Read data:
    //TString tempLine=Form("%s/GammaIsoTree_%s.root",inputdir.Data(),triggerstring.Data());
    TFile* inputfile=TFile::Open(inputfilename.Data(),"READ");
    //Put input tree in tchain
    //TString treeName = GetTreeName(tempLine.Data());
    TChain *chain = new TChain("CaloTree_V1");
    if ((!inputfile) || inputfile->IsZombie()){
      ERROR(Form("Input file %s is not healthy! Skipping this one...", inputfilename.Data()))
      inputfile->Close();
      return;
    }
    inputfile->Close();
    INFO(Form("Adding tree %s from %s to chain", "CaloTree_V1", inputfilename.Data()))
    chain->Add(inputfilename.Data());
    optns.TreeFormat = listTreeBranches(chain);
    TreeBuffer tree(chain, optns);
//Output:
    //Create output dir:
    TString outputdir=inputfilename;
    if(optns.isMC){
        TString folderreplace="GA_pp_MC_AOD";
        Ssiz_t index = outputdir.Index(folderreplace);
        if(index!=kNPOS){
            outputdir.Replace(index, folderreplace.Length(), "ML_pp_MC_AOD");
        }
        TString filenamereplace=Form("GammaIsoTree_%s.root",optns.trainConfig.Data());
        Ssiz_t indexname = outputdir.Index(filenamereplace);
        if(indexname!=kNPOS){
            TString outputdirfinal=outputdir.Replace(indexname, filenamereplace.Length(), "");
            cout << outputdir <<"\n";
        }
    }else{
        TString folderreplace="GA_pp_AOD";
        Ssiz_t index = outputdir.Index(folderreplace);
        if(index!=kNPOS){
            outputdir.Replace(index, folderreplace.Length(), "ML_pp_AOD");
        }
        TString filenamereplace=Form("GammaIsoTree_%s.root",optns.trainConfig.Data());
        Ssiz_t indexname = outputdir.Index(filenamereplace);
        if(indexname!=kNPOS){
            TString outputdirfinal=outputdir.Replace(indexname, filenamereplace.Length(), "");
            cout << outputdir <<"\n";
        }
    }
    
    
    //createDirectory(outputdir);
    std::filesystem::create_directories(outputdir.Data());
    TFile* outputfile=TFile::Open(Form("%s/GammaIsoTree_%s.root",outputdir.Data(),optns.trainConfig.Data()),"RECREATE");
    if ((!outputfile) || outputfile->IsZombie()){
      ERROR(Form("Input file %s is not healthy! Skipping this one...", inputfilename.Data()))
      inputfile->Close();
      return;
    }
    //Create outputtree:
    TTree* caloclustertree = new TTree(Form("caloclustertree"),Form("caloclustertree"));
    if (!caloclustertree) {
        ERROR(Form("Input tree CaloTree_V1 in file %s does not exist! Skipping this one...", inputfilename.Data()));
        inputfile->Close();
        return;
    }
    //Define Containers for output:
//Event attributes:
    float Event_Rho = 0;
    float Event_ZVtx = 0;
    unsigned short Event_NPrimaryTracks = 0;

    caloclustertree->Branch("Event_Rho",&Event_Rho);
    caloclustertree->Branch("Event_ZVtx",&Event_ZVtx);
    caloclustertree->Branch("Event_NPrimaryTracks",&Event_NPrimaryTracks);

//Cluster attributes:
    float Cluster_E = 0;
    float Cluster_Pt = 0;
    float Cluster_M02 = 0;
    float Cluster_M20 = 0;
    float Cluster_V1SplitMass = 0;
    float Cluster_MinMassDiffToPi0 = 0;
    float Cluster_MinMassDiffToEta = 0;
    float Cluster_EFrac = 0;
    float Cluster_IsoGammaCorrected = 0;
    float Cluster_MatchTrackdEta = 0;
    float Cluster_MatchTrackdPhi = 0;
    float Cluster_MatchTrackP = 0;
    float Cluster_MatchTrackPt = 0;
    float Cluster_DistanceToBadChannel = 0; 
    unsigned short Cluster_NLM = 0;
    unsigned short Cluster_SM = 0;
    unsigned short Cluster_NCells = 0;
    bool Cluster_MatchTrackIsConv = 0;

    caloclustertree->Branch("Cluster_E",&Cluster_E);
    caloclustertree->Branch("Cluster_Pt",&Cluster_Pt);
    caloclustertree->Branch("Cluster_M02",&Cluster_M02);
    caloclustertree->Branch("Cluster_M20",&Cluster_M20);
    caloclustertree->Branch("Cluster_V1SplitMass",&Cluster_V1SplitMass);
    caloclustertree->Branch("Cluster_MinMassDiffToPi0",&Cluster_MinMassDiffToPi0);
    caloclustertree->Branch("Cluster_MinMassDiffToEta",&Cluster_MinMassDiffToEta);
    caloclustertree->Branch("Cluster_EFrac",&Cluster_EFrac);
    caloclustertree->Branch("Cluster_IsoGammaCorrected",&Cluster_IsoGammaCorrected);
    caloclustertree->Branch("Cluster_MatchTrackdEta",&Cluster_MatchTrackdEta);
    caloclustertree->Branch("Cluster_MatchTrackdPhi",&Cluster_MatchTrackdPhi);
    caloclustertree->Branch("Cluster_MatchTrackP",&Cluster_MatchTrackP);
    caloclustertree->Branch("Cluster_MatchTrackPt",&Cluster_MatchTrackPt);
    caloclustertree->Branch("Cluster_DistanceToBadChannel",&Cluster_DistanceToBadChannel);
    caloclustertree->Branch("Cluster_Cluster_NLM",&Cluster_NLM);
    caloclustertree->Branch("Cluster_Cluster_SM",&Cluster_SM);
    caloclustertree->Branch("Cluster_NCells",&Cluster_NCells);
    caloclustertree->Branch("Cluster_SM",&Cluster_SM);
    caloclustertree->Branch("Cluster_NCells",&Cluster_NCells);
    caloclustertree->Branch("Cluster_MatchTrackIsConv",&Cluster_MatchTrackIsConv);

    bool Cluster_isSignal = 0;
    caloclustertree->Branch("Cluster_isSignal",&Cluster_isSignal);

//Load isogamma cuts, necessary to determine if cluster issignal:
    TDirectory* dummy= new TDirectory();
    IsoGammaCuts cuts= IsoGammaCuts(optns, dummy);

//Loop over evets:
    for (int iEvent = 0; iEvent < tree.NEvents; iEvent++){
        //Get event:
        chain->GetEntry(iEvent);
        //Read into physics objects
        Event event(tree, optns);
        std::vector<IsoGamma> clusters;
        saveClustersFromEventInVector(tree, clusters, optns);
        calculateIsolation(clusters, event, false);//HARDCODE!!!!
        for (unsigned long icluster=0; icluster<clusters.size(); icluster++){
            //Fill event attributes:
            Event_Rho=event.Rho;
            Event_ZVtx=event.ZVtx;
            Event_NPrimaryTracks=event.NPrimaryTracks;
            //fill cluster attributes
            Cluster_E=clusters.at(icluster).E;
            Cluster_Pt=clusters.at(icluster).Pt();
            Cluster_M02=clusters.at(icluster).M02;
            Cluster_M20=clusters.at(icluster).M20;
            Cluster_V1SplitMass=clusters.at(icluster).V1SplitMass;
            Cluster_MinMassDiffToPi0=clusters.at(icluster).MinMassDiffToPi0;
            Cluster_MinMassDiffToEta=clusters.at(icluster).MinMassDiffToEta;
            Cluster_EFrac=clusters.at(icluster).EFrac;
            Cluster_IsoGammaCorrected=clusters.at(icluster).IsoChargedCorrected;
            Cluster_MatchTrackdEta=clusters.at(icluster).MatchedTrack.dEta;
            Cluster_MatchTrackdPhi=clusters.at(icluster).MatchedTrack.dPhi;
            Cluster_MatchTrackP=clusters.at(icluster).MatchedTrack.P;
            Cluster_MatchTrackPt=clusters.at(icluster).MatchedTrack.Pt;
            Cluster_DistanceToBadChannel=clusters.at(icluster).DistanceToBadChannel;
            Cluster_NLM=clusters.at(icluster).NLM;
            Cluster_NCells=clusters.at(icluster).NCells;
            Cluster_SM=clusters.at(icluster).SM;

            if(optns.isMC){
                Cluster_isSignal = cuts.isSignal(clusters.at(icluster));
            }

            caloclustertree->Fill();
        }
        clusters.clear();
    }
    INFO(Form("Size of input file: %llu",inputfile->GetSize()));
    inputfile->Close();
    outputfile->Write();
    INFO(Form("Size of output file: %llu",outputfile->GetSize()));
    outputfile->Close();
    delete chain;
    delete dummy;
}