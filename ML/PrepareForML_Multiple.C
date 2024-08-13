#include "../Analysis/Utilities.h"
#include "../Analysis/Cuts.h"
#include "../Analysis/PhysicsObjects.h"
#include "../Analysis/TreeUtilities.h"
#include "../Analysis/HistogramLibrary.h"
#include "../Analysis/ClusterECorrections.h"
#include "TRandom.h"

void RestructureTreeForML_TreeUtils(TString AnalysisDirectory,TString inputfilename,GlobalOptions optns, bool bDoClusterCuts);
void DoClusterCuts(std::vector<IsoGamma> &IsoGammas, GlobalOptions optns);
bool PassedClusterCuts(IsoGamma IsoGamma, float EMCalEtaPhiMinMax[2][2], float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2], 
        float EMin, float EMax, unsigned short NcellsMin, unsigned short NLMMax, float DistanceToBadChannelMin, float MatchDetaMin, float MatchDetaMax, 
        float MatchDphiMin, float MatchDphiMax, float MatchVetoMax, float FplusMax);


void PrepareForML_Multiple(TString AnalysisDirectory,int GroupID, bool bDoClusterCuts){
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
        RestructureTreeForML_TreeUtils(AnalysisDirectory,TString(line.c_str()),optns, bDoClusterCuts);
    }

    inputFileNames.close();

    INFO(Form("Done with %s",Form("%s/%s/InputFiles/InputFiles_group_%i.txt",optns.dataSet.Data(),optns.trainConfig.Data(),GroupID)))
    EXIT
}

void RestructureTreeForML_TreeUtils(TString AnalysisDirectory,TString inputfilename,GlobalOptions optns,bool bDoClusterCuts){
//Input:
    //Read data:
    TFile* inputfile=TFile::Open(inputfilename.Data(),"READ");
    //Put input tree in tchain
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

    TString OutName;
    if(bDoClusterCuts){
        OutName=Form("%s/GammaIsoTree_%s_ClusterCuts.root",outputdir.Data(),optns.trainConfig.Data());
    }else{
        OutName=Form("%s/GammaIsoTree_%s.root",outputdir.Data(),optns.trainConfig.Data());
    }

    TFile* outputfile=TFile::Open(OutName,"RECREATE");
    INFO(Form("Created:%s",OutName.Data()));
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
    float Cluster_SplitFloat = 0;
    float Cluster_NLM = 0;                  //unsigned short
    float Cluster_SM = 0;                   //unsigned short
    float Cluster_NCells = 0;               //unsigned short
    float Cluster_MatchTrackIsConv = 0;      //bool

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
    caloclustertree->Branch("Cluster_NLM",&Cluster_NLM);
    caloclustertree->Branch("Cluster_SM",&Cluster_SM);
    caloclustertree->Branch("Cluster_NCells",&Cluster_NCells);
    caloclustertree->Branch("Cluster_MatchTrackIsConv",&Cluster_MatchTrackIsConv);
    caloclustertree->Branch("Cluster_SplitFloat",&Cluster_SplitFloat);

    bool Cluster_isSignal = 0;
    caloclustertree->Branch("Cluster_isSignal",&Cluster_isSignal);

//Load isogamma cuts, necessary to determine if cluster issignal:
    TDirectory* dummy= new TDirectory();
    IsoGammaCuts cuts= IsoGammaCuts(optns, dummy);

    //For generating random numbers for splitting into training and hyper param training:
    TRandom randGen;
    randGen.SetSeed(42);

//Loop over evets:
    for (int iEvent = 0; iEvent < tree.NEvents; iEvent++){
        //Get event:
        chain->GetEntry(iEvent);
        //Read into physics objects
        Event event(tree, optns);
        std::vector<IsoGamma> clusters;
        saveClustersFromEventInVector(tree, clusters, optns);
        calculateIsolation(clusters, event, false);//HARDCODE!!!!

        //INFO(Form("Before Cluster cuts: %i", clusters.size()))
        if(bDoClusterCuts){
            DoClusterCuts(clusters, optns);
        }
        //INFO(Form("After Cluster cuts: %i", clusters.size()))
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
            Cluster_SplitFloat=randGen.Uniform(0, 1);
            
            if(clusters.at(icluster).MatchedTrack.P>0){
                Cluster_MatchTrackdEta=TMath::Abs(clusters.at(icluster).MatchedTrack.dEta);
                Cluster_MatchTrackdPhi=TMath::Abs(clusters.at(icluster).MatchedTrack.dPhi);
                Cluster_MatchTrackP=clusters.at(icluster).MatchedTrack.P;
                Cluster_MatchTrackPt=clusters.at(icluster).MatchedTrack.Pt;
            }else{
                Cluster_MatchTrackdEta=10;
                Cluster_MatchTrackdPhi=10;
                Cluster_MatchTrackP=10;
                Cluster_MatchTrackPt=10;
            }
            
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

void DoClusterCuts(std::vector<IsoGamma> &IsoGammas, GlobalOptions optns){
//Load cuts
    YAML::Node ycut = YAML::LoadFile("Cuts.yaml");

    if (!ycut[(std::string)optns.cutString])
        FATAL(Form("Cutstring %s not found in YAML file Cuts.yaml", optns.cutString.Data()))

    YAML::Node standardcut = ycut["Standard"];
    YAML::Node chosencut = ycut[(std::string)optns.cutString];

    // Acceptance cut
    float EMCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
    float DCalEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
    float DCalHoleEtaPhiMinMax[2][2] = {{0, 0}, {0, 0}};
    // EMcal
    EMCalEtaPhiMinMax[0][0] = chosencut["cluster_min_EMcal_eta"].IsDefined() ? chosencut["cluster_min_EMcal_eta"].as<float>() : standardcut["cluster_min_EMcal_eta"].as<float>();
    EMCalEtaPhiMinMax[0][1] = chosencut["cluster_max_EMcal_eta"].IsDefined() ? chosencut["cluster_max_EMcal_eta"].as<float>() : standardcut["cluster_max_EMcal_eta"].as<float>();
    EMCalEtaPhiMinMax[1][0] = chosencut["cluster_min_EMcal_phi"].IsDefined() ? chosencut["cluster_min_EMcal_phi"].as<float>() : standardcut["cluster_min_EMcal_phi"].as<float>();
    EMCalEtaPhiMinMax[1][1] = chosencut["cluster_max_EMcal_phi"].IsDefined() ? chosencut["cluster_max_EMcal_phi"].as<float>() : standardcut["cluster_max_EMcal_phi"].as<float>();
    // DCal
    DCalEtaPhiMinMax[0][0] = chosencut["cluster_min_Dcal_eta"].IsDefined() ? chosencut["cluster_min_Dcal_eta"].as<float>() : standardcut["cluster_min_Dcal_eta"].as<float>();
    DCalEtaPhiMinMax[0][1] = chosencut["cluster_max_Dcal_eta"].IsDefined() ? chosencut["cluster_max_Dcal_eta"].as<float>() : standardcut["cluster_max_Dcal_eta"].as<float>();
    DCalEtaPhiMinMax[1][0] = chosencut["cluster_min_Dcal_phi"].IsDefined() ? chosencut["cluster_min_Dcal_phi"].as<float>() : standardcut["cluster_min_Dcal_phi"].as<float>();
    DCalEtaPhiMinMax[1][1] = chosencut["cluster_max_Dcal_phi"].IsDefined() ? chosencut["cluster_max_Dcal_phi"].as<float>() : standardcut["cluster_max_Dcal_phi"].as<float>();
    // DCal hole
    DCalHoleEtaPhiMinMax[0][0] = chosencut["cluster_min_DcalHole_eta"].IsDefined() ? chosencut["cluster_min_DcalHole_eta"].as<float>() : standardcut["cluster_min_DcalHole_eta"].as<float>();
    DCalHoleEtaPhiMinMax[0][1] = chosencut["cluster_max_DcalHole_eta"].IsDefined() ? chosencut["cluster_max_DcalHole_eta"].as<float>() : standardcut["cluster_max_DcalHole_eta"].as<float>();
    DCalHoleEtaPhiMinMax[1][0] = chosencut["cluster_min_DcalHole_phi"].IsDefined() ? chosencut["cluster_min_DcalHole_phi"].as<float>() : standardcut["cluster_min_DcalHole_phi"].as<float>();
    DCalHoleEtaPhiMinMax[1][1] = chosencut["cluster_max_DcalHole_phi"].IsDefined() ? chosencut["cluster_max_DcalHole_phi"].as<float>() : standardcut["cluster_max_DcalHole_phi"].as<float>();
    // Load E cut
    float EMin = chosencut["cluster_min_E"].IsDefined() ? chosencut["cluster_min_E"].as<float>() : standardcut["cluster_min_E"].as<float>();
    float EMax = chosencut["cluster_max_E"].IsDefined() ? chosencut["cluster_max_E"].as<float>() : standardcut["cluster_max_E"].as<float>();

    // Load min Ncells
    unsigned short NcellsMin = chosencut["cluster_min_Nc"].IsDefined() ? chosencut["cluster_min_Nc"].as<unsigned short>() : standardcut["cluster_min_Nc"].as<unsigned short>();

    // Load NLM cut
    unsigned short NLMMax = chosencut["cluster_max_NLM"].IsDefined() ? chosencut["cluster_max_NLM"].as<unsigned short>() : standardcut["cluster_max_NLM"].as<unsigned short>();

    // Load min dist to bad channel cut
    float DistanceToBadChannelMin = chosencut["cluster_min_distbadch"].IsDefined() ? chosencut["cluster_min_distbadch"].as<float>() : standardcut["cluster_min_distbadch"].as<float>();

    // Load track matching cut
    float MatchDetaMin = chosencut["cluster_min_MatchDeta"].IsDefined() ? chosencut["cluster_min_MatchDeta"].as<float>() : standardcut["cluster_min_MatchDeta"].as<float>();
    float MatchDetaMax = chosencut["cluster_max_MatchDeta"].IsDefined() ? chosencut["cluster_max_MatchDeta"].as<float>() : standardcut["cluster_max_MatchDeta"].as<float>();
    float MatchDphiMin = chosencut["cluster_min_MatchDeta"].IsDefined() ? chosencut["cluster_min_MatchDeta"].as<float>() : standardcut["cluster_min_MatchDeta"].as<float>();
    float MatchDphiMax = chosencut["cluster_max_MatchDeta"].IsDefined() ? chosencut["cluster_max_MatchDeta"].as<float>() : standardcut["cluster_max_MatchDeta"].as<float>();
    float MatchVetoMax = chosencut["cluster_max_MatchVeto"].IsDefined() ? chosencut["cluster_max_MatchVeto"].as<float>() : standardcut["cluster_max_MatchVeto"].as<float>();

    // Load FM+
    float FplusMax = chosencut["cluster_max_Fplus"].IsDefined() ? chosencut["cluster_max_Fplus"].as<float>() : standardcut["cluster_max_Fplus"].as<float>();

//Apply cuts:
    std::vector<IsoGamma>::iterator iter;
    for (iter = IsoGammas.begin(); iter != IsoGammas.end();)
    {
        if (!PassedClusterCuts(*iter,EMCalEtaPhiMinMax, DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax, 
        EMin, EMax, NcellsMin, NLMMax, DistanceToBadChannelMin, MatchDetaMin, MatchDetaMax, 
        MatchDphiMin, MatchDphiMax, MatchVetoMax, FplusMax))
        {
            iter = IsoGammas.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

bool PassedClusterCuts(IsoGamma IsoGamma, float EMCalEtaPhiMinMax[2][2], float DCalEtaPhiMinMax[2][2], float DCalHoleEtaPhiMinMax[2][2], 
        float EMin, float EMax, unsigned short NcellsMin, unsigned short NLMMax, float DistanceToBadChannelMin, float MatchDetaMin, float MatchDetaMax, 
        float MatchDphiMin, float MatchDphiMax, float MatchVetoMax, float FplusMax)
{
  bool passed = true;
  // Check cluster acceptance
  if (!IsoGamma.isInEMCalAcceptance(EMCalEtaPhiMinMax) && !IsoGamma.isInDCalAcceptance(DCalEtaPhiMinMax, DCalHoleEtaPhiMinMax))
  {
    passed = false;
  }
  // check cluster energy
  if (IsoGamma.E < EMin || IsoGamma.E > EMax)
  {
    passed = false;
  }
  // check cells pr. cluster
  //if (IsoGamma.NCells < NcellsMin)
  //{
  //  passed = false;
  //}
  // Check number of local maxima
  //if (IsoGamma.NLM > NLMMax)
  //{
  //  passed = false;
  //}
  // Check distance to bad channel
  //if (IsoGamma.DistanceToBadChannel < DistanceToBadChannelMin)
  //{
  //  passed = false;
  //}
   //Check track matching
  if (IsoGamma.MatchedTrack.P>0 && IsoGamma.MatchedTrack.dEta < MatchDetaMax && IsoGamma.MatchedTrack.dEta > MatchDetaMin && IsoGamma.MatchedTrack.dPhi > MatchDphiMin && IsoGamma.MatchedTrack.dPhi < MatchDphiMax && (IsoGamma.E / IsoGamma.MatchedTrack.P) < MatchVetoMax)
  {
    passed = false;
  }
  // Check Fplus
  //if (IsoGamma.EFrac > FplusMax)
  //{
  //  passed = false;
  //}
  return passed;
}