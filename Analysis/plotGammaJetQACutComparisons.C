#include "PlottingClass.h"

std::vector<Float_t> jetPt = {10,20,30,40,50,60,70,80,90,100};


std::pair<TString,TString> extractCentAndRInfoFromFile(TString filename){
    TString cent = "";
    TString r = "";
    TObjArray *tokens = filename.Tokenize("/");
    for (int i = 0; i < tokens->GetEntries(); i++){
        TString token = ((TObjString*)tokens->At(i))->GetString();
        if (token.Contains("Run3-")){
            // cent is everything minus the first 5 characters
            cent = token.Remove(0,5);
        }
        if (token.Contains("JetRadius_")){
            r = token.Remove(0,10);
        }
        if (token.Contains("Standard")){
            r = "R04";
        }
        if (token.Contains("PerpCone")){
            r = "R09"; // dummy
        }
    }
    // convert string to int
    float fR  = std::stof(r.Remove(0,1).Data());
    r = Form("%1.1f",fR/10);
    return std::make_pair(cent,r);
}
void plotGammaJetQACutComparisons(TString inputfolder="../LHC23zzm_pass4", bool normalize = true){

    // define colors
    Int_t GoetheColor[4][2];
    Float_t r[4][2] = {{77, 228}, {115, 165}, {0, 72}, {134, 173}};
    Float_t g[4][2] = {{75, 227}, {124, 171}, {97, 143}, {0, 59}};
    Float_t b[4][2] = {{70, 221}, {69, 82}, {143, 218}, {71, 118}};
    for (Int_t iC = 0; iC < 4; iC++)
    {
        for (Int_t iH = 0; iH < 2; iH++)
        {
            Int_t ColorCode = TColor::GetFreeColorIndex();
            TColor *color = new TColor(ColorCode, r[iC][iH] / 255., g[iC][iH] / 255., b[iC][iH] / 255.);
            GoetheColor[iC][iH] = ColorCode;
        }
    }

    // define fancy colors from hex codes
//[#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"]



TString fancyColorHex[8] = {"#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff"};
Int_t fancyColors[8];
for (int i = 0; i < 8; i++){
    fancyColors[i] = TColor::GetColor(fancyColorHex[i]);
}

    std::vector<TString> excludeList = {"NoOccupancy","GammaV1","SigTT","RestrictedAcceptance","SigTT"};

    TString outputFolder = Form("%s/GammaJetComparisonPlots",inputfolder.Data());
    // create directory if it does not exist
    if (!std::filesystem::exists(outputFolder.Data())){
        std::filesystem::create_directory(outputFolder.Data());
    }

    std::vector<std::vector<TFile*>> files;
    // loop over all folders in inputfolder that contain train configurations (one per folder)
    // then loop over all folders inside that folder, to get all Cut configurations.
    // Inside the Cut configuration folder, there should be a file called "HistosFromTree.root which should be pushed into a vector of TFiles.
    
    // convert inputfolder to absolute path

    for (const auto & entry : std::filesystem::directory_iterator(inputfolder.Data())){
        if (entry.is_directory()){ // centrality
            // if folder does not start with Run3, continue
            if (entry.path().filename().string().find("Run3") == std::string::npos) continue;
            std::vector<TFile*> tempfiles;
            for (const auto & entry2 : std::filesystem::directory_iterator(entry.path())){   // jet radius
                // if it is a directory and pathname is not "InputFiles"
                if (entry2.is_directory() && entry2.path().filename() != "InputFiles"){
                    // check if folder name contains a string from the excludeList
                    bool exclude = false;
                    for (const auto & excludeString : excludeList){
                        if (entry2.path().filename().string().find(excludeString) != std::string::npos) exclude = true;
                    }
                    if (exclude) continue;
                    TFile* file = (TFile*)TFile::Open(Form("%s/HistosFromTree.root",entry2.path().c_str()));
                    // check if file exists, otherwise continue
                    if (!file) continue;
                    tempfiles.push_back(file);
                    
                }
            }
            files.push_back(tempfiles);
        }
    }

    // loop over all files in vector and print path
    cout << "-> I found the following files:" << endl;
    for (auto & filevec : files){
        for (auto & file : filevec){
            cout << file->GetName() << endl;
        }
    }

    // sort the vector inside so that radius is in ascending order
    for (int i = 0; i < files.size(); i++){
        std::sort(files[i].begin(),files[i].end(),[](TFile* a, TFile* b){
            TString rA = extractCentAndRInfoFromFile(a->GetName()).second;
            TString rB = extractCentAndRInfoFromFile(b->GetName()).second;
            return std::stof(rA.Data()) < std::stof(rB.Data());
        });
    }

    // loop over all files in vector and get the THnSparseF called
    // hJetPtEtaPhi from direcotry called JetQA
    THnSparseF* vJetPtEtaPhi[files.size()][files[0].size()];
    std::pair<TString,TString> labelCentR[files.size()][files[0].size()];
    for (int i = 0; i < files.size(); i++){
        for (int j = 0; j < files[i].size(); j++){
            files[i][j]->cd("JetQA");
            vJetPtEtaPhi[i][j] = (THnSparseF*)gDirectory->Get("hJetPtEtaPhi");
            labelCentR[i][j] = extractCentAndRInfoFromFile(files[i][j]->GetName());
            if(labelCentR[i][j].second == "0.9"){
                labelCentR[i][j].second = "0.4-PerpConeRho";
            }
        }
    }

    // make sure that R (which is a float stored as string in the second element of the pair) is sorted ascending, together with vJetPtEtaPhi


    // do a comparison plot of the eta distribution (integrated in pt) for all centralities and same R

    // loop over all centralities
    Plotting1D pEta[files.size()];
    for (int i = 0; i < files.size(); i++){ // loop over all centralities
        // loop over all R
        for (int j = 0; j < files[i].size(); j++){
            // get the eta distribution for all jets
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            if(normalize) hEta->Scale(1./hEta->Integral());
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
            pEta[i].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEta[i].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEta[i].New(hEta,Form("R = %s",labelCentR[i][j].second.Data()), 1, 1, fancyColors[j], "hist");
        }
        pEta[i].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; %s%% centrality",labelCentR[i][0].first.Data()), StdTextSize, 0.04);
        pEta[i].Plot(Form("%s/EtaComparison_%s.pdf",outputFolder.Data(),labelCentR[i][0].first.Data()));

    }

    Plotting1D pEtaCent[files[0].size()];
    // do the same plot as above but for all R and same centrality
    for (int j = 0; j < files[0].size(); j++){
        for (int i = 0; i < files.size(); i++){
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            if(normalize) hEta->Scale(1./hEta->Integral());
            hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
            pEtaCent[j].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEtaCent[j].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEtaCent[j].New(hEta,labelCentR[i][j].first, 1, 1, fancyColors[i], "hist");
        }
        pEtaCent[j].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; anti k_{T} R = %s",labelCentR[0][j].second.Data()), StdTextSize, 0.04);
        pEtaCent[j].Plot(Form("%s/EtaComparisonRadius_%s.pdf",outputFolder.Data(),labelCentR[0][j].second.Data()));
    }
    
    Double_t highPt[2] = {40,200};

    Double_t lowPt[2] = {10,20};
    // do the same figure but one only for high pt jets pt > 40 GeV/c
    Plotting1D pEtaHighPt[files.size()];
    for (int i = 0; i < files.size(); i++){ // loop over all centralities
        // loop over all R
        for (int j = 0; j < files[i].size(); j++){
            // get the eta distribution for all jets
            vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(highPt[0],highPt[1]);
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            if(normalize) hEta->Scale(1./hEta->Integral());
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            // pEtaHighPt[i].SetAxisRange(-0.81,0.81);
            pEtaHighPt[i].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEtaHighPt[i].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEtaHighPt[i].New(hEta,Form("R = %s",labelCentR[i][j].second.Data()), 1, 1, fancyColors[j], "hist");
        }
        pEtaHighPt[i].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; %s%% centrality; %.0f < #it{p}_{T} < %.0f GeV/#it{c}",labelCentR[i][0].first.Data(),highPt[0],highPt[1]), StdTextSize, 0.04);
        pEtaHighPt[i].Plot(Form("%s/EtaComparisonHighPt_%.0f_%.0f_%s.pdf",outputFolder.Data(), highPt[0],highPt[1],labelCentR[i][0].first.Data()));

    }

    Plotting1D pEtaCentHighPt[files[0].size()];
    for (int j = 0; j < files[0].size(); j++){
        for (int i = 0; i < files.size(); i++){
            vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(highPt[0],highPt[1]);
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
            if(normalize) hEta->Scale(1./hEta->Integral());
            pEtaCentHighPt[j].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEtaCentHighPt[j].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEtaCentHighPt[j].New(hEta,labelCentR[i][j].first, 1, 1, fancyColors[i], "hist");
        }
        pEtaCentHighPt[j].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; anti k_{T} R = %s; %.0f < #it{p}_{T} < %.0f GeV/#it{c}",labelCentR[0][j].second.Data(),highPt[0],highPt[1]), StdTextSize, 0.04);
        pEtaCentHighPt[j].Plot(Form("%s/EtaComparisonRadiusHighPt_%s.pdf",outputFolder.Data(),labelCentR[0][j].second.Data()));
    }

    // do the same figure but one only for low pt jets
    Plotting1D pEtaLowPt[files.size()];
    for (int i = 0; i < files.size(); i++){ // loop over all centralities
        // loop over all R
        for (int j = 0; j < files[i].size(); j++){
            // get the eta distribution for all jets
            vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(lowPt[0],lowPt[1]);
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
            if(normalize) hEta->Scale(1./hEta->Integral());
            pEtaLowPt[i].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEtaLowPt[i].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEtaLowPt[i].New(hEta,Form("R = %s",labelCentR[i][j].second.Data()), 1, 1, fancyColors[j], "hist");
        }
        pEtaLowPt[i].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; %s%% centrality; %.0f < #it{p}_{T} < %.0f GeV/#it{c}",labelCentR[i][0].first.Data(),lowPt[0],lowPt[1]), StdTextSize, 0.04);
        pEtaLowPt[i].Plot(Form("%s/EtaComparisonLowPt_%.0f_%.0f_%s.pdf",outputFolder.Data(), lowPt[0],lowPt[1],labelCentR[i][0].first.Data()));

    }

    Plotting1D pEtaCentLowPt[files[0].size()];
    for (int j = 0; j < files[0].size(); j++){
        for (int i = 0; i < files.size(); i++){
            vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(lowPt[0],lowPt[1]);
            TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
            hEta->SetName(Form("hEta_%s_%s",labelCentR[i][j].first.Data(),labelCentR[i][j].second.Data()));
            if(normalize) hEta->Scale(1./hEta->Integral());
            hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
            pEtaCentLowPt[j].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
            pEtaCentLowPt[j].SetLegend(0.15, 0.3, 0.70, 0.95, true);
            pEtaCentLowPt[j].New(hEta,labelCentR[i][j].first, 1, 1, fancyColors[i], "hist");
        }
        pEtaCentLowPt[j].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; anti k_{T} R = %s; %.0f < #it{p}_{T} < %.0f GeV/#it{c}",labelCentR[0][j].second.Data(),lowPt[0],lowPt[1]), StdTextSize, 0.04);
        pEtaCentLowPt[j].Plot(Form("%s/EtaComparisonRadiusLowPt_%s.pdf",outputFolder.Data(),labelCentR[0][j].second.Data()));
    }

    // Create eta distributions in different pT intervals
    const int nPtBins = 5;
    double ptBins[nPtBins+1] = {20, 40, 60, 80, 100,150}; // Define pT intervals

    // One plot per jet radius and centrality
    Plotting1D pEtaPtBins[files[0].size()][files.size()];

    // Loop over jet radii
    for (int j = 0; j < files[0].size(); j++) {
        // Loop over centralities  
        for (int i = 0; i < files.size(); i++) {
            // Loop over pT bins
            for (int iPt = 0; iPt < nPtBins; iPt++) {
                vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(ptBins[iPt], ptBins[iPt+1]);
                TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
                hEta->SetName(Form("hEta_%s_%s_pt%.0f_%.0f",labelCentR[i][j].first.Data(),
                                 labelCentR[i][j].second.Data(), ptBins[iPt], ptBins[iPt+1]));
                hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
                if(normalize) hEta->Scale(1./hEta->Integral());
                
                pEtaPtBins[j][i].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
                pEtaPtBins[j][i].SetLegend(0.15, 0.3, 0.70, 0.95, true);
                pEtaPtBins[j][i].New(hEta, Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", 
                                     ptBins[iPt], ptBins[iPt+1]), 1, 1, fancyColors[iPt], "hist");
            }
            
            pEtaPtBins[j][i].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; %s%% centrality; R = %s",
                                     labelCentR[i][j].first.Data(), labelCentR[i][j].second.Data()), 
                                     StdTextSize, 0.04);
            pEtaPtBins[j][i].Plot(Form("%s/EtaComparisonPtBins_%s_%s.pdf",outputFolder.Data(),
                                 labelCentR[i][j].second.Data(), labelCentR[i][j].first.Data()));
        }
    }

    // do the same plot but only for very low pt jets
    const int nPtBinsLowPt = 6;
    double ptBinsLowPt[nPtBinsLowPt+1] = {-10,-5,0,5,10,15,20}; // Define pT intervals
    // One plot per jet radius and centrality for low pT
    Plotting1D pEtaPtBinsLowPt[files[0].size()][files.size()];

    // Loop over jet radii
    for (int j = 0; j < files[0].size(); j++) {
        // Loop over centralities  
        for (int i = 0; i < files.size(); i++) {
            // Loop over low pT bins
            for (int iPt = 0; iPt < nPtBinsLowPt; iPt++) {
                vJetPtEtaPhi[i][j]->GetAxis(0)->SetRangeUser(ptBinsLowPt[iPt], ptBinsLowPt[iPt+1]);
                TH1F* hEta = (TH1F*)vJetPtEtaPhi[i][j]->Projection(1);
                hEta->SetName(Form("hEta_%s_%s_lowpt%.0f_%.0f",labelCentR[i][j].first.Data(),
                                 labelCentR[i][j].second.Data(), ptBinsLowPt[iPt], ptBinsLowPt[iPt+1]));
                hEta->GetXaxis()->SetRangeUser(-0.71,0.71);
                if(normalize) hEta->Scale(1./hEta->Integral());
                
                pEtaPtBinsLowPt[j][i].SetAxisLabel("Jet #bf{#eta}", "#bf{norm. counts}");
                pEtaPtBinsLowPt[j][i].SetLegend(0.15, 0.3, 0.70, 0.95, true);
                pEtaPtBinsLowPt[j][i].New(hEta, Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", 
                                     ptBinsLowPt[iPt], ptBinsLowPt[iPt+1]), 1, 1, fancyColors[iPt], "hist");
            }
            
            pEtaPtBinsLowPt[j][i].NewLatex(0.9, 0.9, Form("ALICE work in progress; Pb-Pb LHC23zzm; %s%% centrality; R = %s",
                                     labelCentR[i][j].first.Data(), labelCentR[i][j].second.Data()), 
                                     StdTextSize, 0.04);
            pEtaPtBinsLowPt[j][i].Plot(Form("%s/EtaComparisonPtBinsLowPt_%s_%s.pdf",outputFolder.Data(),
                                 labelCentR[i][j].second.Data(), labelCentR[i][j].first.Data()));
        }
    }








    
    


 }