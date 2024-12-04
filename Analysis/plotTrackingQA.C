#include "PlottingClass.h"
void plotTrackingQA(TString inputfile=" /alf/data/calorimetrees/JE_PbPb_AO2D/000006_LHC23zzm_JetPt0_NewTM_PbPb/AnalysisResults_merged.root", TString outputfolder = "../LHC23zzm_pass4/TrackingQA")
{
    // define track pt bins
    std::vector<Float_t> trackPtBins = {0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 50, 60, 80, 100, 120, 160, 200};

    TFile *input = new TFile(inputfile.Data(), "READ");

    // get TDirectory called gamma-jet-tree-producer
    TDirectory *dir = (TDirectory *)input->Get("gamma-jet-tree-producer");

    // get THnSparseF called trackPtEtaPhi 
    // axis 0: track pT axis 1: track eta axis 2: track phi
    THnSparseF *trackPtEtaPhi = (THnSparseF *)dir->Get("trackPtEtaPhi");

    // create output folder if it does not exist
    if (!std::filesystem::exists(outputfolder.Data()))
    {
        std::filesystem::create_directory(outputfolder.Data());
    }
    
    // plot 2D distribution of tack pt vs track eta
    TH2F *hTrackPtEta = (TH2F*)trackPtEtaPhi->Projection(1,0);
    hTrackPtEta->SetName("hTrackPtEta");
    Plotting2D pTrackPtEta;
    pTrackPtEta.New(hTrackPtEta, "COLZ");
    pTrackPtEta.SetAxisLabel("Track #bf{#eta}", "Track #bf{p}_{T} (GeV/#it{c})");
    pTrackPtEta.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central", StdTextSize*0.3, 0.1);
    pTrackPtEta.Plot(Form("%s/TrackPtEta.pdf",outputfolder.Data()));

    // plot 2D distribution of tack pt vs track phi
    TH2F *hTrackPtPhi = (TH2F*)trackPtEtaPhi->Projection(2,0);
    hTrackPtPhi->SetName("hTrackPtPhi");
    Plotting2D pTrackPtPhi;
    pTrackPtPhi.New(hTrackPtPhi, "COLZ");
    pTrackPtPhi.SetAxisLabel("Track #bf{#phi}", "Track #bf{p}_{T} (GeV/#it{c})");
    pTrackPtPhi.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central", StdTextSize*0.3, 0.1);
    pTrackPtPhi.Plot(Form("%s/TrackPtPhi.pdf",outputfolder.Data()));

    // PLot 1D projections of track eta in bins of track pt
    std::vector<TH1F*> hTrackEtaProjections(trackPtBins.size()-1);
    for (int i = 0; i < trackPtBins.size()-1; i++){
        int binMin = trackPtEtaPhi->GetAxis(0)->FindBin(trackPtBins[i]);
        int binMax = trackPtEtaPhi->GetAxis(0)->FindBin(trackPtBins[i+1]);
        trackPtEtaPhi->GetAxis(0)->SetRange(binMin, binMax);
        hTrackEtaProjections[i] = (TH1F*)trackPtEtaPhi->Projection(1);
        hTrackEtaProjections[i]->SetName(Form("hTrackEta_%d_%d",int(trackPtBins[i]),int(trackPtBins[i+1])));
        hTrackEtaProjections[i]->GetXaxis()->SetRangeUser(-0.9,0.9);
        hTrackEtaProjections[i]->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}",trackPtBins[i],trackPtBins[i+1]));
    }

    PlottingGrid pTrackEtaProjections;
    pTrackEtaProjections.SetMargins(0.12, 0.12, 0.05, 0.125);
    pTrackEtaProjections.SetAxisLabel("Track #bf{#eta}", "Counts");
    for (int i = 0; i < hTrackEtaProjections.size(); i++){
        pTrackEtaProjections.New(hTrackEtaProjections[i],"",-1,1,kBlack,"hist");
        pTrackEtaProjections.NextPad(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}",trackPtBins[i],trackPtBins[i+1]));
    }
    pTrackEtaProjections.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central; #it{E}_{clus} > 8 GeV", StdTextSize*0.3, 0.1);
    pTrackEtaProjections.Plot(Form("%s/TrackEtaProjections.pdf",outputfolder.Data()));
    // ignore
    

    // plot 2D histogram eta phi
    TH2F *hTrackEtaPhi = (TH2F*)trackPtEtaPhi->Projection(1,2);
    hTrackEtaPhi->SetName("hTrackEtaPhi");
    Plotting2D pTrackEtaPhi;
    pTrackEtaPhi.New(hTrackEtaPhi, "COLZ");
    pTrackEtaPhi.SetAxisLabel("Track #bf{#phi}", "Track #bf{#eta}");
    pTrackEtaPhi.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central", StdTextSize*0.3, 0.1);
    pTrackEtaPhi.Plot(Form("%s/TrackEtaPhi.pdf",outputfolder.Data()));

    // plot phi vs pt
    TH2F *hTrackPhiPt = (TH2F*)trackPtEtaPhi->Projection(2,0);
    hTrackPhiPt->SetName("hTrackPhiPt");
    Plotting2D pTrackPhiPt;
    pTrackPhiPt.New(hTrackPhiPt, "COLZ");
    pTrackPhiPt.SetAxisLabel("Track #bf{#phi}", "Track #bf{p}_{T} (GeV/#it{c})");
    pTrackPhiPt.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central", StdTextSize*0.3, 0.1);
    pTrackPhiPt.Plot(Form("%s/TrackPhiPt.pdf",outputfolder.Data()));

    // Plot 1D projections of track phi in bins of track pt
    std::vector<TH1F*> hTrackPhiProjections(trackPtBins.size()-1);
    for (int i = 0; i < trackPtBins.size()-1; i++){
        int binMin = trackPtEtaPhi->GetAxis(0)->FindBin(trackPtBins[i]);
        int binMax = trackPtEtaPhi->GetAxis(0)->FindBin(trackPtBins[i+1]);
        trackPtEtaPhi->GetAxis(0)->SetRange(binMin, binMax);
        hTrackPhiProjections[i] = (TH1F*)trackPtEtaPhi->Projection(2);
        hTrackPhiProjections[i]->SetName(Form("hTrackPhi_%d_%d",int(trackPtBins[i]),int(trackPtBins[i+1])));
        hTrackPhiProjections[i]->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
        hTrackPhiProjections[i]->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}",trackPtBins[i],trackPtBins[i+1]));
    }

    PlottingGrid pTrackPhiProjections;
    pTrackPhiProjections.SetMargins(0.12, 0.12, 0.05, 0.125);
    pTrackPhiProjections.SetAxisLabel("Track #bf{#phi}", "Counts");
    for (int i = 0; i < hTrackPhiProjections.size(); i++){
        pTrackPhiProjections.New(hTrackPhiProjections[i],"",-1,1,kBlack,"hist");
        pTrackPhiProjections.NextPad(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}",trackPtBins[i],trackPtBins[i+1]));
    }
    pTrackPhiProjections.NewLatex(0.9, 0.9, "#bf{ALICE work in progress}; Pb-Pb LHC23zzm; 0-100% central; #it{E}_{clus} > 8 GeV", StdTextSize*0.3, 0.1);
    pTrackPhiProjections.Plot(Form("%s/TrackPhiProjections.pdf",outputfolder.Data()));


   

}