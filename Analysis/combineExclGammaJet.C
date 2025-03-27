#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "yaml-cpp/yaml.h"
#include "analyseExclGammaJet.h"

Int_t fancyColors[8];
TString fancyColorHex[8] = {"#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff"};

TGraphAsymmErrors* GetRecoilJetPtIntegralPhi(CombineExclGammaJetOptions options, CentralityEnum centrality, JetRadiusEnum jetRadius, double recoilJetPtMin, double recoilJetPtMax,bool isMergedPi0s = true){
    if (isMergedPi0s){
        return (TGraphAsymmErrors*) options.GetRootFile(centrality, jetRadius)->Get(Form("RecoilJet_MergedPi0s/gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtMin, recoilJetPtMax));
    } else {
        return (TGraphAsymmErrors*) options.GetRootFile(centrality, jetRadius)->Get(Form("RecoilJet_Photons/gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtMin, recoilJetPtMax));
    }
}

void drawDeltaPhiDistributions(CombineExclGammaJetOptions options, bool mergedPi0 = true){
    // loop over all jet radius folders and centrality folders
    std::vector<CentralityEnum> centralities = {k0_10, k10_30, k30_50, k50_90};
    std::vector<TString> centralityLabels = {"0-10%", "10-30%", "30-50%", "50-90%"};
    std::vector<JetRadiusEnum> jetRadii = {kR02, kR03, kR04, kR05, kR06};
    std::vector<TString> jetRadiusLabels = {"R=0.2", "R=0.3", "R=0.4", "R=0.5", "R=0.6"};
    std::vector<TString> jetRadiusFileStrings = {"R02", "R03", "R04", "R05", "R06"};
    Plotting1D PRecoilPhi[options.jetRadiusFolders.size()][recoilJetPtBinsIntegration.size()];
    int iJetRadius = 0;
   
    for (auto jetRadius : jetRadii){
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        ExclusiveTriggerParticleSelection exclusiveTriggerPhotonSelection(varOptions);
        double sigTriggerPtMin = exclusiveTriggerPhotonSelection.getPtMinSignal();
        double sigTriggerPtMax = exclusiveTriggerPhotonSelection.getPtMaxSignal();
        double refTriggerPtMin = exclusiveTriggerPhotonSelection.getPtMinReference();
        double refTriggerPtMax = exclusiveTriggerPhotonSelection.getPtMaxReference();
        TString triggerLabel = mergedPi0 ? Form("merged #pi^{0} trigger TT{%d,%d} - TT{%d,%d}", (int)sigTriggerPtMin, (int)sigTriggerPtMax, (int)refTriggerPtMin, (int)refTriggerPtMax) : Form("iso. photon trigger TT{%d,%d} - TT{%d,%d}", (int)sigTriggerPtMin, (int)sigTriggerPtMax, (int)refTriggerPtMin, (int)refTriggerPtMax);
        int iRecoilJetPtBin = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration){
            // PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetMargin 

            // loop over centralities
            double maxY = 0;
            for (auto centrality : centralities){
                TGraphAsymmErrors* histo = GetRecoilJetPtIntegralPhi(options, centrality, jetRadius, recoilJetPtBin.first, recoilJetPtBin.second, mergedPi0);
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].New(histo, centralityLabels[centrality].Data(), 20, 2, fancyColors[centrality], "pe");
                // find maximum y value in TGraphAsymmErrors
                for (int i = 0; i < histo->GetN(); i++){
                    if (histo->GetY()[i] > maxY){
                        maxY = histo->GetY()[i];
                    }
                }
            }
            
            // Set plot labels and ranges
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetLegend(0.25, 0.5, 0.4, 0.6, true);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisRange(0, TMath::Pi(), 0, maxY*1.1);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T},%s;%d < p_{T,jet}^{ch} < %d GeV/c", options.dataset.Data(), triggerLabel.Data(), 
                                                                            jetRadiusLabels[iJetRadius].Data(), (int)recoilJetPtBin.first, (int)recoilJetPtBin.second));
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisLabel("#Delta#varphi", "1/N_{trig} dN/d#Delta#varphi", 1, 1.4);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisRange(1.38, 1.05*TMath::Pi(), 0, maxY*1.15);

            // Draw dashed line at pi

            
            // Create output directory if it doesn't exist
            TString outputPath = Form("%s/ExclGammaJetCombined/DeltaPhiDistributions/%s", options.dataset.Data(), jetRadiusFileStrings[iJetRadius].Data());
            gSystem->mkdir(outputPath, kTRUE);
            
            // Plot and save
            if (mergedPi0){
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/RecoilJetPhi_mergedPi0s_%s_%d_%d.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[iJetRadius].Data(),
                                                             (int)recoilJetPtBin.first, 
                                                             (int)recoilJetPtBin.second), 
                                                        false, false);
            } else {
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/RecoilJetPhi_isoPhoton_%s_%d_%d.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[iJetRadius].Data(),
                                                             (int)recoilJetPtBin.first, 
                                                             (int)recoilJetPtBin.second), 
                                                        false, false);
            }
            iRecoilJetPtBin++;
        }
        iJetRadius++;
}
}
void combineExclGammaJet(TString AnalysisDirectory="", TString configYaml = "RunConfigFlorian.yaml"){
    std::vector<TString> folders;
    for (int i = 0; i < 8; i++){
        fancyColors[i] = TColor::GetColor(fancyColorHex[i]);
    }
    CombineExclGammaJetOptions options(AnalysisDirectory, configYaml);
    options.printOptions();
    
    // find a specific file
    cout << "File path: " << options.GetInputFilePath(k0_10, kR02) << endl;
    cout << "Root file: " << options.GetRootFile(k0_10, kR02) << endl;

    drawDeltaPhiDistributions(options, true);
    drawDeltaPhiDistributions(options, false);


    
}