#include "PlottingClass.h"
#include <sys/stat.h>
#include <sys/types.h>
#include "Utilities.h"
#include "Cuts.h"
#include "yaml-cpp/yaml.h"
#include "analyseExclGammaJet.h"

Int_t fancyColors[8];
TString fancyColorHex[8] = {"#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff"};



template <typename T>
T* GetHistogram(const char* histoName, CombineExclGammaJetOptions options, CentralityEnum centrality, JetRadiusEnum jetRadius, bool isMergedPi0s = true) {
    TString folderName = isMergedPi0s ? "RecoilJet_MergedPi0s" : "RecoilJet_Photons";
    TString fullPath = Form("%s/%s", folderName.Data(), histoName);
    
    T* histo = (T*) options.GetRootFile(centrality, jetRadius)->Get(fullPath.Data());
    if (!histo) {
        std::cerr << "Warning: Could not find histogram " << fullPath.Data() 
                  << " in file for centrality " << centrality 
                  << " and jet radius " << jetRadius << std::endl;
        return nullptr;
    }
    
    return histo;
}
template <typename T>
T* GetHistogramPPRef(const char* histoName, CombineExclGammaJetOptions options, JetRadiusEnum jetRadius, bool isMergedPi0s = true) {
    TString folderName = isMergedPi0s ? "RecoilJet_MergedPi0s" : "RecoilJet_Photons";
    TString fullPath = Form("%s/%s", folderName.Data(), histoName);
    
    T* histo = (T*) options.GetRootFile(jetRadius)->Get(fullPath.Data());
    if (!histo) {
        std::cerr << "Warning: Could not find pp ref histogram " << fullPath.Data() 
                  << " in file for jet radius " << jetRadius << std::endl;
        return nullptr;
    }
    
    return histo;
}


// TODO actually make it treat asymmetric errors. For now everything assumed to be symmetric
TGraphAsymmErrors* CalculateRatioTGraphAsymmErrors(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TString name){
    // check if graphs exist
    if (!g1 || !g2) {
        std::cerr << "Error: One or both graphs are nullptr" << std::endl;
        return nullptr;
    }

    Double_t* x1 = g1->GetX();
    Double_t* ex1 = g1->GetEXlow();
    Double_t* y1 = g1->GetY();
    Double_t* ey1 = g1->GetEYlow();

    Double_t* x2 = g2->GetX();
    Double_t* ex2 = g2->GetEXlow();
    Double_t* y2 = g2->GetY();
    Double_t* ey2 = g2->GetEYlow();

    // check that same number of points and x values are the same
    if (g1->GetN() != g2->GetN()){
        std::cerr << "Error: Different number of points in graphs" << std::endl;
        return nullptr;
    }   
    for (int i = 0; i < g1->GetN(); i++){
        if (x1[i] != x2[i]){
            std::cerr << "Error: Different x values in graphs" << std::endl;
            return nullptr;
        }
    }

    double ratio[g1->GetN()];
    double ratioErr[g1->GetN()];
    for (int i = 0; i < g1->GetN(); i++){
        // make sure y1 and y2 are not zero, if they are, set the point and error to zero
        if(y1[i] == 0 || y2[i] == 0){
            ratio[i] = 0;
            ratioErr[i] = 0;
        } else {
            ratio[i] = y1[i] / y2[i];
            ratioErr[i] = ratio[i] * TMath::Sqrt(TMath::Power(ey1[i]/y1[i], 2) + TMath::Power(ey2[i]/y2[i], 2));
        }
    }


    TGraphAsymmErrors* ratioGraph = new TGraphAsymmErrors(g1->GetN(), x1, ratio, ex1,ex1, ratioErr,ratioErr);
    ratioGraph->SetName(name.Data());
    return ratioGraph;
}

TString getTriggerLabel(GlobalOptions options, bool mergedPi0 = true){
    ExclusiveTriggerParticleSelection exclusiveTriggerPhotonSelection(options);
    double sigTriggerPtMin = exclusiveTriggerPhotonSelection.getPtMinSignal();
    double sigTriggerPtMax = exclusiveTriggerPhotonSelection.getPtMaxSignal();
    double refTriggerPtMin = exclusiveTriggerPhotonSelection.getPtMinReference();
    double refTriggerPtMax = exclusiveTriggerPhotonSelection.getPtMaxReference();
    return Form("merged #pi^{0} trigger TT{%d,%d} - TT{%d,%d}", (int)sigTriggerPtMin, (int)sigTriggerPtMax, (int)refTriggerPtMin, (int)refTriggerPtMax);
}

void drawCRefPhiDistributions(CombineExclGammaJetOptions options, bool mergedPi0 = true){
    Plotting1D PCRefPhiDistributions[options.jetRadiusFolders.size()];
    int iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum){
        double maxY = 0;
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        int iCentrality = 0;
        for (auto centrality : options.centralityFoldersEnum){
            TGraphAsymmErrors* gCRefValues = GetHistogram<TGraphAsymmErrors>("gCRefValuesVsPhi", options, centrality, jetRadius, mergedPi0);
            if (!gCRefValues) continue;
            
            PCRefPhiDistributions[iJetRadius].New(gCRefValues, centralityLabels[centrality].Data(), 20, 2, fancyColors[iCentrality], "pe");
            // find maximum y value in TGraphAsymmErrors
            for (int i = 0; i < gCRefValues->GetN(); i++){
                if (gCRefValues->GetY()[i] > maxY){
                    maxY = gCRefValues->GetY()[i];
                }
            }
            iCentrality++;
        }
        // Create output directory if it doesn't exist
        TString outputPath = Form("%s/ExclGammaJetCombined/CRefPhiDistributions/%s", options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data());
        gSystem->mkdir(outputPath, kTRUE);
        // Set plot labels and ranges
        PCRefPhiDistributions[iJetRadius].SetLegend(0.25, 0.5, 0.4, 0.6, true);
        PCRefPhiDistributions[iJetRadius].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        PCRefPhiDistributions[iJetRadius].SetAxisRange(1.38, TMath::Pi(), 0.9, 1.1);
        // draw dashed line at unity
        PCRefPhiDistributions[iJetRadius].NewLine(1.38, 1., TMath::Pi(), 1., 2, 2, kGray);
        PCRefPhiDistributions[iJetRadius].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T},%s", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data()));
        PCRefPhiDistributions[iJetRadius].SetAxisLabel("#varphi", "c_{ref}");
        PCRefPhiDistributions[iJetRadius].Plot(Form("%s/CRefPhiDistributions_%s_%s_%s.pdf", outputPath.Data(), options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data(), mergedPi0 ? "mergedPi0s" : "isoPhoton"), false, false);
        
        iJetRadius++;
    }

    // plot the slope of pol1 fit to cRef vs phi
    Plotting1D PCRefSlopePhiDistributions[options.jetRadiusFolders.size()];
    iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum){
        double maxY = -99999;
        double minY = 99999;
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        int iCentrality = 0;
        for (auto centrality : options.centralityFoldersEnum){
            TGraphAsymmErrors* gCRefSlopeValues = GetHistogram<TGraphAsymmErrors>("gCRefSlopeValuesVsPhi", options, centrality, jetRadius, mergedPi0);
            if (!gCRefSlopeValues) continue;
            
            PCRefSlopePhiDistributions[iJetRadius].New(gCRefSlopeValues, centralityLabels[centrality].Data(), 20, 2, fancyColors[iCentrality], "pe");
            // find maximum and minimum y values in TGraphAsymmErrors
            for (int i = 0; i < gCRefSlopeValues->GetN(); i++){
                if (gCRefSlopeValues->GetY()[i] > maxY){
                    maxY = gCRefSlopeValues->GetY()[i];
                }
                if (gCRefSlopeValues->GetY()[i] < minY){
                    minY = gCRefSlopeValues->GetY()[i];
                }
            }
            iCentrality++;
        }
        // Create output directory if it doesn't exist
        TString outputPath = Form("%s/ExclGammaJetCombined/CRefSlopePhiDistributions/%s", options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data());
        gSystem->mkdir(outputPath, kTRUE);
        // Set plot labels and ranges
        PCRefSlopePhiDistributions[iJetRadius].SetLegend(0.25, 0.5, 0.4, 0.6, true);
        PCRefSlopePhiDistributions[iJetRadius].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        PCRefSlopePhiDistributions[iJetRadius].SetAxisRange(1.38, TMath::Pi(), minY*1.1, maxY*1.1);
        // draw dashed line at zero
        PCRefSlopePhiDistributions[iJetRadius].NewLine(1.38, 0., TMath::Pi(), 0., 2, 2, kGray);
        PCRefSlopePhiDistributions[iJetRadius].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T},%s", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data()));
        PCRefSlopePhiDistributions[iJetRadius].SetAxisLabel("#varphi", "c_{ref} slope");
        PCRefSlopePhiDistributions[iJetRadius].Plot(Form("%s/CRefSlopePhiDistributions_%s_%s_%s.pdf", outputPath.Data(), options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data(), mergedPi0 ? "mergedPi0s" : "isoPhoton"), false, false);
        
        iJetRadius++;
    }
    
}


void drawDeltaPhiDistributions(CombineExclGammaJetOptions options, bool mergedPi0 = true){
    // loop over all jet radius folders and centrality folders
    Plotting1D PRecoilPhi[options.jetRadiusFolders.size()][recoilJetPtBinsIntegration.size()];
    int iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum){
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        int iRecoilJetPtBin = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration){
            // PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetMargin 

            // loop over options.centralityFoldersEnum
            double maxY = 0;
            for (auto centrality : options.centralityFoldersEnum){
                // TGraphAsymmErrors* histo = GetRecoilJetPtIntegralPhi(options, centrality, jetRadius, recoilJetPtBin.first, recoilJetPtBin.second, mergedPi0);
                // get histogram from root file using GetHistogram function instead of GetRecoilJetPtIntegralPhi
                TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
                TGraphAsymmErrors* histo = GetHistogram<TGraphAsymmErrors>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].New(histo, centralityLabels[centrality].Data(), 20, 2, fancyColors[centrality], "pe");
                // find maximum y value in TGraphAsymmErrors
                for (int i = 0; i < histo->GetN(); i++){
                    if (histo->GetY()[i] > maxY){
                        maxY = histo->GetY()[i];
                    }
                }
            }

            if(options.ppRefAvailable){
                // plot ratio to pp reference
                TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
                TGraphAsymmErrors* refHisto = GetHistogramPPRef<TGraphAsymmErrors>(histoName.Data(), options, jetRadius, mergedPi0);
                if (!refHisto) continue;
                refHisto->SetName(Form("%s_%f_%f_ref", histoName.Data(), recoilJetPtBin.first, recoilJetPtBin.second));
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].New(refHisto, "pp", 20, 2, kBlack, "pe");

                for (int i = 0; i <= refHisto->GetN(); i++){
                    if (refHisto->GetY()[i] > maxY){
                        maxY = refHisto->GetY()[i];
                    }
                }
            }
            
            // Set plot labels and ranges
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetLegend(0.25, 0.5, 0.4, 0.6, true);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisRange(0, TMath::Pi(), 0, maxY*1.1);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T},%s;%d < p_{T,jet}^{ch} < %d GeV/c", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data(), (int)recoilJetPtBin.first, (int)recoilJetPtBin.second));
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisLabel("#Delta#varphi", "1/N_{trig} dN/d#Delta#varphi", 1, 1.4);
            PRecoilPhi[iJetRadius][iRecoilJetPtBin].SetAxisRange(1.38, 1.05*TMath::Pi(), 0, maxY*1.15);

            // Draw dashed line at pi

            
            // Create output directory if it doesn't exist
            TString outputPath = Form("%s/ExclGammaJetCombined/DeltaPhiDistributions/%s", options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data());
            gSystem->mkdir(outputPath, kTRUE);
            
            // Plot and save
            if (mergedPi0){
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/RecoilJetPhi_mergedPi0s_%s_%d_%d.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[jetRadius].Data(),
                                                             (int)recoilJetPtBin.first, 
                                                             (int)recoilJetPtBin.second), 
                                                        false, false);

                // also plot with log y axis
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/RecoilJetPhi_mergedPi0s_%s_%d_%d_log.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[jetRadius].Data(),
                                                             (int)recoilJetPtBin.first, 
                                                             (int)recoilJetPtBin.second), 
                                                        false, true);
            } else {
                PRecoilPhi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/RecoilJetPhi_isoPhoton_%s_%d_%d.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[jetRadius].Data(),
                                                             (int)recoilJetPtBin.first, 
                                                             (int)recoilJetPtBin.second), 
                                                        false, false);
            }
            iRecoilJetPtBin++;
        }
        iJetRadius++;
    }

    // 
}

void drawRecoilJetPtDistributions(CombineExclGammaJetOptions options, JetRadiusEnum jetRefRadius, bool mergedPi0 = true){
    // the histograms are hRecoilJetPtSignal_Subtracted_%.2f_%.2f", phiBin.first, phiBin.second and we want the one that is back to back which is the first one of phiBins
    Plotting1D PRecoilJetPt[options.jetRadiusFolders.size()][phiBins.size()];
    int iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum){
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        int iRecoilPhiBin = 0;
        for (auto recoilPhiBin : phiBins){
            double maxY = 0;
            for (auto centrality : options.centralityFoldersEnum){
                TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", recoilPhiBin.first, recoilPhiBin.second);
                TH1F* histo = GetHistogram<TH1F>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
                if (!histo) continue;
                PRecoilJetPt[iJetRadius][iRecoilPhiBin].New(histo, centralityLabels[centrality].Data(), 20, 2, fancyColors[centrality], "pe");
                // find maximum y value in TH1F
                for (int i = 0; i < histo->GetNbinsX(); i++){
                    if (histo->GetBinContent(i) > maxY){
                        maxY = histo->GetBinContent(i);
                    }
                }
            }

            // Set plot labels and ranges
            PRecoilJetPt[iJetRadius][iRecoilPhiBin].SetLegend(0.25, 0.5, 0.4, 0.6, true);
            PRecoilJetPt[iJetRadius][iRecoilPhiBin].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
            PRecoilJetPt[iJetRadius][iRecoilPhiBin].SetAxisRange(0, 150, 0, maxY*1.1);
            TString phiBinString = Form("%.2f < #varphi < %.2f", (float)recoilPhiBin.first, (float)recoilPhiBin.second);
            if(iRecoilPhiBin == 0)
                phiBinString = "|#Delta#varphi-#pi| < 0.6";
            PRecoilJetPt[iJetRadius][iRecoilPhiBin].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T},%s;%s", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data(), phiBinString.Data()));
            PRecoilJetPt[iJetRadius][iRecoilPhiBin].SetAxisLabel("p_{T,jet}^{ch} (GeV/c)", "#Delta_{recoil} p_{T,jet}^{ch}");

            TString outputPath = Form("%s/ExclGammaJetCombined/RecoilJetPtDistributions/%s", options.dataset.Data(), jetRadiusFileStrings[jetRadius].Data());
            gSystem->mkdir(outputPath, kTRUE);

            if(mergedPi0){
                PRecoilJetPt[iJetRadius][iRecoilPhiBin].Plot(Form("%s/RecoilJetPt_mergedPi0s_%s_%.2f_%.2f.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[jetRadius].Data(),
                                                             (float)recoilPhiBin.first, 
                                                             (float)recoilPhiBin.second), 
                                                        false, true);
            } else {
                PRecoilJetPt[iJetRadius][iRecoilPhiBin].Plot(Form("%s/RecoilJetPt_isoPhoton_%s_%.2f_%.2f.pdf", 
                                                             outputPath.Data(), 
                                                             jetRadiusFileStrings[jetRadius].Data(),
                                                             (float)recoilPhiBin.first, 
                                                             (float)recoilPhiBin.second), 
                                                        false, true);
            }
            
            iRecoilPhiBin++;
        }  
        iJetRadius++;
    }

    // Draw the ratio of all radii with respect to R=0.4 for different options.centralityFoldersEnum
    // Draw ratio plots for different radii relative to R=0.4
    Plotting1D PRecoilJetPtRatio[options.jetRadiusFoldersEnum.size()];
    
    // Only look at back-to-back bin (first phi bin)
    auto backToBackBin = phiBins[0];
    
    // Get R=0.4 histograms for each centrality as reference
    std::map<TString, TH1F*> referenceHistos;
    for (auto centrality : options.centralityFoldersEnum) {
        TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
        TH1F* refHisto = GetHistogram<TH1F>(histoName.Data(), options, centrality, jetRefRadius, mergedPi0);
        if (!refHisto) continue;
        
        // Rebin reference histogram using coarse bins
        TH1F* refHistoRebin = (TH1F*)refHisto->Rebin(nCoarseJetPtBins, Form("%s_rebin", refHisto->GetName()), coarseJetPtBins);
        referenceHistos[centrality] = refHistoRebin;
    }

    // Loop over jet radii (skip R=0.4 since it's the reference)
    iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        if (jetRadius == jetRefRadius) {
            iJetRadius++;
            continue;
        }
        
        double maxY = 0;
        // Loop over options.centralityFoldersEnum
        for (auto centrality : options.centralityFoldersEnum) {
            TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
            TH1F* histo = GetHistogram<TH1F>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
            if (!histo) continue;
            
            // Rebin histogram
            TH1F* histoRebin = (TH1F*)histo->Rebin(nCoarseJetPtBins, Form("%s_rebin", histo->GetName()), coarseJetPtBins);
            
            // Calculate ratio to R=0.4
            TH1F* ratioHisto = (TH1F*)histoRebin->Clone(Form("%s_ratio", histoRebin->GetName()));
            ratioHisto->Divide(referenceHistos[centrality]);
            
            PRecoilJetPtRatio[iJetRadius].New(ratioHisto, centralityLabels[centrality].Data(), 20, 2, fancyColors[centrality], "pe");
            
            // Find maximum y value
            for (int i = 1; i <= ratioHisto->GetNbinsX(); i++) {
                if (ratioHisto->GetBinContent(i) > maxY) {
                    maxY = ratioHisto->GetBinContent(i);
                }
            }
        }

        // Set plot labels and ranges
        PRecoilJetPtRatio[iJetRadius].SetLegend(0.35, 0.55, 0.5, 0.65, true);
        PRecoilJetPtRatio[iJetRadius].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        PRecoilJetPtRatio[iJetRadius].SetAxisRange(0, 150, 0, 4.);
        PRecoilJetPtRatio[iJetRadius].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T};|#Delta#varphi-#pi| < 0.6", options.dataset.Data(), triggerLabel.Data()));
        PRecoilJetPtRatio[iJetRadius].SetAxisLabel("p_{T,jet}^{ch} (GeV/c)", Form("#Delta_{recoil}(%s)/#Delta_{recoil}(%s)", jetRadiusLabels[jetRadius].Data(), jetRadiusLabels[jetRefRadius].Data()));

        // draw dashed line at 1
        PRecoilJetPtRatio[iJetRadius].NewLine(0, 1., 150, 1., 2, 2, kGray);

        TString outputPath = Form("%s/ExclGammaJetCombined/RecoilJetPtRatios", options.dataset.Data());
        gSystem->mkdir(outputPath, kTRUE);

        if(mergedPi0) {
            PRecoilJetPtRatio[iJetRadius].Plot(Form("%s/RecoilJetPtRadiusRatio_mergedPi0s_%s_Ref_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), jetRadiusFileStrings[jetRefRadius].Data()), false, false);
        } else {
            PRecoilJetPtRatio[iJetRadius].Plot(Form("%s/RecoilJetPtRadiusRatio_isoPhoton_%s_Ref_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), jetRadiusFileStrings[jetRefRadius].Data()), false, false);
        }
        
        iJetRadius++;
    }
    
    
}

void drawICPDistributions(CombineExclGammaJetOptions options, CentralityEnum refCentrality,  bool mergedPi0 = true){
    
     // Draw the ratio of all radii with respect to R=0.4 for different options.centralityFoldersEnum
    // Draw ratio plots for different radii relative to R=0.4
    Plotting1D PICP_Pt[options.jetRadiusFoldersEnum.size()];
    
    // Only look at back-to-back bin (first phi bin)
    auto backToBackBin = phiBins[0];
    
    // Get reference histogram for centrality refCentrality for all radii
    std::map<TString, TH1F*> referenceHistos;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
        TH1F* refHisto = GetHistogram<TH1F>(histoName.Data(), options, refCentrality, jetRadius, mergedPi0);
        if (!refHisto) continue;
        
        // Rebin reference histogram using coarse bins
        TH1F* refHistoRebin = (TH1F*)refHisto->Rebin(nCoarseJetPtBins, Form("%s_rebin", refHisto->GetName()), coarseJetPtBins);
        referenceHistos[jetRadius] = refHistoRebin;
    }

    // Loop over jet radii (skip R=0.4 since it's the reference)
    int iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        
        double maxY = 0;
        // Loop over options.centralityFoldersEnum
        for (auto centrality : options.centralityFoldersEnum) {
            // skip reference centrality
            if (centrality >= refCentrality) continue;

            TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
            TH1F* histo = GetHistogram<TH1F>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
            if (!histo) continue;
            
            // Rebin histogram
            TH1F* histoRebin = (TH1F*)histo->Rebin(nCoarseJetPtBins, Form("%s_rebin", histo->GetName()), coarseJetPtBins);
            
            // Calculate ratio to peripheral
            TH1F* ratioHisto = (TH1F*)histoRebin->Clone(Form("%s_ratio", histoRebin->GetName()));
            ratioHisto->Divide(referenceHistos[jetRadius]);
            
            PICP_Pt[iJetRadius].New(ratioHisto, Form("%s/%s", centralityLabels[centrality].Data(), centralityLabels[refCentrality].Data()), 20, 2, fancyColors[centrality], "pe");
            
            // Find maximum y value
            for (int i = 1; i <= ratioHisto->GetNbinsX(); i++) {
                if (ratioHisto->GetBinContent(i) > maxY) {
                    maxY = ratioHisto->GetBinContent(i);
                }
            }
        }

        // Set plot labels and ranges
        PICP_Pt[iJetRadius].SetLegend(0.25, 0.45, 0.2, 0.35, true);
        PICP_Pt[iJetRadius].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        PICP_Pt[iJetRadius].SetAxisRange(0, 205, 0, 1.8);
        PICP_Pt[iJetRadius].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T} %s;|#Delta#varphi-#pi| < 0.6", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data()));
        PICP_Pt[iJetRadius].SetAxisLabel("p_{T,jet}^{ch} (GeV/c)", Form("I_{CP} (X / %s)",centralityLabels[refCentrality].Data()));

        // draw dashed line at 1
        PICP_Pt[iJetRadius].NewLine(0, 1., 205, 1., 2, 2, kGray);

        TString outputPath = Form("%s/ExclGammaJetCombined/ICP_Pt", options.dataset.Data());
        gSystem->mkdir(outputPath, kTRUE);

        if(mergedPi0) {
            PICP_Pt[iJetRadius].Plot(Form("%s/ICP_Pt_mergedPi0s_%s_Ref_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), centralityLabels[refCentrality].Data()), false, false);
        } else {
            PICP_Pt[iJetRadius].Plot(Form("%s/ICP_Pt_isoPhoton_%s_Ref_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), centralityLabels[refCentrality].Data()), false, false);
        }
        
        iJetRadius++;
    }
    // do the same figure but as a function of delta phi
    Plotting1D PICP_Phi[options.jetRadiusFoldersEnum.size()][recoilJetPtBinsIntegration.size()]; 
    iJetRadius = 0;
    std::map<std::pair<int, int>, TGraphAsymmErrors*> referenceHistosPhi;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);

        // reference histogram for this radius and each jet pt
        int irpt = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration) {
            TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
            TGraphAsymmErrors* refHisto = GetHistogram<TGraphAsymmErrors>(histoName.Data(), options, refCentrality, jetRadius, mergedPi0);
            if (!refHisto) continue;
            refHisto->SetName(Form("%s_%f_%f_ref", histoName.Data(), recoilJetPtBin.first, recoilJetPtBin.second));
            referenceHistosPhi[std::make_pair(iJetRadius, irpt)] = refHisto;
            irpt++;
        }
        // loop over jetintegration pt
        int iRecoilJetPtBin = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration) {
            // loop over options.centralityFoldersEnum
            double maxY = 0;
            for (auto centrality : options.centralityFoldersEnum) {
                // skip reference centrality
                if (centrality >= refCentrality) continue;

                TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
                TGraphAsymmErrors* histo = GetHistogram<TGraphAsymmErrors>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
                if (!histo) continue;

                // calculate ratio to peripheral

                TGraphAsymmErrors* ratioHisto = CalculateRatioTGraphAsymmErrors(histo, referenceHistosPhi[std::make_pair(iJetRadius, iRecoilJetPtBin)], Form("%s_%f_%f_ratio", histo->GetName(), recoilJetPtBin.first, recoilJetPtBin.second));
                PICP_Phi[iJetRadius][iRecoilJetPtBin].New(ratioHisto, Form("%s/%s", centralityLabels[centrality].Data(), centralityLabels[refCentrality].Data()), 20, 2, fancyColors[centrality], "pe");

                for (int i = 0; i < ratioHisto->GetN(); i++){
                    if (ratioHisto->GetY()[i] > maxY){
                        maxY = ratioHisto->GetY()[i];
                    }
                }   
            }

            // Set plot labels and ranges
            PICP_Phi[iJetRadius][iRecoilJetPtBin].SetLegend(0.25, 0.45, 0.2, 0.35, true);
            PICP_Phi[iJetRadius][iRecoilJetPtBin].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
            PICP_Phi[iJetRadius][iRecoilJetPtBin].SetAxisRange(0, TMath::Pi(), 0, 1.8);
            PICP_Phi[iJetRadius][iRecoilJetPtBin].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T} %s;%.1f < p_{T,jet}^{ch} < %.1f GeV/c;|#Delta#varphi-#pi| < 0.6", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second));
            PICP_Phi[iJetRadius][iRecoilJetPtBin].SetAxisLabel("#Delta#varphi", Form("I_{CP} (X / %s)",centralityLabels[refCentrality].Data()));

            // draw dashed line at 1
            PICP_Phi[iJetRadius][iRecoilJetPtBin].NewLine(0, 1., TMath::Pi(), 1., 2, 2, kGray);

            TString outputPath = Form("%s/ExclGammaJetCombined/ICP_Phi", options.dataset.Data());
            gSystem->mkdir(outputPath, kTRUE);  

            if(mergedPi0) {
                PICP_Phi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/ICP_Phi_mergedPi0s_%s_Ref_%s_%.1f_%.1f.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), centralityLabels[refCentrality].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second), false, false);
            } else {
                PICP_Phi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/ICP_Phi_isoPhoton_%s_Ref_%s_%.1f_%.1f.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), centralityLabels[refCentrality].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second), false, false);
            }
            
            iRecoilJetPtBin++;
        }
        iJetRadius++;
    }
        
}

// draw IAA distributions
void drawIAAdistributions(CombineExclGammaJetOptions options, bool mergedPi0 = true){

    Plotting1D IAA_Pt[options.jetRadiusFoldersEnum.size()];

    auto backToBackBin = phiBins[0];

    // get reference histogram for all radii
    std::map<TString, TH1F*> referenceHistos;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
        TH1F* refHisto = GetHistogramPPRef<TH1F>(histoName.Data(), options, jetRadius, mergedPi0);
        if (!refHisto) continue;

        // rebin reference histogram using coarse bins
        TH1F* refHistoRebin = (TH1F*)refHisto->Rebin(nCoarseJetPtBins, Form("%s_rebin", refHisto->GetName()), coarseJetPtBins);
        referenceHistos[jetRadius] = refHistoRebin;
    }

    int iJetRadius = 0;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);
        
        double maxY = 0;
        // Loop over options.centralityFoldersEnum
        for (auto centrality : options.centralityFoldersEnum) {
            TString histoName = Form("hRecoilJetPtSignal_Subtracted_%.2f_%.2f", backToBackBin.first, backToBackBin.second);
            TH1F* histo = GetHistogram<TH1F>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
            if (!histo) continue;
            
            // Rebin histogram
            TH1F* histoRebin = (TH1F*)histo->Rebin(nCoarseJetPtBins, Form("%s_rebin", histo->GetName()), coarseJetPtBins);
            
            // Calculate ratio to peripheral
            TH1F* ratioHisto = (TH1F*)histoRebin->Clone(Form("%s_ratio", histoRebin->GetName()));
            ratioHisto->Divide(referenceHistos[jetRadius]);
            
            IAA_Pt[iJetRadius].New(ratioHisto, Form("%s/pp", centralityLabels[centrality].Data()), 20, 2, fancyColors[centrality], "pe");
            
            // Find maximum y value
            for (int i = 1; i <= ratioHisto->GetNbinsX(); i++) {
                if (ratioHisto->GetBinContent(i) > maxY) {
                    maxY = ratioHisto->GetBinContent(i);
                }
            }
        }

        // Set plot labels and ranges
        IAA_Pt[iJetRadius].SetLegend(0.25, 0.45, 0.2, 0.35, true);
        IAA_Pt[iJetRadius].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
        IAA_Pt[iJetRadius].SetAxisRange(0, 205, 0, 1.8);
        IAA_Pt[iJetRadius].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T} %s;|#Delta#varphi-#pi| < 0.6", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data()));
        IAA_Pt[iJetRadius].SetAxisLabel("p_{T,jet}^{ch} (GeV/c)", "I_{AA}");

        // draw dashed line at 1
        IAA_Pt[iJetRadius].NewLine(0, 1., 205, 1., 2, 2, kGray);

        TString outputPath = Form("%s/ExclGammaJetCombined/IAA_Pt", options.dataset.Data());
        gSystem->mkdir(outputPath, kTRUE);

        if(mergedPi0) {
            IAA_Pt[iJetRadius].Plot(Form("%s/IAA_Pt_mergedPi0s_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data()), false, false);
        } else {
            IAA_Pt[iJetRadius].Plot(Form("%s/IAA_Pt_isoPhoton_%s.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data()), false, false);
        }
        
        iJetRadius++;
    }

    // IAA vs delta phi
    Plotting1D IAA_Phi[options.jetRadiusFoldersEnum.size()][recoilJetPtBinsIntegration.size()];
    iJetRadius = 0;
     std::map<std::pair<int, int>, TGraphAsymmErrors*> referenceHistosPhi;
    for (auto jetRadius : options.jetRadiusFoldersEnum) {
        GlobalOptions varOptions = options.jetRadiusOptions[iJetRadius];
        TString triggerLabel = getTriggerLabel(varOptions, mergedPi0);

        // reference histogram for this radius and each jet pt
        int irpt = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration) {
            TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
            TGraphAsymmErrors* refHisto = GetHistogramPPRef<TGraphAsymmErrors>(histoName.Data(), options, jetRadius, mergedPi0);
            if (!refHisto) continue;
            refHisto->SetName(Form("%s_%f_%f_ref", histoName.Data(), recoilJetPtBin.first, recoilJetPtBin.second));
            referenceHistosPhi[std::make_pair(iJetRadius, irpt)] = refHisto;
            irpt++;
        }

        int iRecoilJetPtBin = 0;
        for (auto recoilJetPtBin : recoilJetPtBinsIntegration) {
            // loop over options.centralityFoldersEnum
            double maxY = 0;
            for (auto centrality : options.centralityFoldersEnum) {
                // skip reference centrality

                TString histoName = Form("gRecoilJetPtIntegralPhi_%.2f_%.2f", recoilJetPtBin.first, recoilJetPtBin.second);
                TGraphAsymmErrors* histo = GetHistogram<TGraphAsymmErrors>(histoName.Data(), options, centrality, jetRadius, mergedPi0);
                if (!histo) continue;

                // calculate ratio to peripheral

                TGraphAsymmErrors* ratioHisto = CalculateRatioTGraphAsymmErrors(histo, referenceHistosPhi[std::make_pair(iJetRadius, iRecoilJetPtBin)], Form("%s_%f_%f_ratio", histo->GetName(), recoilJetPtBin.first, recoilJetPtBin.second));
                IAA_Phi[iJetRadius][iRecoilJetPtBin].New(ratioHisto, Form("%s/pp", centralityLabels[centrality].Data()), 20, 2, fancyColors[centrality], "pe");

                for (int i = 0; i < ratioHisto->GetN(); i++){
                    if (ratioHisto->GetY()[i] > maxY){
                        maxY = ratioHisto->GetY()[i];
                    }
                }   
            }

            // Set plot labels and ranges
            IAA_Phi[iJetRadius][iRecoilJetPtBin].SetLegend(0.25, 0.45, 0.2, 0.35, true);
            IAA_Phi[iJetRadius][iRecoilJetPtBin].SetMargins(0.12, 0.17, 0.025, 0.025, 2000, 1500, 0.33);
            IAA_Phi[iJetRadius][iRecoilJetPtBin].SetAxisRange(0, TMath::Pi(), 0, 1.8);
            IAA_Phi[iJetRadius][iRecoilJetPtBin].NewLatex(0.25, 0.9, Form("#bf{ALICE work in progress};%s;%s;ch. jets, anti-k_{T} %s;%.1f < p_{T,jet}^{ch} < %.1f GeV/c;|#Delta#varphi-#pi| < 0.6", options.dataset.Data(), triggerLabel.Data(), jetRadiusLabels[jetRadius].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second));
            IAA_Phi[iJetRadius][iRecoilJetPtBin].SetAxisLabel("#Delta#varphi", "I_{AA}");

            // draw dashed line at 1
            IAA_Phi[iJetRadius][iRecoilJetPtBin].NewLine(0, 1., TMath::Pi(), 1., 2, 2, kGray);

            TString outputPath = Form("%s/ExclGammaJetCombined/IAA_Phi", options.dataset.Data());
            gSystem->mkdir(outputPath, kTRUE);  

            if(mergedPi0) {
                IAA_Phi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/IAA_Phi_mergedPi0s_%s_%.1f_%.1f.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second), false, false);
            } else {
                IAA_Phi[iJetRadius][iRecoilJetPtBin].Plot(Form("%s/IAA_Phi_isoPhoton_%s_%.1f_%.1f.pdf", outputPath.Data(), jetRadiusFileStrings[jetRadius].Data(), (float)recoilJetPtBin.first, (float)recoilJetPtBin.second), false, false);
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
    // cout << "File path: " << options.GetInputFilePath(k0_10, kR02) << endl;
    // cout << "Root file: " << options.GetRootFile(k0_10, kR02) << endl;

    drawDeltaPhiDistributions(options, true);
    drawDeltaPhiDistributions(options, false);

    drawCRefPhiDistributions(options, true);
    drawCRefPhiDistributions(options, false);

    drawRecoilJetPtDistributions(options, kR04, true);
    drawRecoilJetPtDistributions(options, kR04, false);
    drawRecoilJetPtDistributions(options, kR05, true);
    drawRecoilJetPtDistributions(options, kR05, false);

    // drawICPDistributions(options, k30_50, true);
    // drawICPDistributions(options, k30_50, false);
    drawICPDistributions(options, k50_90, true);
    drawICPDistributions(options, k50_90, false);

    if(options.ppRefAvailable){
        drawIAAdistributions(options, true);
        drawIAAdistributions(options, false);
    }


    
}