#include "PlottingClass.h"

void plotEfficiencyPurityAcceptance(string filepath1="PuritiesEfficienciesResults/MergedGJ18DandJJ18D/102_Charged/Standard", bool doefficiency=true, bool dopurity=true, bool doacceptance=true, bool doABCD=true, bool doPurityComp=true, bool doAllEfficiencies=true){
    //Setting plot linewidths. (global setting).
    gStyle->SetLineWidth(2);
    TFile* file1=(TFile*)TFile::Open(Form("%s/AnlysedHistosFromTree.root",filepath1.c_str()));

    TString TriggerString="";
    TString ColorString="";
    TString DataString="";

    if(filepath1.find("100_Charged")!=string::npos){
        TriggerString="Trigger: kINT7";
        ColorString="#b30000";
    }else if(filepath1.find("101_Charged")!=string::npos){
        TriggerString="Trigger: EG2";
        ColorString="#1a53ff";
    }else if(filepath1.find("102_Charged")!=string::npos){
        TriggerString="Trigger: EG1";
        ColorString="#007700";
    }else{
        WARN("No Trigger registered in folder name.")
        TriggerString="?";
        ColorString="kBLACK";
    }
    if(filepath1.find("JJ18D")!=string::npos||filepath1.find("GJ18D")!=string::npos||filepath1.find("MergedGJ18DandJJ18D")!=string::npos){
        if(filepath1.find("JJ18D")!=string::npos){
            DataString="JJ";
        }if(filepath1.find("GJ18D")!=string::npos){
            DataString="#gamma-J";
        }if(filepath1.find("MergedGJ18DandJJ18D")!=string::npos){
            DataString="#gamma-J + JJ";
    }
    }else{
        WARN("No dataset registered in folder name.")
        DataString="?";
    }

    //Do efficency plot
    if(doefficiency){
        //Load hist:
        TH1F* effplot= (TH1F*)file1->Get("FullEfficiency");
        //initiate canvas
        TCanvas *ceff = new TCanvas("Efficiency","Efficiency",800,600);
        ceff->Draw();
        ceff->cd();
        //Initiate TPad and set style:
        TPad* padeff = new TPad("pad", "pad",0,1,0,1,0);
        padeff->SetLeftMargin(0.125);
        padeff->SetFillStyle(4000);
        padeff->SetBorderSize(2);
        padeff->SetLineWidth(2);
        padeff->Draw();
        padeff->cd();
        //SetFigure style:
        effplot->SetFillStyle(0);
        effplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        effplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        effplot->SetMarkerStyle(kFullSquare);
        effplot->SetMarkerSize(2);
        effplot->SetLineWidth(2);
        //Set y-axis
        effplot->GetYaxis()->SetTitleFont(43);
        effplot->GetYaxis()->SetLabelFont(43);
        effplot->GetYaxis()->SetTitleSize(24);
        effplot->GetYaxis()->SetLabelSize(18);
        effplot->GetYaxis()->SetRangeUser(0.22,0.60);
        //Set x-axis
        effplot->GetXaxis()->SetTitleFont(43);
        effplot->GetXaxis()->SetLabelFont(43);
        effplot->GetXaxis()->SetTitleSize(24);
        effplot->GetXaxis()->SetLabelSize(18);
        effplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        effplot->SetTitle(";#bf{#it{p}_{T}};#it{#bf{#epsilon}}");
        effplot->SetStats(0);
        effplot->Draw();

        //Text:
        TLatex * efftext = new TLatex (40,0.3,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        efftext->SetTextFont(43);
        efftext->SetTextSize(24);
        efftext->Draw("same");
        TLatex * efftext2 = new TLatex (40,0.4,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        efftext2->SetTextFont(43);
        efftext2->SetTextSize(24);
        efftext2->Draw("same");

        ceff->SaveAs(Form("%s/Efficiency.png",filepath1.c_str()));
        if(doAllEfficiencies){
            //Load hist:
            TH1F* effclustercutsplot= (TH1F*)file1->Get("EfficiencyClusterCuts");
            TH1F* effclusterisolationcutsplot= (TH1F*)file1->Get("EfficiencyClusterIsolationCuts");
            //initiate canvas
            TCanvas *ceffcomp = new TCanvas("EfficiencyComp","EfficiencyComp",800,600);
            ceffcomp->Draw();
            ceffcomp->cd();
            //Initiate TPad and set style:
            TPad* padeffcomp = new TPad("padeffcomp", "padeffcomp",0,1,0,1,0);
            padeffcomp->SetLeftMargin(0.125);
            padeffcomp->SetFillStyle(4000);
            padeffcomp->SetBorderSize(2);
            padeffcomp->SetLineWidth(2);
            padeffcomp->Draw();
            padeffcomp->cd();

            //SetFigure style:
            effclustercutsplot->SetFillStyle(0);
            effclustercutsplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
            effclustercutsplot->SetLineColor(TColor::GetColor(ColorString.Data()));
            effclustercutsplot->SetMarkerStyle(kOpenCircle);
            effclustercutsplot->SetMarkerSize(2);
            effclustercutsplot->SetLineWidth(2);

            effclusterisolationcutsplot->SetFillStyle(0);
            effclusterisolationcutsplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
            effclusterisolationcutsplot->SetLineColor(TColor::GetColor(ColorString.Data()));
            effclusterisolationcutsplot->SetMarkerStyle(kOpenDiamond);
            effclusterisolationcutsplot->SetMarkerSize(2);
            effclusterisolationcutsplot->SetLineWidth(2);

            effplot->GetYaxis()->SetRangeUser(0,1);
            effplot->DrawCopy();
            effclustercutsplot->Draw("same");
            effclusterisolationcutsplot->Draw("same");

            TLatex * efftext3 = new TLatex (50,0.2,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
            efftext3->SetTextFont(43);
            efftext3->SetTextSize(24);
            efftext3->Draw("same");
            TLatex * efftext4 = new TLatex (50,0.1,Form("%s, %s", DataString.Data(),TriggerString.Data()));
            efftext4->SetTextFont(43);
            efftext4->SetTextSize(24);
            efftext4->Draw("same");

            TLegend *effcompleg=new TLegend();
            effcompleg->SetTextFont(43);
            effcompleg->SetBorderSize(0);
            effcompleg->SetTextSize(24);
            effcompleg->AddEntry(effplot,"Full","flp");
            effcompleg->AddEntry(effclustercutsplot,"Cluster Cuts","flp");
            effcompleg->AddEntry(effclusterisolationcutsplot,"Isolation Cuts","flp");

            effcompleg->Draw();

            ceffcomp->SaveAs(Form("%s/EfficiencyComp.png",filepath1.c_str()));

        }

    }
    //Do acceptance plot
    if(doacceptance){
        TH1F* accplot= (TH1F*)file1->Get("Acceptance;1");
        TCanvas *cacc = new TCanvas("Acceptance","Acceptance",800,600);
        cacc->cd();

         //Initiate TPad and set style:
        TPad* padacc = new TPad("padacc", "padacc",0,1,0,1,0);
        padacc->SetLeftMargin(0.125);
        padacc->SetFillStyle(4000);
        padacc->SetBorderSize(2);
        padacc->SetLineWidth(2);
        padacc->Draw();
        padacc->cd();
        //SetFigure style:
        accplot->SetFillStyle(0);
        accplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        accplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        accplot->SetMarkerStyle(kFullSquare);
        accplot->SetMarkerSize(2);
        accplot->SetLineWidth(2);
        //Set y-axis
        accplot->GetYaxis()->SetTitleFont(43);
        accplot->GetYaxis()->SetLabelFont(43);
        accplot->GetYaxis()->SetTitleSize(24);
        accplot->GetYaxis()->SetLabelSize(18);
        //Set x-axis
        accplot->GetXaxis()->SetTitleFont(43);
        accplot->GetXaxis()->SetLabelFont(43);
        accplot->GetXaxis()->SetTitleSize(24);
        accplot->GetXaxis()->SetLabelSize(18);
        accplot->GetXaxis()->SetRangeUser(10.,80.);
        //Drawplot
        accplot->SetTitle(";#bf{#it{p}_{T}};#bf{Acceptance}");
        accplot->GetYaxis()->SetRangeUser(0.,1.0);
        accplot->SetStats(0);
        accplot->Draw();

        //Text:
        TLatex * acctext = new TLatex (40,0.8,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        acctext->SetTextFont(43);
        acctext->SetTextSize(24);
        acctext->Draw("same");
        TLatex * acctext2 = new TLatex (40,0.9,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        acctext2->SetTextFont(43);
        acctext2->SetTextSize(24);
        acctext2->Draw("same");
        cacc->SaveAs(Form("%s/acceptance.png",filepath1.c_str()));
    }
    //Do purity plot
    if(dopurity){
        TH1F* purplot= (TH1F*)file1->Get("Purity");
        TCanvas *cpur = new TCanvas("Purity","Purity",800,600);
        cpur->cd();
        //TPad* pad = new TPad("pad", "pad", 0., 0., 0.97, 0.97);
        //Figure::SetPadStyle(pad);
        //pad->cd();
        purplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        purplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        purplot->SetMarkerStyle(kFullSquare);
        purplot->SetMarkerSize(2);
        //Set y-axis
        purplot->GetYaxis()->SetTitleFont(43);
        purplot->GetYaxis()->SetLabelFont(43);
        purplot->GetYaxis()->SetTitleSize(24);
        purplot->GetYaxis()->SetLabelSize(18);

        //Set x-axis
        purplot->GetXaxis()->SetTitleFont(43);
        purplot->GetXaxis()->SetLabelFont(43);
        purplot->GetXaxis()->SetTitleSize(24);
        purplot->GetXaxis()->SetLabelSize(18);
        purplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        purplot->SetTitle(";#bf{#it{p}_{T}};#bf{Purity}");
        purplot->SetStats(0);
        purplot->Draw();

        //Text:
        TLatex * purtext = new TLatex (40,0.35,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        purtext->SetTextFont(43);
        purtext->SetTextSize(24);
        purtext->Draw("same");
        TLatex * purtext2 = new TLatex (40,0.4,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        purtext2->SetTextFont(43);
        purtext2->SetTextSize(24);
        purtext2->Draw("same");
        //Save plot
        cpur->SaveAs(Form("%s/purity.png",filepath1.c_str()));

    }
    //Do ABCD plot
    if(doABCD){
        //PABCDraw plot
        TDirectory* ABCDdir=(TDirectory*)file1->Get("ABCD");
        TH1F* PABCDrawplot= (TH1F*)ABCDdir->Get("hPABCDraw");
        TCanvas *cPABCDrawplot = new TCanvas("PABCDraw","PABCDraw",800,600);
        cPABCDrawplot->cd();
        //TPad* pad = new TPad("pad", "pad", 0., 0., 0.97, 0.97);
        //Figure::SetPadStyle(pad);
        //pad->cd();
        PABCDrawplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        PABCDrawplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        PABCDrawplot->SetMarkerStyle(kFullSquare);
        PABCDrawplot->SetMarkerSize(2);
        //Set y-axis
        PABCDrawplot->GetYaxis()->SetTitleFont(43);
        PABCDrawplot->GetYaxis()->SetLabelFont(43);
        PABCDrawplot->GetYaxis()->SetTitleSize(24);
        PABCDrawplot->GetYaxis()->SetLabelSize(18);
        PABCDrawplot->GetYaxis()->SetRangeUser(-0.5,1.2);
        //Set x-axis
        PABCDrawplot->GetXaxis()->SetTitleFont(43);
        PABCDrawplot->GetXaxis()->SetLabelFont(43);
        PABCDrawplot->GetXaxis()->SetTitleSize(24);
        PABCDrawplot->GetXaxis()->SetLabelSize(18);
        PABCDrawplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDrawplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{P_{ABCD,raw}}}");
        PABCDrawplot->SetStats(0);
        PABCDrawplot->Draw();

        //Text:
        TLatex * PABCDrawtext = new TLatex (45,0.18,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDrawtext->SetTextFont(43);
        PABCDrawtext->SetTextSize(24);
        PABCDrawtext->Draw("same");
        TLatex * PABCDrawtext2 = new TLatex (45,0.0,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        PABCDrawtext2->SetTextFont(43);
        PABCDrawtext2->SetTextSize(24);
        PABCDrawtext2->Draw("same");
        cPABCDrawplot->SaveAs(Form("%s/ABCDraw.png",filepath1.c_str()));

        //alfa
        //PABCDraw plot
        TH1F* PABCDalfaplot= (TH1F*)ABCDdir->Get("hPABCDalfa");
        TCanvas *cPABCDalfaplot = new TCanvas("PABCDalfa","PABCDalfa",800,600);
        cPABCDalfaplot->cd();
        //TPad* pad = new TPad("pad", "pad", 0., 0., 0.97, 0.97);
        //Figure::SetPadStyle(pad);
        //pad->cd();
        PABCDalfaplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        PABCDalfaplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        PABCDalfaplot->SetMarkerStyle(kFullSquare);
        PABCDalfaplot->SetMarkerSize(2);
        //Set y-axis
        PABCDalfaplot->GetYaxis()->SetTitleFont(43);
        PABCDalfaplot->GetYaxis()->SetLabelFont(43);
        PABCDalfaplot->GetYaxis()->SetTitleSize(24);
        PABCDalfaplot->GetYaxis()->SetLabelSize(18);
        PABCDalfaplot->GetYaxis()->SetRangeUser(1,2.5);
        //Set x-axis
        PABCDalfaplot->GetXaxis()->SetTitleFont(43);
        PABCDalfaplot->GetXaxis()->SetLabelFont(43);
        PABCDalfaplot->GetXaxis()->SetTitleSize(24);
        PABCDalfaplot->GetXaxis()->SetLabelSize(18);
        PABCDalfaplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDalfaplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{#alpha}}");
        PABCDalfaplot->SetStats(0);
        PABCDalfaplot->Draw();

        //Text:
        TLatex * PABCDalfatext = new TLatex (45,2.2,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDalfatext->SetTextFont(43);
        PABCDalfatext->SetTextSize(24);
        PABCDalfatext->Draw("same");
        TLatex * PABCDalfatext2 = new TLatex (45,2,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        PABCDalfatext2->SetTextFont(43);
        PABCDalfatext2->SetTextSize(24);
        PABCDalfatext2->Draw("same");
        cPABCDalfaplot->SaveAs(Form("%s/ABCDalfa.png",filepath1.c_str()));

        //PABCD
        TH1F* PABCDplot= (TH1F*)ABCDdir->Get("hPABCD");
        TCanvas *cPABCDplot = new TCanvas("PABCD","PABCD",800,600);
        cPABCDplot->cd();
        
        PABCDplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        PABCDplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        PABCDplot->SetMarkerStyle(kFullSquare);
        PABCDplot->SetMarkerSize(2);
        //Set y-axis
        PABCDplot->GetYaxis()->SetTitleFont(43);
        PABCDplot->GetYaxis()->SetLabelFont(43);
        PABCDplot->GetYaxis()->SetTitleSize(24);
        PABCDplot->GetYaxis()->SetLabelSize(18);
        PABCDplot->GetYaxis()->SetRangeUser(-0.5,1.2);
        //Set x-axis
        PABCDplot->GetXaxis()->SetTitleFont(43);
        PABCDplot->GetXaxis()->SetLabelFont(43);
        PABCDplot->GetXaxis()->SetTitleSize(24);
        PABCDplot->GetXaxis()->SetLabelSize(18);
        PABCDplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{P_{ABCD}}}");
        PABCDplot->SetStats(0);
        PABCDplot->DrawCopy();

        //Text:
        TLatex * PABCDtext = new TLatex (45,0.18,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDtext->SetTextFont(43);
        PABCDtext->SetTextSize(24);
        PABCDtext->Draw("same");
        TLatex * PABCDtext2 = new TLatex (45,0,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        PABCDtext2->SetTextFont(43);
        PABCDtext2->SetTextSize(24);
        PABCDtext2->Draw("same");
        cPABCDplot->SaveAs(Form("%s/ABCD.png",filepath1.c_str()));
    }

    if(doABCD && dopurity && doPurityComp){
        TH1F* purplot= (TH1F*)file1->Get("Purity");

        purplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        purplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        purplot->SetMarkerStyle(kFullSquare);
        purplot->SetMarkerSize(2);
        //Set y-axis
        purplot->GetYaxis()->SetTitleFont(43);
        purplot->GetYaxis()->SetLabelFont(43);
        purplot->GetYaxis()->SetTitleSize(24);
        purplot->GetYaxis()->SetLabelSize(18);

        //Set x-axis
        purplot->GetXaxis()->SetTitleFont(43);
        purplot->GetXaxis()->SetLabelFont(43);
        purplot->GetXaxis()->SetTitleSize(24);
        purplot->GetXaxis()->SetLabelSize(18);
        purplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        purplot->SetTitle(";#bf{#it{p}_{T}};#bf{Purity}");
        purplot->SetStats(0);

        TDirectory* ABCDdir=(TDirectory*)file1->Get("ABCD");
        TH1F* PABCDplot= (TH1F*)ABCDdir->Get("hPABCD");
        
        PABCDplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        PABCDplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        PABCDplot->SetMarkerStyle(kOpenSquare);
        PABCDplot->SetMarkerSize(2);
        //Set y-axis
        PABCDplot->GetYaxis()->SetTitleFont(43);
        PABCDplot->GetYaxis()->SetLabelFont(43);
        PABCDplot->GetYaxis()->SetTitleSize(24);
        PABCDplot->GetYaxis()->SetLabelSize(18);
        PABCDplot->GetYaxis()->SetRangeUser(-0.5,1.2);
        //Set x-axis
        PABCDplot->GetXaxis()->SetTitleFont(43);
        PABCDplot->GetXaxis()->SetLabelFont(43);
        PABCDplot->GetXaxis()->SetTitleSize(24);
        PABCDplot->GetXaxis()->SetLabelSize(18);
        PABCDplot->GetXaxis()->SetRangeUser(10.,80.);
        PABCDplot->SetStats(0);

        TLegend *purcompleg=new TLegend();
        purcompleg->SetTextFont(43);
        purcompleg->SetBorderSize(0);
        purcompleg->SetTextSize(24);
        purcompleg->AddEntry(PABCDplot,"#it{P_{abcd}}","flp");
        purcompleg->AddEntry(purplot,"#it{Purity(MC)}","flp");



        TCanvas *cpuritycomp = new TCanvas("PurityComparison","PurityComparison",800,600);
        cpuritycomp->cd();
        PABCDplot->DrawCopy();
        purplot->DrawCopy("same");
        purcompleg->Draw();

        TLatex * purcomplegtext = new TLatex (45,0.2,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        purcomplegtext->SetTextFont(43);
        purcomplegtext->SetTextSize(24);
        purcomplegtext->Draw("same");
        TLatex * purcomplegtext2 = new TLatex (45,0.05,Form("%s, %s", DataString.Data(),TriggerString.Data()));
        purcomplegtext2->SetTextFont(43);
        purcomplegtext2->SetTextSize(24);
        purcomplegtext2->Draw("same");

        cpuritycomp->SaveAs(Form("%s/PurityComparison.png",filepath1.c_str()));

    }


}