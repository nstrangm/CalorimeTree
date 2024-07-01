#include "PlottingClass.h"

void plotEfficiencyPurityAcceptance(string filepath1="PuritiesEfficienciesResults/102_Charged/Standard", bool doefficiency=true, bool dopurity=true, bool doacceptance=true, bool doABCD=true){
    TFile* file1=(TFile*)TFile::Open(Form("%s/AnlysedHistosFromTree.root",filepath1.c_str()));

    TString TriggerString="";
    TString ColorString=""; 

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

    //Do efficency plot
    if(doefficiency){
        TH1F* effplot= (TH1F*)file1->Get("Efficiency");
        TCanvas *ceff = new TCanvas("Efficiency","Efficiency",800,600);
        ceff->cd();
        //TPad* pad = new TPad("pad", "pad", 0., 0., 0.97, 0.97);
        //Figure::SetPadStyle(pad);
        //pad->cd();
        effplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        effplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        effplot->SetMarkerStyle(kFullSquare);
        effplot->SetMarkerSize(2);
        //Set y-axis
        effplot->GetYaxis()->SetTitleFont(43);
        effplot->GetYaxis()->SetLabelFont(43);
        effplot->GetYaxis()->SetTitleSize(24);
        effplot->GetYaxis()->SetLabelSize(18);
        //Set x-axis
        effplot->GetXaxis()->SetTitleFont(43);
        effplot->GetXaxis()->SetLabelFont(43);
        effplot->GetXaxis()->SetTitleSize(24);
        effplot->GetXaxis()->SetLabelSize(18);
        effplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        effplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{#epsilon}}");
        effplot->Draw();

        //Text:
        TLatex * efftext = new TLatex (40,0.0057,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        efftext->SetTextFont(43);
        efftext->SetTextSize(24);
        efftext->Draw("same");
        TLatex * efftext2 = new TLatex (40,0.0051,Form("%s",TriggerString.Data()));
        efftext2->SetTextFont(43);
        efftext2->SetTextSize(24);
        efftext2->Draw("same");

    }
    //Do acceptance plot
    if(doacceptance){
        TH1F* accplot= (TH1F*)file1->Get("Acceptance");
        TCanvas *cacc = new TCanvas("Acceptance","Acceptance",800,600);
        cacc->cd();
        //TPad* pad = new TPad("pad", "pad", 0., 0., 0.97, 0.97);
        //Figure::SetPadStyle(pad);
        //pad->cd();
        accplot->SetMarkerColor(TColor::GetColor(ColorString.Data()));
        accplot->SetLineColor(TColor::GetColor(ColorString.Data()));
        accplot->SetMarkerStyle(kFullSquare);
        accplot->SetMarkerSize(2);
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
        accplot->GetYaxis()->SetRangeUser(0.96,1.0);
        accplot->Draw();

        //Text:
        TLatex * acctext = new TLatex (40,0.975,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        acctext->SetTextFont(43);
        acctext->SetTextSize(24);
        acctext->Draw("same");
        TLatex * acctext2 = new TLatex (40,0.972,Form("%s",TriggerString.Data()));
        acctext2->SetTextFont(43);
        acctext2->SetTextSize(24);
        acctext2->Draw("same");
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
        purplot->Draw();

        //Text:
        TLatex * purtext = new TLatex (20,0.38,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        purtext->SetTextFont(43);
        purtext->SetTextSize(24);
        purtext->Draw("same");
        TLatex * purtext2 = new TLatex (20,0.35,Form("%s",TriggerString.Data()));
        purtext2->SetTextFont(43);
        purtext2->SetTextSize(24);
        purtext2->Draw("same");

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
        //Set x-axis
        PABCDrawplot->GetXaxis()->SetTitleFont(43);
        PABCDrawplot->GetXaxis()->SetLabelFont(43);
        PABCDrawplot->GetXaxis()->SetTitleSize(24);
        PABCDrawplot->GetXaxis()->SetLabelSize(18);
        PABCDrawplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDrawplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{P_{ABCD,raw}}}");
        PABCDrawplot->Draw();

        //Text:
        TLatex * PABCDrawtext = new TLatex (50,0.45,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDrawtext->SetTextFont(43);
        PABCDrawtext->SetTextSize(24);
        PABCDrawtext->Draw("same");
        TLatex * PABCDrawtext2 = new TLatex (50,0.41,Form("%s",TriggerString.Data()));
        PABCDrawtext2->SetTextFont(43);
        PABCDrawtext2->SetTextSize(24);
        PABCDrawtext2->Draw("same");

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
        //Set x-axis
        PABCDalfaplot->GetXaxis()->SetTitleFont(43);
        PABCDalfaplot->GetXaxis()->SetLabelFont(43);
        PABCDalfaplot->GetXaxis()->SetTitleSize(24);
        PABCDalfaplot->GetXaxis()->SetLabelSize(18);
        PABCDalfaplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDalfaplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{P_{ABCD,raw}}}");
        PABCDalfaplot->Draw();

        //Text:
        TLatex * PABCDalfatext = new TLatex (30,2.1,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDalfatext->SetTextFont(43);
        PABCDalfatext->SetTextSize(24);
        PABCDalfatext->Draw("same");
        TLatex * PABCDalfatext2 = new TLatex (30,2.02,Form("%s",TriggerString.Data()));
        PABCDalfatext2->SetTextFont(43);
        PABCDalfatext2->SetTextSize(24);
        PABCDalfatext2->Draw("same");

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
        //Set x-axis
        PABCDplot->GetXaxis()->SetTitleFont(43);
        PABCDplot->GetXaxis()->SetLabelFont(43);
        PABCDplot->GetXaxis()->SetTitleSize(24);
        PABCDplot->GetXaxis()->SetLabelSize(18);
        PABCDplot->GetXaxis()->SetRangeUser(10.,80.);

        //Drawplot
        PABCDplot->SetTitle(";#bf{#it{p}_{T}};#bf{#it{P_{ABCD,raw}}}");
        PABCDplot->Draw();

        //Text:
        TLatex * PABCDtext = new TLatex (50,0.21,"MC, #it{pp} #sqrt{#it{s}_{NN}} = 13 TeV");
        PABCDtext->SetTextFont(43);
        PABCDtext->SetTextSize(24);
        PABCDtext->Draw("same");
        TLatex * PABCDtext2 = new TLatex (50,0.18,Form("%s",TriggerString.Data()));
        PABCDtext2->SetTextFont(43);
        PABCDtext2->SetTextSize(24);
        PABCDtext2->Draw("same");
    }


}