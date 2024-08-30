#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h" 
#include "TFractionFitter.h"

void PlotsAndTemplateFitPtBinned(const char* Data, const char* Template, const char* Output, TString LegendName = "Data", Double_t BDTcut=-std::numeric_limits<float>::max(), Double_t RandFloat=0){
    //gStyle->SetStats(0);
//Load data
    TFile* fData = TFile::Open(Data,"READ");
    TTree* tData = (TTree*)fData->Get("caloclustertree;1");
//Load template data
    TFile* fTemplate = TFile::Open(Template,"READ");
    TTree* tTemplate = (TTree*)fTemplate->Get("caloclustertree;1");

//Create output file
    TString OutName = Form("%s_TF.root",Output);
    if(BDTcut!=-std::numeric_limits<float>::max()){
        OutName = Form("%s_%0.2f_TF.root",Output,BDTcut);
    }
    TFile* fOutput = TFile::Open(OutName,"RECREATE");
//Vectors for holding fit, integrals and pTbins in output file
    std::vector<Double_t> vecpar0;
    std::vector<Double_t> vecpar0err;
    std::vector<Double_t> vecpar1;
    std::vector<Double_t> vecpar1err;
    std::vector<Double_t> vecSignalIntegral;
    std::vector<Double_t> vecSignalIntegralerr;
    std::vector<Double_t> vecFullSignalIntegral;
    std::vector<Double_t> vecFullSignalIntegralerr;
    std::vector<Double_t> vecTotIntegral;
    std::vector<Double_t> vecTotIntegralerr;
    std::vector<Double_t> vecpTbins;
    std::vector<Double_t> vecpurity;
    std::vector<Double_t> vecpurityerr;
    std::vector<Double_t> vecefficiency;
    std::vector<Double_t> vecefficiencyerr;
    

//Begin loop over pTbins:
    int npTbinEdges=23;
    float pTbinEdges[23] = {3,4,4.5,5,5.5,6,6.5,7,7.5,8,9,10,11,12.5,15,16.25,17.5,18.25,20,25,30,35,40};
    //int npTbinEdges=11;
    //float pTbinEdges[14] = {8,9,10,11,12.5,15,16.25,17.5,18.25,20,25,30,35,40};
    //int npTbinEdges=10;
    //float pTbinEdges[10] = {3,4,4.5,5,5.5,6,6.5,7,7.5,8};
    
    Double_t purity[npTbinEdges-1];
    Double_t purityerr[npTbinEdges-1];
    Double_t purityratio[npTbinEdges-1];
    Double_t x[npTbinEdges-1];
    Double_t purityMC[npTbinEdges-1];

    //Not actual efficiency, has to be scaled by efficiency of IsoCut.
    Double_t efficiency[npTbinEdges-1];
    Double_t efficiencyerr[npTbinEdges-1];
    



    for(int ipt=0; ipt<npTbinEdges-1; ipt++){
        x[ipt]=pTbinEdges[ipt]+(pTbinEdges[ipt+1]-pTbinEdges[ipt])/2;
        cout<<"Current pT-bin:"<<pTbinEdges[ipt]<<"<pT<"<<pTbinEdges[ipt+1]<<"\n";
    //Initiate histograms:
        TH1D* hBDTData = new TH1D(Form("hBDTData%.1f",pTbinEdges[ipt]),Form("hBDTData%.1fGeV<%.1fGeV",pTbinEdges[ipt],pTbinEdges[ipt+1]),100,-10,10);
        float vBDTData;
        float vPtData;
        float vRandfloat;
        float vWeightData;
        tData->SetBranchAddress("BDT", &vBDTData);
        tData->SetBranchAddress("Cluster_SplitFloat", &vRandfloat);
        tData->SetBranchAddress("Cluster_Pt",&vPtData);
        //if(ApplyDataWeights){
        //    cout<<"Applying weights!"<<"\n";
        //    tData->SetBranchAddress("Event_Weight",&vWeightData);
        //}
        tData->SetBranchAddress("Event_Weight",&vWeightData);
        
        
        cout<<"Data"<<"\n";
        int NDataAccepted=0;
        Long64_t nDataentries = tData->GetEntries();
        for (Long64_t i = 0; i < nDataentries; ++i) {
            tData->GetEntry(i);
            if(vRandfloat>RandFloat){
                //cout<<pTbinEdges[i]<<"<"<<vPtData<<"\n";
                //cout<<pTbinEdges[i+1]<<">"<<vPtData<<"\n";
                if(pTbinEdges[ipt]<vPtData && pTbinEdges[ipt+1]>vPtData){
                    NDataAccepted+=1;
                    //if(ApplyDataWeights){
                    //    hBDTData->Fill(vBDTData,vWeightData);
                    //}
                    hBDTData->Fill(vBDTData,vWeightData);
                }
            }
        }
        cout<<"Ndata accepted="<<NDataAccepted<<"/"<<nDataentries<<"\n";
        TH1D* hBDTTemplateSignal = new TH1D(Form("hBDTTemplateSignalpT%.3f",pTbinEdges[ipt]),Form("hBDTTemplateSignalpT%.3f",pTbinEdges[ipt]),100,-10,10);
        TH1D* hBDTTemplateBKG = new TH1D(Form("hBDTTemplateBKGpT%.3f",pTbinEdges[ipt]),Form("hBDTTemplateBKGpT%.3f",pTbinEdges[ipt]),100,-10,10);

        float vBDTTemplate;
        float vPtTemplate;
        float vWeight;
        bool vIsSignalTemplate;
        tTemplate->SetBranchAddress("BDT", &vBDTTemplate);
        tTemplate->SetBranchAddress("Cluster_isSignal", &vIsSignalTemplate);
        tTemplate->SetBranchAddress("Cluster_Pt",&vPtTemplate);
        tTemplate->SetBranchAddress("Event_Weight",&vWeight);
        cout<<"Template"<<"\n";
        Long64_t nTemplateentries = tTemplate->GetEntries();
        for (Long64_t i = 0; i < nTemplateentries; ++i) {
            tTemplate->GetEntry(i);
            if(pTbinEdges[ipt]<vPtTemplate && pTbinEdges[ipt+1]>vPtTemplate){
                if(vIsSignalTemplate==true){
                    hBDTTemplateSignal->Fill(vBDTTemplate,vWeight);
                }
                if(vIsSignalTemplate==false){
                    hBDTTemplateBKG->Fill(vBDTTemplate,vWeight);
                }
            }
            
        }

        Int_t Nbins=hBDTData->GetXaxis()->GetNbins();

        //Normalise histograms:
        Double_t ndat = hBDTData->Integral(0, Nbins);
        Double_t nsig = hBDTTemplateSignal->Integral(0, Nbins);
        Double_t nbac = hBDTTemplateBKG->Integral(0, Nbins);

        cout<<"sigfrac="<<nsig/(nsig+nbac)<<"\n";
        cout<<"bkgfrac="<<nbac/(nsig+nbac)<<"\n";

        purityMC[ipt]=nsig/(nsig+nbac);

        //hBDTTemplateSignal->Scale(ndat / nsig);
        //hBDTTemplateBKG->Scale(ndat / nbac);

        hBDTTemplateSignal->Write();
        hBDTTemplateBKG->Write();
        hBDTData->Write();

        //template fit:
        TObjArray *mcFit = new TObjArray(2);
        mcFit->Add(hBDTTemplateSignal);
        mcFit->Add(hBDTTemplateBKG);
        TFractionFitter* fit = new TFractionFitter(hBDTData, mcFit);

        fit->Constrain(0,0.0,1.0);
        //fit->Constrain(1,0.0,1.0);
        fit->SetRangeX(1,Nbins);

        for (Int_t b = 0; b <Nbins; b++){
           Double_t error =hBDTData->GetBinError(b);
           Double_t errorsig =1;//hBDTTemplateSignal->GetBinError(b);
           Double_t errorbkg =1;//hBDTTemplateBKG->GetBinError(b);
           if(error==0){//||errorsig==0||errorbkg==0
                fit->ExcludeBin(b);
                cout<<"Excluded bin:"<<b<<"/"<<Nbins<<"\n";
                cout<<"Error:"<<error<<"\n";
                cout<<"Signal Error:"<<errorsig<<"\n";
                cout<<"Background Error:"<<errorbkg<<"\n";
           }
        }

        Int_t status = fit->Fit();
        std::cout << "fit status: " << status << std::endl;
        if (status == 0) {                       // check on fit status
        //Plot fit:
            TCanvas *c=new TCanvas("c","c",800,600);
            c->Draw();
            c->cd();
            TPad *p=new TPad("p","p",0,0.4,1,1);
            p->SetBottomMargin(0);
            p->Draw();
            p->cd();

            hBDTData->SetStats(0);
            hBDTData->SetDirectory(0);
            hBDTData->SetMarkerStyle(kFullSquare);
            hBDTData->SetMarkerColor(kBlue);
            hBDTData->SetMarkerSize(1.6);
            hBDTData->GetXaxis()->SetTitle("BDT");
            hBDTData->GetYaxis()->SetTitle("Counts");
            hBDTData->SetTitle(Form("hBDTDFit%.1f;BDT;Counts",pTbinEdges[ipt]));
            hBDTData->Draw("P");

            TH1F* result = (TH1F*) fit->GetPlot();

            Double_t par0;
            Double_t par0err;
            fit->GetResult(0,par0,par0err);

            hBDTTemplateSignal->Scale((ndat / nsig)*par0);
            hBDTTemplateSignal->SetMarkerStyle(kFullSquare);
            hBDTTemplateSignal->SetMarkerColor(kGreen+3);
            hBDTTemplateSignal->Draw("P same");

            Double_t par1;
            Double_t par1err;
            fit->GetResult(1,par1,par1err);

            hBDTTemplateBKG->Scale((ndat / nbac)*par1);
            hBDTTemplateBKG->SetMarkerStyle(kFullSquare);
            hBDTTemplateBKG->SetMarkerColor(kViolet+3);
            hBDTTemplateBKG->Draw("P same");

            TH1D* hBDTfit = (TH1D*)hBDTTemplateSignal->Clone();
            hBDTfit->SetName(Form("hBDTDFit%.1f",pTbinEdges[ipt]));
            hBDTfit->Add(hBDTTemplateBKG);
            hBDTfit->SetMarkerSize(0);
            hBDTfit->SetLineColor(kRed);
            hBDTfit->SetLineWidth(2);
            hBDTfit->Draw("1c same");

            TLegend* legend = new TLegend(0.6,0.7,0.9,0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(Form("hBDTDFit%.1f",pTbinEdges[ipt]),"fit","lep");
            legend->AddEntry(Form("hBDTData%.1f",pTbinEdges[ipt]),"Data","lep");
            legend->AddEntry(Form("hBDTTemplateSignalpT%.3f",pTbinEdges[ipt]),Form("Signal:%.3f#pm%.3f",par0,par0err),"lep");
            legend->AddEntry(Form("hBDTTemplateBKGpT%.3f",pTbinEdges[ipt]),Form("Background:%.3f#pm%.3f",par1,par1err),"lep");
            legend->Draw("NB");

            //Plot ratio
            c->cd();
            TPad *p2=new TPad("p2","p2",0,0.,1,0.4);
            p2->Draw();
            p2->cd();
            p2->SetTopMargin(0.);
            p2->SetBottomMargin(0.3);
            TH1D* hBDTratio=(TH1D*)hBDTfit->Clone();
            hBDTratio->Divide(hBDTData);
            hBDTratio->SetTitle(";BDT;Fit/Data");
            //hBDTratio->SetLabelFont(43);
            //hBDTratio->SetLabelFont(26);
            hBDTratio->Draw();

            c->Print(Form("%s_TF_ResultPlotPt%f.pdf",Output,pTbinEdges[ipt]));

            result->Write();
            delete c;

            //Integrating and calculating uncertainty:
            Double_t FullSignalIntegralErr;
            Double_t SignalIntegralErr;
            Double_t BkgIntegralErr;
            Double_t TotIntegralErr;

            Double_t FullSignalIntegral = hBDTTemplateSignal->IntegralAndError(1, Nbins,FullSignalIntegralErr);
            Double_t SignalIntegral = hBDTTemplateSignal->IntegralAndError(hBDTData->FindBin(BDTcut), Nbins,SignalIntegralErr);
            Double_t BkgIntegral = hBDTTemplateBKG->IntegralAndError(hBDTData->FindBin(BDTcut), Nbins,BkgIntegralErr);
            Double_t TotIntegral = result->IntegralAndError(hBDTData->FindBin(BDTcut), Nbins,TotIntegralErr);

            //Calculating purity:
            purity[ipt]=SignalIntegral/(SignalIntegral+BkgIntegral);//(par1*hBDTTemplateBKG->Integral(hBDTData->FindBin(BDTcut), Nbins)+par0*hBDTTemplateSignal->Integral(hBDTData->FindBin(BDTcut), Nbins));
            purityratio[ipt]=purity[ipt]/purityMC[ipt];

            //Errorpropagation:
            Double_t A=SignalIntegralErr/TotIntegral;
            Double_t B=SignalIntegral*TotIntegralErr/(TotIntegral*TotIntegral);
            purityerr[ipt]=TMath::Sqrt(A*A+B*B);

            //Calculating Efficiency: (Not true efficiency, must be scaled by ptIso cut efficiency):
            efficiency[ipt]=SignalIntegral/FullSignalIntegral;
            //Errorpropagation:
            Double_t A2=SignalIntegralErr/FullSignalIntegral;
            Double_t B2=SignalIntegral*FullSignalIntegralErr/(FullSignalIntegral*FullSignalIntegral);
            efficiencyerr[ipt]=TMath::Sqrt(A2*A2+B2*B2);
 
            //cout<<"Find Bin:"<<hBDTData->FindBin(BDTcut)<<"\n";
            //cout<<"Purity:"<<purity[ipt]<<"\n";
            //cout<<"Signal Integral:"<<SignalIntegral<<"#pm"<<SignalIntegralErr<<"\n";
            //cout<<"Result Integral:"<<TotIntegral<<"#pm"<<TotIntegralErr<<"\n";

            //Saving results in vectors:
            vecpar0.push_back(par0);
            vecpar0err.push_back(par0err);
            vecpar1.push_back(par1);
            vecpar1err.push_back(par1err);
            vecFullSignalIntegral.push_back(FullSignalIntegral);
            vecFullSignalIntegralerr.push_back(FullSignalIntegralErr);
            vecSignalIntegral.push_back(SignalIntegral);
            vecSignalIntegralerr.push_back(SignalIntegralErr);
            vecTotIntegral.push_back(TotIntegral);
            vecTotIntegralerr.push_back(TotIntegralErr);
            vecpTbins.push_back(pTbinEdges[ipt]);
            //If last iteration, save last pTbinedge:
            if(ipt==npTbinEdges-2){
                vecpTbins.push_back(pTbinEdges[ipt+1]);
            }
            vecpurity.push_back(purity[ipt]);
            vecpurityerr.push_back(purityerr[ipt]);
            vecefficiency.push_back(efficiency[ipt]);
            vecefficiencyerr.push_back(efficiencyerr[ipt]);
        }
        //Release memory:
        delete hBDTTemplateSignal;
        delete hBDTTemplateBKG;
        delete hBDTData;
        delete fit;

    }

    //Write results to outputfile:
    fOutput->WriteObject(&vecpar0,"fitresult_par0");
    fOutput->WriteObject(&vecpar0err,"fitresult_par0err");
    fOutput->WriteObject(&vecpar0,"fitresult_par1");
    fOutput->WriteObject(&vecpar0err,"fitresult_par1err");
    fOutput->WriteObject(&vecSignalIntegral,"Integral_Signal");
    fOutput->WriteObject(&vecSignalIntegralerr,"Integralerr_Signal");
    fOutput->WriteObject(&vecTotIntegral,"Integral_Tot");
    fOutput->WriteObject(&vecTotIntegralerr,"Integralerr_Tot");
    fOutput->WriteObject(&vecpTbins,"pTbinedges");
    fOutput->WriteObject(&vecpurity,"purity");
    fOutput->WriteObject(&vecpurity,"purityerr");

    


//Plot purity:
    TCanvas *c1=new TCanvas("c","c",800,600);
    c1->cd();
    TGraphErrors* gpurity=new TGraphErrors(npTbinEdges-1,x,purity,0,purityerr);
    gpurity->SetTitle("MC purity;#it{p_{T}};#it{P_{ML}}");
    gpurity->SetMarkerStyle(kFullSquare);
    gpurity->SetMarkerColor(kViolet);
    gpurity->SetLineColor(kViolet);
    gpurity->Draw("APE");

    TGraph* gpurityMC=new TGraph(npTbinEdges-1,x,purityMC);
    gpurityMC->SetMarkerStyle(kFullDiamond);
    gpurityMC->SetMarkerColor(kRed);
    gpurityMC->Draw("p same");

    TLegend* gpuritylegend= new TLegend();
    gpuritylegend->AddEntry(gpurity,LegendName.Data());
    gpuritylegend->AddEntry(gpurityMC,LegendName.Data());


    if(BDTcut==-std::numeric_limits<float>::max()){
        c1->Print(Form("%s_No_Cut_TF_ResultPlotPurity.pdf",Output));
    }else{
        c1->Print(Form("%s_%f_TF_ResultPlotPurity.pdf",Output,BDTcut));
    }
//Plot purity ratio:
    TCanvas *c2=new TCanvas("c2","c2",800,600);
    c2->cd();
    TGraphErrors* gpurityratio=new TGraphErrors(npTbinEdges-1,x,purityratio,0,purityerr);
    gpurityratio->SetTitle("purityratio;#it{p_{T}};#it{P_{ML}}");
    gpurityratio->SetMarkerStyle(kFullSquare);
    gpurityratio->SetMarkerColor(kViolet);
    gpurityratio->SetLineColor(kViolet);
    gpurityratio->Draw("APE");

    c2->Print(Form("%s_PurityRatio.pdf",Output));
    
    fOutput->Close();
    fData->Close();
    fTemplate->Close();

    //Plot Efficiency:
    TCanvas *c3=new TCanvas("c","c",800,600);
    c3->cd();

    TGraphErrors* gefficiency=new TGraphErrors(npTbinEdges-1,x,efficiency,0,efficiencyerr);
    gefficiency->SetTitle("ML efficiency;#it{p_{T}};#it{#epsilon_{ML}}");
    gefficiency->SetMarkerStyle(kFullSquare);
    gefficiency->SetMarkerColor(kViolet);
    gefficiency->SetLineColor(kViolet);
    gefficiency->Draw("APE");

    TLegend* gefficiencylegend= new TLegend();
    gefficiencylegend->AddEntry(gefficiency,LegendName.Data());

    if(BDTcut==-std::numeric_limits<float>::max()){
        c3->Print(Form("%s_No_Cut_TF_ResultPlotEfficiency.pdf",Output));
    }else{
        c3->Print(Form("%s_%f_TF_ResultPlotEfficiency.pdf",Output,BDTcut));
    }

}