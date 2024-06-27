#ifndef _ANALYSEHISTOSFROMTREE_
#define _ANALYSEHISTOSFROMTREE_

#include "Utilities.h"
#include "HistogramLibrary.h"

class ABCD
{
public:
    //ABCD values with uncertainties:
    //Data:
    Double_t AData = 0;//N isolated & narrow
    Double_t BData = 0;//N isolated & wide
    Double_t CData = 0;//N antiisolated & narrow
    Double_t DData = 0;//N antiisolated & wide
    Double_t ADataErr = 0;
    Double_t BDataErr = 0;
    Double_t CDataErr = 0;
    Double_t DDataErr = 0;
    //MC:
    Double_t AMC = 0;//N isolated & narrow (From background)
    Double_t BMC = 0;//N isolated & wide
    Double_t CMC = 0;//N antiisolated & narrow
    Double_t DMC = 0;//N antiisolated & wide
    Double_t AMCErr = 0;
    Double_t BMCErr = 0;
    Double_t CMCErr = 0;
    Double_t DMCErr = 0;
    //PABCDraw:
    Double_t PABCDraw = 0;
    Double_t PABCDrawErr = 0;
    //PABCD:
    Double_t PABCD = 0;
    Double_t PABCDErr = 0;
};

void CalculateABCDData(std::vector<ABCD> &ABCDs ,TDirectory* dIsoGammaData, TDirectory* dIsoGammaMC, int nbins , float pTbinEdges[], float pTIsoABCDCuts[], float M02ABCDCuts[]){  
    //Load THnSparse:
    //Load MC data:
    THnSparse *PtIsovsM02vsPtMC = (THnSparse*)dIsoGammaMC->Get("hIsoGammaIsovsM02vsPtnoCuts");
    THnSparse *PtIsovsM02vsPtMCBackground = (THnSparse*)dIsoGammaMC->Get("hIsoGammaIsovsM02vsPtnoCutsBackground");
    //Load actual data:
    THnSparse *PtIsovsM02vsPtData = (THnSparse*)dIsoGammaData->Get("hIsoGammaIsovsM02vsPtnoCuts");

    //Loop over pT bins and calculate A,B,C and D:
    for(int bin = 1; bin<=nbins; ++bin){
        ABCD ABCDtemp;
        //Calculating the correction factor alfa
        //For projecting THnsparse
        TH2F *h2M02vsPtprojectionMC;
        TH2F *h2M02vsPtprojectionMCBackground;

        //Restrict pt-range
        PtIsovsM02vsPtMC->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        PtIsovsM02vsPtMCBackground->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        h2M02vsPtprojectionMC = (TH2F*)PtIsovsM02vsPtMC->Projection(0,1,"E");
        h2M02vsPtprojectionMCBackground = (TH2F*)PtIsovsM02vsPtMCBackground->Projection(0,1,"E");
        //find bins
        int pTIsoMinBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(pTIsoABCDCuts[0]);
        int pTIsoMaxBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(pTIsoABCDCuts[1]);
        int pTIsoMinBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(pTIsoABCDCuts[2]);
        int pTIsoMaxBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(pTIsoABCDCuts[3]);

        int M02MinBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(M02ABCDCuts[0]);
        int M02MinBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(M02ABCDCuts[1]);
        int M02MaxBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(M02ABCDCuts[2]);
        int M02MaxBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindFixBin(M02ABCDCuts[3]);

        ABCDtemp.AMC=h2M02vsPtprojectionMCBackground->IntegralAndError(pTIsoMinBin1,pTIsoMaxBin1,M02MinBin1,M02MaxBin1,ABCDtemp.AMCErr);
        ABCDtemp.BMC=h2M02vsPtprojectionMC->IntegralAndError(pTIsoMinBin1,pTIsoMaxBin1,M02MinBin2,M02MaxBin2,ABCDtemp.BMCErr);
        ABCDtemp.CMC=h2M02vsPtprojectionMC->IntegralAndError(pTIsoMinBin2,pTIsoMaxBin2,M02MinBin1,M02MaxBin1,ABCDtemp.CMCErr);
        ABCDtemp.DMC=h2M02vsPtprojectionMC->IntegralAndError(pTIsoMinBin2,pTIsoMaxBin2,M02MinBin2,M02MaxBin2,ABCDtemp.DMCErr);

        //Freeing memory:
        delete h2M02vsPtprojectionMC;
        delete h2M02vsPtprojectionMCBackground;

        //For claculating the PABCD-raw:
        TH2F *h2M02vsPtprojectionData;
        //Restrict Pt-range
        PtIsovsM02vsPtData->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        //Project onto pt-range:
        h2M02vsPtprojectionData = (TH2F*)PtIsovsM02vsPtData->Projection(0,1);
        //Integrate
        ABCDtemp.AData=h2M02vsPtprojectionData->IntegralAndError(pTIsoMinBin1,pTIsoMaxBin1,M02MinBin1,M02MaxBin1,ABCDtemp.ADataErr,"");
        ABCDtemp.BData=h2M02vsPtprojectionData->IntegralAndError(pTIsoMinBin1,pTIsoMaxBin1,M02MinBin2,M02MaxBin2,ABCDtemp.BDataErr,"");
        ABCDtemp.CData=h2M02vsPtprojectionData->IntegralAndError(pTIsoMinBin2,pTIsoMaxBin2,M02MinBin1,M02MaxBin1,ABCDtemp.CDataErr,"");
        ABCDtemp.DData=h2M02vsPtprojectionData->IntegralAndError(pTIsoMinBin2,pTIsoMaxBin2,M02MinBin2,M02MaxBin2,ABCDtemp.DDataErr,"");
        //Freeing memory:
        delete h2M02vsPtprojectionData;

        cout<<ABCDtemp.AData<<"\n";
        cout<<ABCDtemp.BData<<"\n";
        cout<<ABCDtemp.CData<<"\n";
        cout<<ABCDtemp.DData<<"\n";

        //Fill purity:
        ABCDtemp.PABCDraw = 1-(ABCDtemp.AData/ABCDtemp.CData)/(ABCDtemp.BData/ABCDtemp.DData);
        cout<<ABCDtemp.PABCDraw<<"\n";
        ABCDtemp.PABCDrawErr = 0;
        //PABCD:
        ABCDtemp.PABCD = 1-(ABCDtemp.AData/ABCDtemp.CData)/(ABCDtemp.BData/ABCDtemp.DData)*(ABCDtemp.CMC/ABCDtemp.AMC)/(ABCDtemp.DMC/ABCDtemp.DMC);
        cout<<ABCDtemp.PABCD<<"\n";
        ABCDtemp.PABCDErr = 0;

        ABCDs.push_back(ABCDtemp);
    }
}

TDirectory *DefineABCDHistos(TFile *f,int nbins, float pTbinEdges[]){
    TDirectory *dir = f->mkdir("ABCD");
    dir->cd();
    TH1D *hPABCD=new TH1D("hPABCD","P_{ABCD}",nbins, pTbinEdges);
    dir->Add(hPABCD);
    TH1D *hPABCDraw=new TH1D("hPABCDraw","P_{ABCD,RAW}",nbins, pTbinEdges);
    dir->Add(hPABCDraw);

    return dir;
}

template<typename T>
void FillABCDHistos(T obj, TDirectory *dir){
    if constexpr (std::is_same<T, std::vector<ABCD>>::value)
    {
        for(unsigned long i=0; i < obj.size(); i++)
        {
            ((TH1D*)dir->FindObject("hPABCD"))->Fill(i,obj.at(i).PABCD);
            cout<<obj.at(i).PABCD<<"\n";
            ((TH1D*)dir->FindObject("hPABCDraw"))->Fill(i,obj.at(i).PABCDraw);
            cout<<obj.at(i).PABCDraw<<"\n";
        }
    }else{

    } 
    return;
}

#endif // _ANALYSEHISTOSFROMTREE_

