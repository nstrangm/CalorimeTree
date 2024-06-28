#ifndef _ANALYSEHISTOSFROMTREE_
#define _ANALYSEHISTOSFROMTREE_

#include "Utilities.h"
#include "HistogramLibrary.h"
#include "PlottingClass.h"

const char* suffix="png";

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
    //PABCDalfa:
    Double_t PABCDalfa = 0;
    Double_t PABCDalfaErr = 0;
    //PABCD:
    Double_t PABCD = 0;
    Double_t PABCDErr = 0;
};

void CalculateABCDData(std::vector<ABCD> &ABCDs ,TDirectory* dGammaData, TDirectory* dGammaMC, TDirectory* ABCDdir, int nbins , float pTbinEdges[], float pTIsoABCDCuts[], float M02ABCDCuts[]){  
    //Load THnSparse:
    //Load MC data:
    THnSparse *PtIsovsM02vsPtMC = (THnSparse*)dGammaMC->Get("hIsoGammaIsovsM02vsPt");
    THnSparse *PtIsovsM02vsPtMCBackground = (THnSparse*)dGammaMC->Get("hIsoGammaIsovsM02vsPtBackground");
    //Load actual data:
    THnSparse *PtIsovsM02vsPtData = (THnSparse*)dGammaData->Get("hIsoGammaIsovsM02vsPt");

    //Loop over pT bins and calculate A,B,C and D:
    for(int bin = 1; bin<=nbins; ++bin){
        //ABCD object for storing output
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
        int pTIsoMinBin1 = h2M02vsPtprojectionMC->GetYaxis()->FindBin(pTIsoABCDCuts[0]);//Isolated min
        int pTIsoMaxBin1 = h2M02vsPtprojectionMC->GetYaxis()->FindBin(pTIsoABCDCuts[1]);//Isolated max
        int pTIsoMinBin2 = h2M02vsPtprojectionMC->GetYaxis()->FindBin(pTIsoABCDCuts[2]);//Antiisolated min
        int pTIsoMaxBin2 = h2M02vsPtprojectionMC->GetYaxis()->FindBin(pTIsoABCDCuts[3]);//Antiisolated max
        int M02MinBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindBin(M02ABCDCuts[0]);//Narrow min
        int M02MaxBin1 = h2M02vsPtprojectionMC->GetXaxis()->FindBin(M02ABCDCuts[1]);//Narrow max
        int M02MinBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindBin(M02ABCDCuts[2]);//wide min
        int M02MaxBin2 = h2M02vsPtprojectionMC->GetXaxis()->FindBin(M02ABCDCuts[3]);//wide max
        //Integrate MC-data
        ABCDtemp.AMC=h2M02vsPtprojectionMCBackground->IntegralAndError(M02MinBin1,M02MaxBin1-1,pTIsoMinBin1,pTIsoMaxBin1-1,ABCDtemp.AMCErr,"");
        ABCDtemp.BMC=h2M02vsPtprojectionMC->IntegralAndError(M02MinBin2,M02MaxBin2-1,pTIsoMinBin1,pTIsoMaxBin1-1,ABCDtemp.BMCErr,"");
        ABCDtemp.CMC=h2M02vsPtprojectionMC->IntegralAndError(M02MinBin1,M02MaxBin1-1,pTIsoMinBin2,pTIsoMaxBin2-1,ABCDtemp.CMCErr,"");
        ABCDtemp.DMC=h2M02vsPtprojectionMC->IntegralAndError(M02MinBin2,M02MaxBin2-1,pTIsoMinBin2,pTIsoMaxBin2-1,ABCDtemp.DMCErr,"");
        cout<<"MC:"<<"\n";
        cout<<"Integrateing MC background in M02: from "<<h2M02vsPtprojectionMCBackground->GetXaxis()->GetBinLowEdge(M02MinBin1)<<" To "<<h2M02vsPtprojectionMCBackground->GetXaxis()->GetBinUpEdge(M02MaxBin1-1)<<"\n";
        cout<<"Integrateing MC background in pTiso: from "<<h2M02vsPtprojectionMCBackground->GetYaxis()->GetBinLowEdge(pTIsoMinBin1)<<" To "<<h2M02vsPtprojectionMCBackground->GetYaxis()->GetBinUpEdge(pTIsoMaxBin1-1)<<"\n";
        cout<<"Yielded AMC="<<ABCDtemp.AMC<<"\n";

        cout<<"Integrateing MC in M02: from "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinLowEdge(M02MinBin2)<<" To "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinUpEdge(M02MaxBin2-1)<<"\n";
        cout<<"Integrateing MC in pTiso: from "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinLowEdge(pTIsoMinBin1)<<" To "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinUpEdge(pTIsoMaxBin1-1)<<"\n";
        cout<<"Yielded BMC="<<ABCDtemp.BMC<<"\n";

        cout<<"Integrateing MC in M02: from "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinLowEdge(M02MinBin1)<<" To "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinUpEdge(M02MaxBin1-1)<<"\n";
        cout<<"Integrateing MC in pTiso: from "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinLowEdge(pTIsoMinBin2)<<" To "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinUpEdge(pTIsoMaxBin2-1)<<"\n";
        cout<<"Yielded CMC="<<ABCDtemp.CMC<<"\n";

        cout<<"Integrateing MC in M02: from "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinLowEdge(M02MinBin2)<<" To "<<h2M02vsPtprojectionMC->GetXaxis()->GetBinUpEdge(M02MaxBin2-1)<<"\n";
        cout<<"Integrateing MC in pTiso: from "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinLowEdge(pTIsoMinBin2)<<" To "<<h2M02vsPtprojectionMC->GetYaxis()->GetBinUpEdge(pTIsoMaxBin2-1)<<"\n";
        cout<<"Yielded DMC="<<ABCDtemp.DMC<<"\n";


        ABCDdir->cd();
        TString M02vsPtIsoPTMCName=Form("M02vsPtIsoPTMC%f_%f",pTbinEdges[bin-1],pTbinEdges[bin]);
        h2M02vsPtprojectionMC->SetName(M02vsPtIsoPTMCName);
        h2M02vsPtprojectionMC->SetTitle(";#sigma_{long};p_{T}^{iso}[GeV/c]");
        h2M02vsPtprojectionMC->Write();

        TString M02vsPtIsoPTMCBAckgroundName=Form("M02vsPtIsoPTMCBackground%f_%f",pTbinEdges[bin-1],pTbinEdges[bin]);
        h2M02vsPtprojectionMCBackground->SetName(M02vsPtIsoPTMCBAckgroundName);
        h2M02vsPtprojectionMCBackground->SetTitle(";#sigma_{long};p_{T}^{iso}[GeV/c]");
        h2M02vsPtprojectionMCBackground->Write();



        //Freeing memory:
        delete h2M02vsPtprojectionMC;
        delete h2M02vsPtprojectionMCBackground;

        //For calculating the PABCD-raw:
        TH2F *h2M02vsPtprojectionData;
        //Restrict Pt-range
        PtIsovsM02vsPtData->GetAxis(2)->SetRangeUser(pTbinEdges[bin-1],pTbinEdges[bin]);
        //Project onto pt-range:
        h2M02vsPtprojectionData = (TH2F*)PtIsovsM02vsPtData->Projection(0,1);

        //Integrate
        ABCDtemp.AData=h2M02vsPtprojectionData->IntegralAndError(M02MinBin1,M02MaxBin1-1,pTIsoMinBin1,pTIsoMaxBin1-1,ABCDtemp.ADataErr,"");
        ABCDtemp.BData=h2M02vsPtprojectionData->IntegralAndError(M02MinBin2,M02MaxBin2-1,pTIsoMinBin1,pTIsoMaxBin1-1,ABCDtemp.BDataErr,"");
        ABCDtemp.CData=h2M02vsPtprojectionData->IntegralAndError(M02MinBin1,M02MaxBin1-1,pTIsoMinBin2,pTIsoMaxBin2-1,ABCDtemp.CDataErr,"");
        ABCDtemp.DData=h2M02vsPtprojectionData->IntegralAndError(M02MinBin2,M02MaxBin2-1,pTIsoMinBin2,pTIsoMaxBin2-1,ABCDtemp.DDataErr,"");
        
        cout<<"Data:"<<"\n";
        cout<<"Integrateing Data in M02: from "<<h2M02vsPtprojectionData->GetXaxis()->GetBinLowEdge(M02MinBin1)<<" To "<<h2M02vsPtprojectionData->GetXaxis()->GetBinUpEdge(M02MaxBin1-1)<<"\n";
        cout<<"Integrateing Data in pTiso: from "<<h2M02vsPtprojectionData->GetYaxis()->GetBinLowEdge(pTIsoMinBin1)<<" To "<<h2M02vsPtprojectionData->GetYaxis()->GetBinUpEdge(pTIsoMaxBin1-1)<<"\n";
        cout<<"Yielded AData="<<ABCDtemp.AData<<"\n";

        cout<<"Integrateing Data in M02: from "<<h2M02vsPtprojectionData->GetXaxis()->GetBinLowEdge(M02MinBin2)<<" To "<<h2M02vsPtprojectionData->GetXaxis()->GetBinUpEdge(M02MaxBin2-1)<<"\n";
        cout<<"Integrateing Data in pTiso: from "<<h2M02vsPtprojectionData->GetYaxis()->GetBinLowEdge(pTIsoMinBin1)<<" To "<<h2M02vsPtprojectionData->GetYaxis()->GetBinUpEdge(pTIsoMaxBin1-1)<<"\n";
        cout<<"Yielded BData="<<ABCDtemp.BData<<"\n";

        cout<<"Integrateing Data in M02: from "<<h2M02vsPtprojectionData->GetXaxis()->GetBinLowEdge(M02MinBin1)<<" To "<<h2M02vsPtprojectionData->GetXaxis()->GetBinUpEdge(M02MaxBin1-1)<<"\n";
        cout<<"Integrateing Data in pTiso: from "<<h2M02vsPtprojectionData->GetYaxis()->GetBinLowEdge(pTIsoMinBin2)<<" To "<<h2M02vsPtprojectionData->GetYaxis()->GetBinUpEdge(pTIsoMaxBin2-1)<<"\n";
        cout<<"Yielded CData="<<ABCDtemp.CData<<"\n";

        cout<<"Integrateing Data in M02: from "<<h2M02vsPtprojectionData->GetXaxis()->GetBinLowEdge(M02MinBin2)<<" To "<<h2M02vsPtprojectionData->GetXaxis()->GetBinUpEdge(M02MaxBin2-1)<<"\n";
        cout<<"Integrateing Data in pTiso: from "<<h2M02vsPtprojectionData->GetYaxis()->GetBinLowEdge(pTIsoMinBin2)<<" To "<<h2M02vsPtprojectionData->GetYaxis()->GetBinUpEdge(pTIsoMaxBin2-1)<<"\n";
        cout<<"Yielded DData="<<ABCDtemp.DData<<"\n";
        
        
        //Freeing memory:

        //Add distribution plot to putput:
        ABCDdir->cd();
        TString M02vsPtIsoPTName=Form("M02vsPtIsoPT%f_%f",pTbinEdges[bin-1],pTbinEdges[bin]);
        h2M02vsPtprojectionData->SetName(M02vsPtIsoPTName);
        h2M02vsPtprojectionData->SetTitle(";#sigma_{long};p_{T}^{iso}[GeV/c]");
        h2M02vsPtprojectionData->Write();

        


        delete h2M02vsPtprojectionData;


        //Error propagation(PABCDraw):
        float Atermraw=(ABCDtemp.CData*ABCDtemp.BData*ABCDtemp.ADataErr)/(ABCDtemp.DData*ABCDtemp.AData*ABCDtemp.AData);
        float Btermraw=(ABCDtemp.CData*ABCDtemp.BDataErr)/(ABCDtemp.DData*ABCDtemp.AData);
        float Ctermraw=(ABCDtemp.CDataErr*ABCDtemp.BData)/(ABCDtemp.DData*ABCDtemp.AData);
        float Dtermraw=(ABCDtemp.CData*ABCDtemp.BData*ABCDtemp.DDataErr)/(ABCDtemp.DData*ABCDtemp.DData*ABCDtemp.AData);
        float RDataerr=TMath::Sqrt(Atermraw*Atermraw+Btermraw*Btermraw+Ctermraw*Ctermraw+Dtermraw*Dtermraw);

        //Fill purity PABCDraw:
        float Rraw=(ABCDtemp.CData/ABCDtemp.AData)/(ABCDtemp.DData/ABCDtemp.BData);
        ABCDtemp.PABCDraw = 1-Rraw;
        ABCDtemp.PABCDrawErr = TMath::Sqrt(Atermraw*Atermraw+Btermraw*Btermraw+Ctermraw*Ctermraw+Dtermraw*Dtermraw);

        //Error propagation alfa:
        float AtermAlfa=(ABCDtemp.DMC*ABCDtemp.AMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC);
        float BtermAlfa=(ABCDtemp.AMC*ABCDtemp.DMC*ABCDtemp.BMCErr)/(ABCDtemp.BMC*ABCDtemp.BMC*ABCDtemp.CMC);
        float CtermAlfa=(ABCDtemp.AMC*ABCDtemp.DMC*ABCDtemp.CMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC*ABCDtemp.CMC);
        float DtermAlfa=(ABCDtemp.AMC*ABCDtemp.DMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC);

        //Calculate alfa:
        ABCDtemp.PABCDalfa= (ABCDtemp.AMC/ABCDtemp.CMC)/(ABCDtemp.BMC/ABCDtemp.DMC);
        ABCDtemp.PABCDalfaErr=TMath::Sqrt(AtermAlfa*AtermAlfa+BtermAlfa*BtermAlfa+CtermAlfa*CtermAlfa+DtermAlfa*DtermAlfa);

        //error propagation(PABCD):
        float RtermMC=(ABCDtemp.AMC*ABCDtemp.DMC*ABCDtemp.PABCDrawErr)/(ABCDtemp.BMC*ABCDtemp.CMC);
        float AtermMC=(Rraw*ABCDtemp.DMC*ABCDtemp.AMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC);
        float BtermMC=(Rraw*ABCDtemp.AMC*ABCDtemp.DMC*ABCDtemp.BMCErr)/(ABCDtemp.BMC*ABCDtemp.BMC*ABCDtemp.CMC);
        float CtermMC=(Rraw*ABCDtemp.AMC*ABCDtemp.DMC*ABCDtemp.CMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC*ABCDtemp.CMC);
        float DtermMC=(Rraw*ABCDtemp.AMC*ABCDtemp.DMCErr)/(ABCDtemp.BMC*ABCDtemp.CMC);

        //PABCD:
        ABCDtemp.PABCD = 1-Rraw*ABCDtemp.PABCDalfa;
        ABCDtemp.PABCDErr = TMath::Sqrt(RtermMC*RtermMC+AtermMC*AtermMC+BtermMC*BtermMC+CtermMC*CtermMC+DtermMC*DtermMC);

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
    TH1D *hPABCDalfa=new TH1D("hPABCDalfa","#alpha_{ABCD}",nbins, pTbinEdges);
    dir->Add(hPABCDalfa);

    return dir;
}

template<typename T>
void FillABCDHistos(T obj, TDirectory *dir){
    if constexpr (std::is_same<T, std::vector<ABCD>>::value)
    {
        for(unsigned long i=0; i < obj.size(); i++)
        {
            ((TH1D*)dir->FindObject("hPABCD"))->Fill(pTbinsToplot[i],obj.at(i).PABCD);
            ((TH1D*)dir->FindObject("hPABCD"))->SetBinError(i+1,obj.at(i).PABCDErr);
            ((TH1D*)dir->FindObject("hPABCDraw"))->Fill(pTbinsToplot[i],obj.at(i).PABCDraw);
            ((TH1D*)dir->FindObject("hPABCDraw"))->SetBinError(i+1,obj.at(i).PABCDrawErr);
            ((TH1D*)dir->FindObject("hPABCDalfa"))->Fill(pTbinsToplot[i],obj.at(i).PABCDalfa);
            ((TH1D*)dir->FindObject("hPABCDalfa"))->SetBinError(i+1,obj.at(i).PABCDalfaErr);
        }

    }else{

    } 
    return;
}

#endif // _ANALYSEHISTOSFROMTREE_

