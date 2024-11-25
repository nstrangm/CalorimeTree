void eventsorting(TString inputfile = "input.root",TString outputfile = "output.root", TString cut = "event_centrality>=0 && event_centrality<10"){
    // open file with inputfile
    TFile *file = TFile::Open(inputfile, "READ");
    // get Tree called eventTree
    TTree *eventTree = (TTree *)file->Get("eventTree");
    TFile *outFile = new TFile(outputfile, "RECREATE");
    outFile->cd();
    TTree* outTree = eventTree->CopyTree(cut);
    outTree->Write();
    outFile->Close();
}