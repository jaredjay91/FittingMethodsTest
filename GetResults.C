

void GetResults() {

  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05);
  gStyle->SetStatW(0.3);
  gStyle->SetStatW(0.3);

  TString fileName = "LLSPseudoExpResults.root";
  TFile* theFile = new TFile(fileName);
  TNtuple* ntupleResults = (TNtuple*)theFile->Get("ntuple");

  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,45,1000,400);
  cntuple->Divide(2,1);
  cntuple->cd(1);
  TH1F* shisto = new TH1F("shisto", "distribution of s",100,800,1200);
  ntupleResults->Draw("s>>shisto");
  cntuple->cd(2);
  TH1F* bhisto = new TH1F("bhisto", "distribution of b",100,9,11);
  ntupleResults->Draw("b>>bhisto");

  TCanvas* cpull =  new TCanvas("cpull","pullresults",4,545,1000,400);
  cpull->Divide(2,1);
  TH1F* spull = new TH1F("spull", "pull from s",100,-6,6);
  TH1F* bpull = new TH1F("bpull", "pull from b",100,-0.15,0.15);

  shisto->GetXaxis()->SetTitle("s");
  shisto->GetXaxis()->SetTitleOffset(1.0);
  bhisto->GetXaxis()->SetTitle("b");
  bhisto->GetXaxis()->SetTitleOffset(1.0);
  spull->GetXaxis()->SetTitle("pull");
  spull->GetXaxis()->SetTitleOffset(1.0);
  bpull->GetXaxis()->SetTitle("pull");
  bpull->GetXaxis()->SetTitleOffset(1.0);

  TLeaf *sLeaf = ntupleResults->GetLeaf("s");
  TLeaf *bLeaf = ntupleResults->GetLeaf("b");
  const int NTrials = (int)ntupleResults->GetEntries();
  for (int iTrial=0; iTrial<NTrials; iTrial++) {
    ntuple->GetEntry(iTrial);
    float sval = (float)sLeaf->GetValue();
    float bval = (float)bLeaf->GetValue();
    float serr = sqrt(sval);
    float berr = sqrt(bval);
    spull->Fill((sval-1000)/serr);
    bpull->Fill((bval-10)/berr);
  }

  cpull->cd(1);
  spull->Draw();
  cpull->cd(2);
  bpull->Draw();

  cntuple->SaveAs("LLSDistributions.png");
  cntuple->SaveAs("LLSDistributions.pdf");
  cpull->SaveAs("LLSPulls.png");
  cpull->SaveAs("LLSPulls.pdf");

}
