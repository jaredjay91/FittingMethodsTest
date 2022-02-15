//1. Bias study in likelihood fitting.  Generate your own pseudodata outcomes, 10k events per sample, with 10 events per GeV background (100-1000 GeV, flat) and a Gaussian signal with 1000 events (average), with a mean of 700 and a sigma of 50 GeV.  Use your existing binned likelihood fitting code to fit each pseudoexperiment and plot the deviation of each fit quantity from the true value.

#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "time.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "FitMethods.h"

void FinalProject() {

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  //gStyle->SetLabelSize(0.05,"X");
  //gStyle->SetLabelSize(0.05,"Y");

  //define generating functions with 10 events per GeV background (100-1000 GeV, flat) and a Gaussian signal with 1000 events (average), with a mean of 700 and a sigma of 50 GeV.
  cout << "Defining generating functions..." << endl;
  TF1* sigfunc = new TF1("sigfunc","[2]/(sqrt(2*3.14159)*[1])*exp(-0.5*((x-[0])/[1])^2)",100,1000);
  sigfunc->SetParameters(700,50,1000);
  TF1* bkgfunc = new TF1("bkgfunc","10",100,1000);
  TF1 *genfunc = new TF1("genfunc","sigfunc+bkgfunc",100,1000);

  //Define parameters
  const int NTrials = 100;
  const int NEvents = 10000;
  const int reportEvery = (int)(NTrials/100.0);
  int binwidth = 36;
  int min = 100;
  int max = 1000;
  int delta = max-min;
  const int numbins = delta/binwidth;
  float histmin = (float)min;
  float histmax = (float)max;

  TString outFileName = "PseudoExpResults_FullRange.root";
  TFile outfile (outFileName, "RECREATE");
  TNtuple* ntuple = new TNtuple("ntuple","results of pseudoexperiments","s:b",NTrials);

  //generate pseudodata to fill a histogram with 10k events
  for (int iTrials = 0; iTrials<NTrials; iTrials++) {
    if (iTrials%reportEvery==0) cout << "Completed trials: " << iTrials << "/" << NTrials << " (" << (float)iTrials/(float)NTrials*100 << "%)" << endl;
    TH1F *h1 = new TH1F("h1","Fitted Resonance Data",numbins,histmin,histmax);
    for (int iEvents = 0; iEvents<NEvents; iEvents++) {
      double value = genfunc->GetRandom();
      h1->Fill(value);
    }

    //draw the histogram
    /*TCanvas * c1 = new TCanvas("c1","c1",4,20,500,400);
    h1->Draw("E");
    h1->SetMarkerStyle(8);
    h1->GetYaxis()->SetTitle("count");
    h1->GetYaxis()->SetTitleOffset(1.5);
    h1->GetXaxis()->SetTitle("invariant mass");
    h1->GetXaxis()->SetTitleOffset(1.2);*/

    double sGuess=1000;
    double bGuess=10;
    double sFit;
    double* sFitptr;
    sFitptr = &sFit;
    double bFit;
    double* bFitptr;
    bFitptr = &bFit;

    simplexFit(h1,sGuess,bGuess,sFitptr,bFitptr);

    double sFitted = *sFitptr;
    double bFitted = *bFitptr;
    ntuple->Fill(sFitted,bFitted);
    //cout << "Trial " << iTrials+1 << ": " << "(s,b) = (" << sFitted << "," << bFitted << ")" << endl;
    delete h1;
  }//end of trials loop

  cout << "Plotting fit results..." << endl;
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(2,1);
  cntuple->cd(1);
  TH1F* shisto = new TH1F("shisto", "deviation from s",100,800,1200);
  ntuple->Draw("s>>shisto");
  cntuple->cd(2);
  TH1F* bhisto = new TH1F("bhisto", "deviation from b",100,8,12);
  ntuple->Draw("b>>bhisto");

  cout << "Writing to file..." << endl;
  cntuple->SaveAs("PseudoExpResults_FullRange.png");
  cntuple->SaveAs("PseudoExpResults_FullRange.pdf");
  ntuple->Write();
  cntuple->Write();
  outfile.Close();

  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;
}
