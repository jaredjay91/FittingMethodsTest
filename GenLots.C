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

void GenLots() {

  //gStyle->SetErrorX(0);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");

  //define generating functions with 10 events per GeV background (100-1000 GeV, flat) and a Gaussian signal with 1000 events (average), with a mean of 700 and a sigma of 50 GeV.
  TF1* sigfunc = new TF1("sigfunc","[2]/(sqrt(2*3.14159)*[1])*exp(-0.5*((x-[0])/[1])^2)",100,1000);
  sigfunc->SetParameters(700,50,1000);
  TF1* bkgfunc = new TF1("bkgfunc","10",100,1000);
  TF1 *genfunc = new TF1("genfunc","sigfunc+bkgfunc",100,1000);

  //Define parameters
  const int NTrials = 1e6;
  const int NEvents = 1e7;
  int binwidth = 25;
  int min = 400;
  int max = 1000;
  int delta = max-min;
  const int numbins = delta/binwidth;
  float histmin = (float)min;
  float histmax = (float)max;

  //generate pseudodata to fill a histogram with 10k events
  TH1F *h1 = new TH1F("h1","PseudoData",numbins,histmin,histmax);
  for (int iEvents = 0; iEvents<NEvents; iEvents++) {
    double value = genfunc->GetRandom();
    h1->Fill(value);
  }

  //draw the histogram
  TCanvas * c1 = new TCanvas("c1","c1",4,20,1000,400);
  c1->Divide(2,1);
  c1->cd(1);
  h1->Draw("E");
  h1->SetMarkerStyle(8);
  h1->GetYaxis()->SetTitle("count");
  h1->GetYaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetTitle("invariant mass");
  h1->GetXaxis()->SetTitleOffset(1.2);

  //Draw the fit function on the plot
  double sigmaGaus = 50;
  double ampGaus = 1/(sqrt(2*3.14159)*sigmaGaus)*binwidth;
  double muGaus = 700;
  TF1* fitfunc = new TF1("fitfunc","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",min,max);
  fitfunc->SetParameters(NEvents*ampGaus/10,700,sigmaGaus,NEvents*binwidth/1000);
  fitfunc->Draw("same");

  //make a pull histogram
  TH1F *h2 = new TH1F("h2","PseudoData Pull",numbins,histmin,histmax);
  double avgpull=0;
  for (int iBins = 1; iBins<=numbins; iBins++) {
    double Nobs = h1->GetBinContent(iBins);
    double xval = h1->GetBinCenter(iBins);
    double Nexp = NEvents*ampGaus/10*exp(-0.5*pow((xval-700)/sigmaGaus,2))+NEvents*binwidth/1000;
    double diff = (Nobs-Nexp)/sqrt(Nexp);
    avgpull += diff;
    h2->SetBinContent(iBins,diff);
  }
  avgpull = avgpull/numbins;//NEvents;
  cout << "avgpull = " << avgpull << endl;

  //draw the histogram
  //TCanvas * c2 = new TCanvas("c2","c2",504,20,500,400);
  c1->cd(2);
  h2->Draw();
  h2->SetMarkerStyle(8);
  h2->SetStats(kFALSE);
  h2->Sumw2();
  //h2->Scale(1/sqrt(NEvents));
  //h2->GetYaxis()->SetRangeUser(-13*NEvents/10000,13*NEvents/10000);
  h2->GetYaxis()->SetTitle("pull");
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetXaxis()->SetTitle("invariant mass");
  h2->GetXaxis()->SetTitleOffset(1.2);
  TLine *zeroline = new TLine(histmin,0,histmax,0);
  zeroline->SetLineColor(kRed);
  zeroline->Draw();

  c1->SaveAs("LotsOfPseudoData.png");
  c1->SaveAs("LotsOfPseudoData.pdf");
}
