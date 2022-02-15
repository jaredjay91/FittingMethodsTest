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

void doOneExp() {

  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");

  //define generating functions with 10 events per GeV background (100-1000 GeV, flat) and a Gaussian signal with 1000 events (average), with a mean of 700 and a sigma of 50 GeV.
  cout << "Defining generating functions..." << endl;
  TF1* sigfunc = new TF1("sigfunc","[2]/(sqrt(2*3.14159)*[1])*exp(-0.5*((x-[0])/[1])^2)",100,1000);
  sigfunc->SetParameters(700,50,1000);
  TF1* bkgfunc = new TF1("bkgfunc","10",100,1000);
  TF1 *genfunc = new TF1("genfunc","sigfunc+bkgfunc",100,1000);

  //Define parameters
  const int NEvents = 10000;
  int binwidth = 25;
  int min = 400;
  int max = 1000;
  int delta = max-min;
  const int numbins = delta/binwidth;
  float histmin = (float)min;
  float histmax = (float)max;

  //generate pseudodata to fill a histogram with 10k events
  TH1F *h1 = new TH1F("h1","Fitted Resonance Data",numbins,histmin,histmax);
  for (int iEvents = 0; iEvents<NEvents; iEvents++) {
    double value = genfunc->GetRandom();
    h1->Fill(value);
  }

  //draw the histogram
  TCanvas * c1 = new TCanvas("c1","c1",4,20,500,400);
  h1->Draw("E");
  h1->SetMarkerStyle(8);
  h1->GetYaxis()->SetTitle("count");
  h1->GetYaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetTitle("invariant mass");
  h1->GetXaxis()->SetTitleOffset(1.2);

  //Define arrays
  double xval, Nexp, Nobs;
  double sigmaGaus = 50;
  double ampGaus = 1/(sqrt(2*3.14159)*sigmaGaus)*binwidth;
  double muGaus = 700;
  double sParam[3] = {600,800,1000};
  double bParam[3] = {5,15,10};
  double negLogL[3] = {0,0,0};
  double minNegLogL = 1e6;
  int smallOne, bigOne, otherOne;
  double dev = (sParam[2]-sParam[0])*(sParam[2]-sParam[0]) + (bParam[2]-bParam[0])*(bParam[2]-bParam[0]);

  //Use converging triangle method to find minimum negative log likelihood
  TNtuple* ntuple = new TNtuple("ntuple","list of fitted values","s:b",1);
  TNtuple* ntuple1 = new TNtuple("ntuple1","list of fitted values","s:b",1);
  TNtuple* ntuple2 = new TNtuple("ntuple2","list of fitted values","s:b",1);
  TNtuple* ntuple3 = new TNtuple("ntuple3","list of fitted values","s:b",1);
  double pinch = 0.9;
  double maxOf3 = 0;
  double minOf3 = 0;
  while (dev>0.001) {
    smallOne = 0;
    otherOne = 0;
    bigOne = 0;
    maxOf3 = 0;
    minOf3 = minNegLogL;
    for (int j = 0; j<3; j++) {
      double logL = 0;
      for (int i = 0; i<numbins; i++) {
        xval = h1->GetBinCenter(i+1);
        Nobs = h1->GetBinContent(i+1);
        //cout << "Bin " << i+1 << ": (x,y) = (" << xval << "," << Nobs << ");" << endl;
        Nexp = (float)binwidth*bParam[j] + sParam[j]*ampGaus*exp(-0.5*pow((xval-muGaus)/sigmaGaus,2));
        logL += Nobs*(log(Nexp/Nobs)+1)-Nexp;
      }
      negLogL[j] = -logL;
      if (negLogL[j]<=minOf3) {
        smallOne = j;
        minOf3 = negLogL[j];
      }
      if (negLogL[j]>=maxOf3) {
        bigOne = j;
        maxOf3 = negLogL[j];
      }
      //cout << "(b,s) = " << "(" << bParam[j] << "," << sParam[j] << "): ";
      //cout << "neglogL = " << negLogL[j] << endl;
    }

    otherOne = 3-smallOne-bigOne;

    if (minOf3<minNegLogL) minNegLogL = minOf3;
    //cout << Form("%i,%i,%i",smallOne,otherOne,bigOne) << endl;
    double x1 = bParam[smallOne];
    double y1 = sParam[smallOne];
    double x2 = bParam[otherOne];
    double y2 = sParam[otherOne];
    double x3 = bParam[bigOne];
    double y3 = sParam[bigOne];
    double xnew = x3+pinch*(x1+x2-2*x3);
    double ynew = y3+pinch*(y1+y2-2*y3);
    dev = (y3-y1)*(y3-y1) + (x3-x1)*(x3-x1);
    bParam[bigOne] = xnew;
    sParam[bigOne] = ynew;
    //cout << "s = " << sParam[smallOne] << endl;
    //cout << "b = " << bParam[smallOne] << endl;
    //cout << "minNegLogL = " << minNegLogL << endl;
    ntuple->Fill(sParam[smallOne],bParam[smallOne]);
    ntuple1->Fill(sParam[0],bParam[0]);
    ntuple2->Fill(sParam[1],bParam[1]);
    ntuple3->Fill(sParam[2],bParam[2]);
  }//end of while loop

  double sFitted = sParam[smallOne];
  double bFitted = bParam[smallOne];
  cout << "(s,b) = (" << sFitted << "," << bFitted << ")" << endl;
  cout << "neglogL = " << minOf3 << endl;

  //Draw the fit function on the plot
  TF1* fitfunc = new TF1("fitfunc","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",min,max);
  fitfunc->SetParameters(sFitted*ampGaus,700,sigmaGaus,bFitted*binwidth);
  fitfunc->Draw("same");

  TPaveText *pt = new TPaveText(.1,.7,.3,.9,"brNDC");
  pt->AddText(Form("s = %.3f",sFitted));
  pt->AddText(Form("b = %.4f",bFitted));
  pt->Draw();

  c1->SaveAs("NLLFittedPseudoData.png");
  c1->SaveAs("NLLFittedPseudoData.pdf");

  //Draw path of parameters
  int Nsteps = (int)ntuple->GetEntries();
  const int constNsteps = Nsteps;
  float svals[constNsteps], bvals[constNsteps];
  float svals1[constNsteps], bvals1[constNsteps];
  float svals2[constNsteps], bvals2[constNsteps];
  float svals3[constNsteps], bvals3[constNsteps];
  cout << Nsteps << endl;
  TLeaf *sLeaf = ntuple->GetLeaf("s");
  TLeaf *bLeaf = ntuple->GetLeaf("b");
  TLeaf *sLeaf1 = ntuple1->GetLeaf("s");
  TLeaf *bLeaf1 = ntuple1->GetLeaf("b");
  TLeaf *sLeaf2 = ntuple2->GetLeaf("s");
  TLeaf *bLeaf2 = ntuple2->GetLeaf("b");
  TLeaf *sLeaf3 = ntuple3->GetLeaf("s");
  TLeaf *bLeaf3 = ntuple3->GetLeaf("b");
  for (int istep=0; istep<Nsteps; istep++) {
    ntuple->GetEntry(istep);
    svals[istep] = (float)sLeaf->GetValue();
    bvals[istep] = (float)bLeaf->GetValue();
    ntuple1->GetEntry(istep);
    svals1[istep] = (float)sLeaf1->GetValue();
    bvals1[istep] = (float)bLeaf1->GetValue();
    ntuple2->GetEntry(istep);
    svals2[istep] = (float)sLeaf2->GetValue();
    bvals2[istep] = (float)bLeaf2->GetValue();
    ntuple3->GetEntry(istep);
    svals3[istep] = (float)sLeaf3->GetValue();
    bvals3[istep] = (float)bLeaf3->GetValue();
  }

  TCanvas * c2 = new TCanvas("c2","c2",4,520,500,400);
  c2->cd();
  TGraph *gmin = new TGraph(Nsteps,bvals,svals);
  c2->SetLeftMargin(0.12);
  gmin->SetTitle("Path to minimum NLL");
  gmin->GetXaxis()->SetTitle("Background Parameter b");
  gmin->GetXaxis()->SetTitleOffset(1.2);
  gmin->GetYaxis()->SetTitle("Signal Parameter s");
  gmin->GetYaxis()->SetTitleOffset(1.8);
  gmin->Draw("AC*");
  c2->SaveAs("FitPath.pdf");

  TCanvas * c3 = new TCanvas("c3","c3",504,520,500,400);
  c3->cd();
  TGraph *g1 = new TGraph(Nsteps,bvals1,svals1);
  g1->SetLineColor(kRed);
  TGraph *g2 = new TGraph(Nsteps,bvals2,svals2);
  g2->SetLineColor(kBlue);
  TGraph *g3 = new TGraph(Nsteps,bvals3,svals3);
  g3->SetLineColor(kGreen);
  g1->Draw();
  g1->Draw("same*");
  g2->Draw("sameC*");
  g3->Draw("sameC*");

}
