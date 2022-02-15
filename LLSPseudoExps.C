#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "time.h"
#include "TNtuple.h"
#include "TCanvas.h"

void LLSPseudoExps() {

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  //gStyle->SetLabelSize(0.05,"X");
  //gStyle->SetLabelSize(0.05,"Y");

  //define generating functions with 10 events per GeV background (100-1000 GeV, flat) and a Gaussian signal with 1000 events (average), with a mean of 700 and a sigma of 50 GeV.
  TF1* sigfunc = new TF1("sigfunc","[2]/(sqrt(2*3.14159)*[1])*exp(-0.5*((x-[0])/[1])^2)",100,1000);
  sigfunc->SetParameters(700,50,1000);
  TF1* bkgfunc = new TF1("bkgfunc","10",100,1000);
  TF1 *genfunc = new TF1("genfunc","sigfunc+bkgfunc",100,1000);

  //Define parameters
  const int NTrials = 1e6;
  const int NEvents = 10000;
  int binwidth = 25;
  int min = 400;
  int max = 1000;
  int delta = max-min;
  const int numbins = delta/binwidth;
  float histmin = (float)min;
  float histmax = (float)max;

  TString outFileName = "LLSPseudoExpResults.root";
  TFile outfile (outFileName, "RECREATE");
  TNtuple* ntuple = new TNtuple("ntuple","results of pseudoexperiments","s:b:serr:berr",NTrials);

  //generate pseudodata to fill a histogram with 10k events
  for (int iTrials = 0; iTrials<NTrials; iTrials++) {
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


    //Define arrays
    const int numparams = 2;
    float xArray[numbins], sigma[numbins], yMeasured[numbins];
    float func[numparams][numbins];
    float sigmaGaus = 50;
    float ampGaus = 1/(sqrt(2*3.14159)*sigmaGaus)*binwidth;
    float muGaus = 700;
    for (int i = 0; i<numbins; i++) {
      xArray[i] = h1->GetBinCenter(i+1);
      yMeasured[i] = h1->GetBinContent(i+1);
      //cout << "Bin " << i+1 << ": (x,y) = (" << xArray[i] << "," << yMeasured[i] << ");" << endl;
      func[0][i] = (float)binwidth;
      func[1][i] = ampGaus*exp(-0.5*pow((xArray[i]-muGaus)/sigmaGaus,2));
      sigma[i] = sqrt(xArray[i]);
    }

    //Calculate y vector
    float yVector[numparams] = {0};
    for (int m = 0; m<numparams; m++) {
      for (int n = 0; n<numbins; n++) {
        yVector[m] += func[m][n]*yMeasured[n]/pow(sigma[n],2);
      }
    }

    //Calculate covariance matrix
    float V[numparams][numparams] = {0};
    //cout << endl << "VARIANCE MATRIX" << endl;
    for (int j=0; j<numparams; j++) {
      for (int k=0; k<numparams; k++) {
        for (int l=0; l<numbins; l++) {
          V[j][k] += func[j][l]*func[k][l]/pow(sigma[l],2);
        }
        //cout << V[j][k] << "    ";
      }
      //cout << endl;
    }
  
    //Calculate V-inverse
    float Vinv[numparams][numparams];
    float det = V[0][0]*V[1][1]-V[0][1]*V[1][0];
    Vinv[0][0] = V[1][1]/det;
    Vinv[0][1] = -V[0][1]/det;
    Vinv[1][0] = -V[1][0]/det;
    Vinv[1][1] = V[0][0]/det;
    //cout << endl << "COVARIANCE MATRIX" << endl;
    for (int j=0; j<numparams; j++) {
      for (int k=0; k<numparams; k++) {
        //cout << Vinv[j][k] << "    ";
      }
      //cout << endl;
    }
    //cout << endl;

    //calculate fitted parameters and errors
    float a[numparams];
    float aerr[numparams];
    a[0] = Vinv[0][0]*yVector[0] + Vinv[0][1]*yVector[1];
    a[1] = Vinv[1][0]*yVector[0] + Vinv[1][1]*yVector[1];
    aerr[0] = sqrt(Vinv[0][0]);
    aerr[1] = sqrt(Vinv[1][1]);
    //cout << "a[0] = " << a[0] << " +- " << aerr[0] << endl;
    //cout << "a[1] = " << a[1] << " +- " << aerr[1] << endl;

    //calculate correlation coefficient
    float rho = Vinv[1][0]/(aerr[0]*aerr[1]);
    //cout << "rho = " << rho << endl;

    ntuple->Fill(a[1],a[0],aerr[1],aerr[0]);
    //cout << "Trial " << iTrials+1 << ": " << "(s,b) = (" << a[1] << "," << a[0] << ")" << endl;
    delete h1;
  }//end of trials loop

  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(2,1);
  cntuple->cd(1);
  TH1F* shisto = new TH1F("shisto", "deviation from s",100,min,max);
  ntuple->Draw("s>>shisto");
  cntuple->cd(2);
  TH1F* bhisto = new TH1F("bhisto", "deviation from b",100,min,max);
  ntuple->Draw("b>>bhisto");

  ntuple->Write();
  outfile.Close();

  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;
}
