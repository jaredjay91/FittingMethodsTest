void simplexFit(TH1F* h1, double sGuess=1000, double bGuess=10, double* sFitptr=0, double* bFitptr=0) {
    //Use converging triangle method (simplex) to find minimum negative log likelihood

    double sParam[3] = {sGuess*0.6,sGuess,sGuess*1.4}; //3 initial points
    double bParam[3] = {bGuess*0.6,bGuess*1.4,bGuess};
    int smallOne, bigOne, otherOne; //indices

    double negLogL[3] = {0,0,0};

    double minNegLogL = 1e6; //This will decrease over time as the fitter searches for the minimum.

    double dev = (sParam[2]-sParam[0])*(sParam[2]-sParam[0]) + (bParam[2]-bParam[0])*(bParam[2]-bParam[0]); //Deviation is sum of squares. No square roots here because they take too long.

    double pinch = 0.9;//Determines how quickly the triangle shrinks and converges. If you converge too slowly it will take a long time. If you converge too quickly you will miss the true minimum.

    double tolerance = 0.01;//the largest deviation from the truth that you're willing to tolerate.

    int numbins = h1->GetXaxis()->GetNbins();
    //cout << "numbins = " << numbins << endl;
    float binwidth = h1->GetBinWidth(1);
    //cout << "binwidth = " << binwidth << endl;

    double xval, Nexp, Nobs;
    double sigmaGaus = 50;
    double ampGaus = 1/(sqrt(2*3.14159)*sigmaGaus)*binwidth;
    double muGaus = 700;

    while (dev>tolerance) {
      smallOne = 0;
      otherOne = 0;
      bigOne = 0;
      double maxOf3 = 0;
      double minOf3 = minNegLogL;
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
    }//end of while loop
    *sFitptr = sParam[smallOne];
    *bFitptr = bParam[smallOne];
}

