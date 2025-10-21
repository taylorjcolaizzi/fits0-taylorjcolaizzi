#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGClient.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"


#include <iostream>
using namespace std;

using TMath::Log;

//parms
const double xmin=1;
const double xmax=20;
const int npoints=12; // modify for more points
const double sigma=0.2;
const int nexperiments=3; // modify for more runs
const int nPar=3; // tied to number of parameters below

double f(double x){
  const double a=0.5;
  const double b=1.3;
  const double c=0.5;
  return a+b*Log(x)+c*Log(x)*Log(x);
}

void getX(double *x){
  double step=(xmax-xmin)/npoints;
  for (int i=0; i<npoints; i++){
    x[i]=xmin+i*step;
  }
}

void getY(const double *x, double *y, double *ey){
  static TRandom2 tr(0);
  for (int i=0; i<npoints; i++){
    y[i]=f(x[i])+tr.Gaus(0,sigma);
    ey[i]=sigma;
  }
}


void leastsq(){
  double x[npoints];
  double y[npoints];
  double ey[npoints];
  getX(x);
  getY(x,y,ey);
  auto tg = new TGraphErrors(npoints,x,y,0,ey);
  tg->Draw("alp");
}

TMatrixD SolveLSQ(const TMatrixD &A, const TMatrixD &y){
  TMatrixD AT=(A);
  AT.T(); // A transpose
  TMatrixD ATAi(AT,TMatrixD::kMult,A);
  ATAi.Invert();
  TMatrixD Adag(ATAi,TMatrixD::kMult,AT); // thingy
  TMatrixD theta(Adag,TMatrixD::kMult,y);
  return theta;
}

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  gStyle->SetOptStat(0); // turn off histogram stats box


  TCanvas *tc = new TCanvas("c1","Sample dataset",dw,dh);

  double lx[npoints];
  double ly[npoints];
  double ley[npoints];

  getX(lx);
  getY(lx,ly,ley);
  auto tgl = new TGraphErrors(npoints,lx,ly,0,ley);
  tgl->SetTitle("Pseudoexperiment;x;y");
  
  // An example of one pseudo experiment
  tgl->Draw("alp");
  tc->Draw();

  // here is my code!
  float function_points[npoints];
  float residuals[npoints];
  float residuals_squared[npoints];
  float chi_square = 0.; // placeholder
  float reduced_chi_square = 0.; // placeholder too
  TVectorD x; x.Use(npoints,lx);
  TVectorD x2; x2.Use(npoints,lx); // will be squaring this.
  for (int index = 0; index < npoints; index++) {
    x2(index) = x2(index) * x2(index);
  }
  TVectorD y; y.Use(npoints,ly);
  TVectorD e; e.Use(npoints,ley);
  TMatrixD A(npoints,nPar); // A matrix
  for (int fit = 0; fit < nexperiments; fit++) {
    getX(lx);
    getY(lx, ly, ley);
    for (int i = 0; i < npoints; i++) {
      function_points[i] = f(lx[i]);
      residuals[i] = ly[i] - function_points[i];
      residuals_squared[i] = residuals[i] * residuals[i];
      chi_square += residuals_squared[i] / (sigma + sigma);
    }
    TMatrixDColumn(A,0) = 1.0; // parameter is a constant added
    TMatrixDColumn(A,1) = x;
    TMatrixDColumn(A,2) = x2;
    cout << "A = ";
    A.Print();
  }


  
  // *** modify and add your code here ***

  TH2F *h1 = new TH2F("h1","Parameter b vs a;a;b",100,0,1,100,0,1);
  TH2F *h2 = new TH2F("h2","Parameter c vs a;a;c",100,0,1,100,0,1);
  TH2F *h3 = new TH2F("h3","Parameter c vs b;b;c",100,0,1,100,0,1);
  TH1F *h4 = new TH1F("h4","reduced chi^2;;frequency",100,0,1);

  // perform many least squares fits on different pseudo experiments here
  // fill histograms w/ required data
  
  TCanvas *tc2 = new TCanvas("c2","my study results",200,200,dw,dh);
  tc2->Divide(2,2);
  tc2->cd(1); h1->Draw("colz");
  tc2->cd(2); h2->Draw("colz");
  tc2->cd(3); h3->Draw("colz");
  tc2->cd(4); h4->Draw();
  
  tc2->Draw();

  // **************************************
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
