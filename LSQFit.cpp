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
const int nexperiments=1000; // modify for more runs
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
  double lx2[npoints]; // will be lx * lx later on
  double ly[npoints];
  double ley[npoints];
  double a_parameters[nexperiments];
  double b_parameters[nexperiments];
  double c_parameters[nexperiments];
  double chi_parameters[nexperiments];
  double red_chi_parameters[nexperiments];

  getX(lx);
  getY(lx,ly,ley);
  for (int i = 0; i < npoints; i++) {
    lx2[i] = lx[i] * lx[i];
  }
  auto tgl = new TGraphErrors(npoints,lx,ly,0,ley);
  tgl->SetTitle("Pseudoexperiment;x;y");
  
  // An example of one pseudo experiment
  // tgl->Draw("alp");
  // tc->Draw();

  // here is my code!
  float function_points[npoints];
  float residuals[npoints];
  float residuals_squared[npoints];
  float chi_square = 0.; // placeholder
  float reduced_chi_square = 0.; // placeholder too
  TVectorD x; x.Use(npoints,lx);
  TVectorD x2; x2.Use(npoints,lx2); // will is lx * lx
  TVectorD y; y.Use(npoints,ly);
  TVectorD e; e.Use(npoints,ley);
  TMatrixD A(npoints,nPar); // A matrix
  for (int fit = 0; fit < nexperiments; fit++) {
    getX(lx);
    getY(lx, ly, ley);
    chi_square = 0.0; // need to reset it.
    for (int i = 0; i < npoints; i++) {
      function_points[i] = f(lx[i]);
      residuals[i] = ly[i] - function_points[i];
      residuals_squared[i] = residuals[i] * residuals[i];
      chi_square += residuals_squared[i] / (sigma + sigma);
    }
    reduced_chi_square = chi_square / (npoints - nPar);
    TMatrixDColumn(A,0) = 1.0; // parameter is a constant added
    TMatrixDColumn(A,1) = x;
    TMatrixDColumn(A,2) = x2;
    // cout << "A = ";
    // A.Print();

    // apply weights
    TMatrixD yw(A.GetNrows(),1);

    for (Int_t irow = 0; irow < A.GetNrows(); irow++) {
      TMatrixDRow(A,irow) *= 1/e(irow);
      TMatrixDRow(yw,irow) += y(irow)/e(irow);
    }
    // cout << "A weighted = ";
    // A.Print();
    // cout << "y weighted = ";
    // yw.Print();

    TMatrixD theta=SolveLSQ(A, yw);
    // cout << "Param vector = ";
    // theta.Print();

    a_parameters[fit] = theta[0][0];
    b_parameters[fit] = theta[1][0];
    c_parameters[fit] = theta[2][0];

    chi_parameters[fit] = chi_square;
    red_chi_parameters[fit] = reduced_chi_square;

    TF1 *fn1 = new TF1("fn1","[0] + [1]*x + [2]*x*x", xmin, xmax);

    fn1->SetParameters(theta[0][0], theta[1][0], theta[2][0]);
    tgl->Draw("alp*");
    fn1->Draw("same");
    // tc->Draw();

    // cout << "parameters";
  }

  // *** modify and add your code here ***

  TH2F *h1 = new TH2F("h1","Parameter b vs a;a;b",100,-0.5,0.5,100,0.7,0.9);
  TH2F *h2 = new TH2F("h2","Parameter c vs a;a;c",100,-0.5,0.5,100,-0.025,-0.015);
  TH2F *h3 = new TH2F("h3","Parameter c vs b;b;c",100,0.7,0.9,100,-0.025,-0.015);
  TH1F *h4 = new TH1F("h4","reduced chi^2;;frequency",100,0,4);

  TH1F *h5 = new TH1F("h5","Parameter a;;frequency",100,-1,1);
  TH1F *h6 = new TH1F("h6","Parameter b;;frequency",100,-1,1);
  TH1F *h7 = new TH1F("h7","Parameter c;;frequency",100,-1,1);
  TH1F *h8 = new TH1F("h8","chi^2;;frequency",100,0,npoints);

  // perform many least squares fits on different pseudo experiments here
  // I actually did that earlier. My data is in the a, b, c chi-square arrays up top.

  // fill histograms w/ required data
  for (int i = 0; i < nexperiments; i++) {
    // part 2 hists
    h1->Fill(a_parameters[i], b_parameters[i]);
    h2->Fill(a_parameters[i], c_parameters[i]);
    h3->Fill(b_parameters[i], c_parameters[i]);
    h4->Fill(red_chi_parameters[i]);
    // part 1 hists
    h5->Fill(a_parameters[i]);
    h6->Fill(b_parameters[i]);
    h7->Fill(c_parameters[i]);
    h8->Fill(chi_parameters[i]);
  }
  
  TCanvas *tc2 = new TCanvas("c2","my study results",200,200,dw,dh);
  tc2->Divide(2,2);
  tc2->cd(1); h1->Draw("colz");
  tc2->cd(2); h2->Draw("colz");
  tc2->cd(3); h3->Draw("colz");
  tc2->cd(4); h4->Draw();
  
  tc2->Draw();
  tc2->Print("test_plot.png");

  // **************************************
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
