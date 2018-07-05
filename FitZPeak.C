#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "BW.C"

void FitZPeak(){

  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptFit(111111);

  TFile *f = new TFile("myZPeakCRAB.root","READ");

  TH1F *Z_mass= (TH1F*)f->Get("analyzeBasicPat/mumuMass");

  float massMIN = 88.0; // 80.0
  float massMAX = 94.0; // 100.0


  ////////////////
  //For Gaussian//
  ///////////////
  /*
  TF1 *func = new TF1("mygauss",mygauss,massMIN, massMAX,3); 
  func->SetParameter(0,1.0);   func->SetParName(0,"const");  
  func->SetParameter(2,5.0);   func->SetParName(2,"sigma");  
  func->SetParameter(1,95.0);     func->SetParName(1,"mean");

  Z_mass->Fit("mygauss","QR");
  TF1 *fit = Z_mass->GetFunction("mygauss");
*/
  /////////////////////
  // For Breit-Wigner//
  ////////////////////
  
     TF1 *func = new TF1("mybw",mybw,massMIN, massMAX,3);
     func->SetParameter(0,1.0);   func->SetParName(0,"const");
     func->SetParameter(2,5.0);     func->SetParName(1,"sigma");
     func->SetParameter(1,95.0);    func->SetParName(2,"mean");

     Z_mass->Fit("mybw","QR");
     TF1 *fit = Z_mass->GetFunction("mybw");
  




  fit->SetLineColor(4);
  fit->SetLineWidth(3);

  Z_mass->SetMarkerStyle(20);
  Z_mass->SetMarkerSize(1.0);
  Z_mass->SetMarkerColor(2);
  Z_mass->SetLineWidth(2.0);
  Z_mass->SetXTitle("Z Mass (in GeV/c^{2})"); 
  // Z_mass->GetXaxis()->SetTitleOffset(1.0);


  Z_mass->Draw("PE0");
  // c1->Print("myZmass_BWfitted.eps");
  // c1->Print ("myZmass_Gausfitted.eps");
}
