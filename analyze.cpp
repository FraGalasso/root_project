#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

void analyze() {
  // rootfile
  TFile *f = new TFile("analysis.root");

  // histos from rootfile

  // particles generated
  TH1F *histoparticles = (TH1F *)f->Get("histoparticles");
  std::cout << "histoparticles: " << histoparticles->GetEntries() << '\n';
  // phi distribution
  TH1F *phi_distrib = (TH1F *)f->Get("phi_distrib");
  std::cout << "phi_distrib: " << phi_distrib->GetEntries() << '\n';
  // theta distribution
  TH1F *theta_distrib = (TH1F *)f->Get("theta_distrib");
  std::cout << "theta_distrib: " << theta_distrib->GetEntries() << '\n';
  // momentum distribution
  TH1F *momentum_distrib = (TH1F *)f->Get("momentum_distrib");
  std::cout << "momentum_distrib: " << momentum_distrib->GetEntries() << '\n';
  // transversal momentum distribution
  TH1F *trans_momentum_distrib = (TH1F *)f->Get("trans_momentum_distrib");
  std::cout << "trans_momentum_distrib: "
            << trans_momentum_distrib->GetEntries() << '\n';
  // energy distribution
  TH1F *energy_distrib = (TH1F *)f->Get("energy_distrib");
  std::cout << "energy_distrib: " << energy_distrib->GetEntries() << '\n';
  // invariant mass distribution
  TH1F *inv_mass = (TH1F *)f->Get("inv_mass");
  std::cout << "inv_mass: " << inv_mass->GetEntries() << '\n';
  // different charge invariant mass distribution
  TH1F *diff_inv_mass = (TH1F *)f->Get("diff_inv_mass");
  std::cout << "diff_inv_mass: " << diff_inv_mass->GetEntries() << '\n';
  // same charge invariant mass distribution
  TH1F *same_inv_mass = (TH1F *)f->Get("same_inv_mass");
  std::cout << "same_inv_mass: " << same_inv_mass->GetEntries() << '\n';
  // different charge kaons and pions invariant mass distribution
  TH1F *diff_inv_mass_ka_pi = (TH1F *)f->Get("diff_inv_mass_ka_pi");
  std::cout << "diff_inv_mass_ka_pi: " << diff_inv_mass_ka_pi->GetEntries()
            << '\n';
  // same charge kaons and pions invariant mass distribution
  TH1F *same_inv_mass_ka_pi = (TH1F *)f->Get("same_inv_mass_ka_pi");
  std::cout << "same_inv_mass_ka_pi: " << same_inv_mass_ka_pi->GetEntries()
            << '\n';
  // generated particles invariant mass distribution
  TH1F *k_inv_mass = (TH1F *)f->Get("k_inv_mass");
  std::cout << "k_inv_mass: " << k_inv_mass->GetEntries() << '\n';

  TCanvas *can = new TCanvas("can", "Many Histos", 200, 10, 1200, 800);

  // massive histogram cosmetics session

  histoparticles->GetXaxis()->SetBinLabel(1, "pion+");
  histoparticles->GetXaxis()->SetBinLabel(2, "pion-");
  histoparticles->GetXaxis()->SetBinLabel(3, "kaon+");
  histoparticles->GetXaxis()->SetBinLabel(4, "kaon-");
  histoparticles->GetXaxis()->SetBinLabel(5, "proton");
  histoparticles->GetXaxis()->SetBinLabel(6, "antiproton");
  histoparticles->GetXaxis()->SetBinLabel(7, "k*");
  histoparticles->GetXaxis()->SetTitle("Types of particles");
  histoparticles->GetYaxis()->SetTitle("Occurrences");
  // checking particle types proportions, expected (with equal distribution of
  // charge): 80% pions, 10% kaons, 9% protons, 1% k*(resonance)
  for (int i = 0; i != 7; ++i) {
    std::cout << histoparticles->GetXaxis()->GetBinLabel(i + 1)
              << " generated: ";
    std::cout << histoparticles->GetBinContent(i + 1) << " +/- ";
    std::cout << histoparticles->GetBinError(i + 1) << '\n';
  }
  histoparticles->SetFillColor(kBlue);

  phi_distrib->GetXaxis()->SetTitle("Phi distribution");
  phi_distrib->GetYaxis()->SetTitle("Occurrences");
  phi_distrib->SetLineColor(kBlue);
  phi_distrib->Fit("pol0", "Q");
  TF1 *fun1 = phi_distrib->GetFunction("pol0");
  fun1->SetLineColor(kRed);
  fun1->SetLineWidth(2);
  std::cout << "\nFitting phi distribution:\nParameters: "
            << fun1->GetParameter(0) << " +/- " << fun1->GetParError(0) << '\n';

  theta_distrib->GetXaxis()->SetTitle("Theta distribution");
  theta_distrib->GetYaxis()->SetTitle("Occurrences");
  theta_distrib->SetLineColor(kBlue);
  theta_distrib->Fit("pol0", "Q");
  TF1 *fun2 = theta_distrib->GetFunction("pol0");
  fun2->SetLineColor(kRed);
  fun2->SetLineWidth(2);
  std::cout << "\nFitting theta distribution:\nParameter: "
            << fun2->GetParameter(0) << " +/- " << fun2->GetParError(0)
            << "\n\n";

  momentum_distrib->GetXaxis()->SetTitle("Momentum distribution");
  momentum_distrib->GetYaxis()->SetTitle("Occurrences");
  momentum_distrib->SetLineColor(kBlue);
  momentum_distrib->SetFillColor(kCyan);
  momentum_distrib->Fit("expo", "Q");
  TF1 *fun3 = momentum_distrib->GetFunction("expo");
  fun3->SetLineColor(kRed);
  fun3->SetLineWidth(2);
  // m = 1 / (-p1)
  double fit_mean = 1 / -(fun3->GetParameter(1));
  // dm = (dp1 / (-p1)) * m
  double fit_mean_err =
      (fun3->GetParError(1) / -(fun3->GetParameter(1))) * fit_mean;
  std::cout << "Fitting momentum distribution:\n"
            << "Mean from the histogram: " << momentum_distrib->GetMean()
            << " +/- " << momentum_distrib->GetMeanError() << '\n'
            << "Mean from the fit: " << fit_mean << " +/- " << fit_mean_err
            << "\n\n";

  trans_momentum_distrib->GetXaxis()->SetTitle(
      "Transverse momentum distribution");
  trans_momentum_distrib->GetYaxis()->SetTitle("Occurrences");
  trans_momentum_distrib->SetLineColor(kBlue);
  trans_momentum_distrib->SetFillColor(kCyan);

  energy_distrib->GetXaxis()->SetTitle("Energy distribution");
  energy_distrib->GetYaxis()->SetTitle("Occurrences");
  energy_distrib->SetLineColor(kBlue);
  energy_distrib->SetFillColor(kCyan);

  inv_mass->GetXaxis()->SetTitle("Invariant mass");
  inv_mass->GetYaxis()->SetTitle("Occurrences");
  inv_mass->SetLineColor(kBlue);
  inv_mass->SetFillColor(kCyan);

  diff_inv_mass->GetXaxis()->SetTitle(
      "Invariant mass of particles with different charge");
  diff_inv_mass->GetYaxis()->SetTitle("Occurrences");
  diff_inv_mass->SetLineColor(kBlue);
  diff_inv_mass->SetFillColor(kCyan);

  same_inv_mass->GetXaxis()->SetTitle(
      "Invariant mass of particles with same charge");
  same_inv_mass->GetYaxis()->SetTitle("Occurrences");
  same_inv_mass->SetLineColor(kBlue);
  same_inv_mass->SetFillColor(kCyan);

  diff_inv_mass_ka_pi->GetXaxis()->SetTitle(
      "Invariant mass of kaons and pions with different charge");
  diff_inv_mass_ka_pi->GetYaxis()->SetTitle("Occurrences");
  diff_inv_mass_ka_pi->SetLineColor(kBlue);
  diff_inv_mass_ka_pi->SetFillColor(kCyan);

  same_inv_mass_ka_pi->GetXaxis()->SetTitle(
      "Invariant mass of kaons and pions with same charge");
  same_inv_mass_ka_pi->GetYaxis()->SetTitle("Occurrences");
  same_inv_mass_ka_pi->SetLineColor(kBlue);
  same_inv_mass_ka_pi->SetFillColor(kCyan);

  k_inv_mass->GetXaxis()->SetTitle(
      "Invariant mass of particles from a decayment");
  k_inv_mass->GetYaxis()->SetTitle("Occurrences");
  k_inv_mass->SetLineColor(kBlue);
  k_inv_mass->SetFillColor(kCyan);

  // drawing histos on a canvas
  can->Divide(4, 3);
  can->cd(1);
  histoparticles->DrawCopy("H");
  histoparticles->DrawCopy("E, P, same");
  can->cd(2);
  phi_distrib->DrawCopy("H");
  phi_distrib->DrawCopy("E, P, same");
  can->cd(3);
  theta_distrib->DrawCopy("H");
  theta_distrib->DrawCopy("E, P, same");
  can->cd(4);
  momentum_distrib->DrawCopy("H");
  momentum_distrib->DrawCopy("E, P, same");
  can->cd(5);
  trans_momentum_distrib->DrawCopy("H");
  trans_momentum_distrib->DrawCopy("E, P, same");
  can->cd(6);
  energy_distrib->DrawCopy("H");
  energy_distrib->DrawCopy("E, P, same");
  can->cd(7);
  inv_mass->DrawCopy("H");
  inv_mass->DrawCopy("E, P, same");
  can->cd(8);
  diff_inv_mass->DrawCopy("H");
  diff_inv_mass->DrawCopy("E, P, same");
  can->cd(9);
  same_inv_mass->DrawCopy("H");
  same_inv_mass->DrawCopy("E, p, same");
  can->cd(10);
  diff_inv_mass_ka_pi->DrawCopy("H");
  diff_inv_mass_ka_pi->DrawCopy("E, P, same");
  can->cd(11);
  same_inv_mass_ka_pi->DrawCopy("H");
  same_inv_mass_ka_pi->DrawCopy("E, P, same");
  can->cd(12);
  k_inv_mass->DrawCopy("H");
  k_inv_mass->DrawCopy("E, P, same");

  TH1F *difference1 = new TH1F("difference1",
                               "Difference between opposite and same charge "
                               "kaons and pions distribution",
                               400, 0, 5);
  difference1->Sumw2();
  difference1->Add(diff_inv_mass_ka_pi, same_inv_mass_ka_pi, 1, -1);
  difference1->GetXaxis()->SetTitle("Invariant mass");
  difference1->GetYaxis()->SetTitle("Occurrences");
  difference1->SetLineColor(kBlue);
  difference1->SetFillColor(kCyan);

  TH1F *difference2 = new TH1F(
      "difference2",
      "Difference between opposite and same charge particle distribution",
      400, 0, 5);
  difference2->Sumw2();
  difference2->Add(diff_inv_mass, same_inv_mass, 1, -1);
  difference2->GetXaxis()->SetTitle("Invariant mass");
  difference2->GetYaxis()->SetTitle("Occurrences");
  difference2->SetLineColor(kBlue);
  difference2->SetFillColor(kCyan);

  TF1 *fun4 = new TF1("fun4", "gaus", 0.6, 1.3);
  fun4->SetLineColor(kRed);
  fun4->SetLineWidth(2);
  TF1 *fun5 = new TF1("fun5", "gaus", 0.6, 1.3);
  fun5->SetLineColor(kRed);
  fun5->SetLineWidth(2);
  TF1 *fun6 = new TF1("fun6", "gaus", 0.6, 1.3);
  fun6->SetLineColor(kRed);
  fun6->SetLineWidth(2);

  k_inv_mass->Fit("fun4", "Q");
  std::cout
      << "Fitting invariant mass distribution of particles from decayment:\n";
  for (int i = 0; i < 3; ++i) {
    std::cout << "p[" << i << "]: " << fun4->GetParameter(i) << " +/- "
              << fun4->GetParError(i) << "\n\n";
  }
  difference1->Fit("fun5", "Q");
  std::cout << "Fitting difference between opposite and same charge kaons and "
               "pions distribution:\n";
  for (int i = 0; i < 3; ++i) {
    std::cout << "p[" << i << "]: " << fun5->GetParameter(i) << " +/- "
              << fun5->GetParError(i) << "\n\n";
  }
  difference2->Fit("fun6", "Q");
  std::cout << "Fitting difference between opposite and same charge particle "
               "distribution:\n";
  for (int i = 0; i < 3; ++i) {
    std::cout << "p[" << i << "]: " << fun6->GetParameter(i) << " +/- "
              << fun6->GetParError(i) << "\n\n";
  }

  TCanvas *diff_can = new TCanvas("diff_can", "Differences", 200, 10, 600, 400);
  diff_can->Divide(1, 3);
  diff_can->cd(1);
  k_inv_mass->DrawCopy("H");
  k_inv_mass->DrawCopy("E, P, same");
  diff_can->cd(2);
  difference1->DrawCopy("H");
  difference1->DrawCopy("E, P, same");
  diff_can->cd(3);
  difference1->DrawCopy("H");
  difference1->DrawCopy("E, P, same");
  difference2->DrawCopy("H");
  difference2->DrawCopy("E, P, same");

  f->Close();
}