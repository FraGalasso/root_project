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
  TH1F *part_gen = (TH1F *)f->Get("part_gen");
  std::cout << "part_gen: " << part_gen->GetEntries() << '\n';

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
  std::cout << "\nFitting theta distribution:\nParameters: "
            << fun2->GetParameter(0) << " +/- " << fun2->GetParError(0) << '\n';

  momentum_distrib->GetXaxis()->SetTitle("Momentum distribution");
  momentum_distrib->GetYaxis()->SetTitle("Occurrences");
  momentum_distrib->SetLineColor(kBlue);
  momentum_distrib->SetFillColor(kCyan);
  momentum_distrib->Fit("expo", "Q");
  TF1 *fun3 = momentum_distrib->GetFunction("expo");
  fun3->SetLineColor(kRed);
  fun3->SetLineWidth(2);
  // m = 1 / (-p1)
  double fit_mean = 1/-(fun3->GetParameter(1));
  // dm = (dp1 / (-p1)) * m
  double fit_mean_err = (fun3->GetParError(1)/-(fun3->GetParameter(1)))*fit_mean;
  std::cout << "Fitting momentum distribution:\n"
            << "Mean from the histogram: " << momentum_distrib->GetMean()
            << " +/- " << momentum_distrib->GetMeanError() << '\n'
            << "Mean from the fit: " << fit_mean << " +/- "
            << fit_mean_err << '\n';

  // drawing histos on a canvas
  TCanvas *can = new TCanvas("can", "Many Histos", 200, 10, 600, 400);
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
  part_gen->DrawCopy("H");
  part_gen->DrawCopy("E, P, same");

  f->Close();
}