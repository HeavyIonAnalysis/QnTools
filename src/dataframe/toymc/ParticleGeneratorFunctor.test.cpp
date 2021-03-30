// Flow Vector Correction Framework
//
// Copyright (C) 2020  Lukas Kreis Ilya Selyuzhenkov
// Contact: l.kreis@gsi.de; ilya.selyuzhenkov@gmail.com
// For a full list of contributors please see docs/Credits
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <DataContainer.hpp>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TRandom3.h>
#include <gtest/gtest.h>

#include "ParticleGeneratorFunctor.hpp"

TEST(ParticleGeneratorFunctor, OneEvent) {
  auto gen = Qn::ToyMC::ParticleGenerator({},100);
  auto histo = TH1F("h","h",10,0.,2*M_PI);
  for (int i = 0; i < 10000; ++i) {
    auto phi = gen(0, 0.,1);
    for (auto iphi :phi) histo.Fill(iphi);
    EXPECT_GT(phi.at(0), 0.);
    EXPECT_LT(phi.at(0), 2 * M_PI);
  }
  for (int i = 0; i < histo.GetNbinsX(); ++i) {
    std::cout << histo.GetBinContent(i) << std::endl;
  }
}

TEST(ParticleGeneratorFunctor, RotatedEvents) {
  auto histo = TH1F("h","h",100,-2.*M_PI,4.*M_PI);
  auto histopsi = TH1F("sh","hs",100,0.,2*M_PI);
  auto histo2 = TH1F("h2", "h2", 100, -1., 1.);
  std::mt19937_64 random_engine(1);
  TRandom3 random;
  std::uniform_real_distribution psigen(0.,2*M_PI);
  std::vector<double> vns= {0.5};
  auto PhiPdf = [vns](double phi, double psi) {
    double value = 1.;
    for (unsigned int n = 1; n < vns.size() + 1; ++n) {
      value += 2 * vns[n - 1] * std::cos(n * (phi - psi));
    }
    return value;
  };
  auto wrap_phi = [](double phi) {
    phi = fmod(phi, 2*M_PI);
    if (phi < 0.0) phi += 2*M_PI;
    return phi;
  };
  auto phigen = std::piecewise_linear_distribution<>{
      10000, 0, 2 * M_PI, [&](const double x) { return PhiPdf(x, 0); }};
  std::vector<double> phis(2000);
  for (auto &phi : phis) {
    phi = random.Gaus(M_PI,0.5);// phigen(random_engine);
  }
  std::vector<TH1F*> histos;
  for (int iev = 0; iev < 10000; ++iev) {
    auto psi = psigen(random_engine);
    std::vector<double> tphis(phis.size());
    std::copy(phis.begin(),phis.end(), tphis.begin());
    for (auto &phi : tphis) {
      phi = phi + psi;
    }
    double qx = 0.;
    for (const auto &iphi : tphis) {
//      auto phi = wrap_phi(iphi - 0.);
      histo.Fill(iphi);
      qx += std::cos(iphi);
    }
    histo2.Fill(qx / phis.size());
    histopsi.Fill(psi);
  }
  auto c1 = TCanvas("c1","C1",800,600);
//  c1.Divide(3);
//  c1.cd(1);
  histo.Draw();
//  c1.cd(2);
//  histo2.Draw();
//  c1.cd(3);
//  histopsi.Draw();
  c1.SaveAs("test.pdf");
//  for (int i = 1; i < histo.GetNbinsX()+1; ++i) {
//    std::cout << histo.GetBinContent(i) << std::endl;
//  }
}

TEST(ParticleGeneratorFunctor, RDataFrameOneParticle) {
  ROOT::RDataFrame df(100000);
  auto phi_generator = Qn::ToyMC::ParticleGenerator({0.},100);
  auto df_phi = df.Define("Psi",[](){return 0.;})
                  .Define("Mult",[]()->int{return 1;})
                  .DefineSlot("Phis",phi_generator, {"Psi","Mult"})
                  .Define("Phi", [](ROOT::RVec<double> phis){return phis[0];}, {"Phis"});
  auto phis = df_phi.Take<ROOT::RVec<double>>("Phis");
  auto phisv = phis.GetValue();
  auto mean = df_phi.Mean("Phi");
  auto stddev = df_phi.StdDev("Phi");
  auto min = df_phi.Min("Phi");
  auto max = df_phi.Max("Phi");
  const auto expected_variance = 4./12*M_PI*M_PI;
  EXPECT_NEAR(stddev.GetValue()*stddev.GetValue(), expected_variance, expected_variance*0.01);
  EXPECT_NEAR(mean.GetValue(), M_PI, M_PI*0.01);
}

TEST(ParticleGeneratorFunctor, RDataFrameManyParticles) {
  ROOT::RDataFrame df(100000);
  auto phi_generator = Qn::ToyMC::ParticleGenerator({0.},100);
  auto df_phi = df.Define("Psi",[](){return 0.;})
      .Define("Mult",[]()->int{return 1000;})
      .DefineSlot("Phis",phi_generator, {"Psi","Mult"})
                    .Define("PhiMean", [](ROOT::RVec<double> phis){
                      auto acc = std::accumulate(phis.begin(), phis.end(), 0.);
                      return acc / phis.size();}, {"Phis"});
  auto phis = df_phi.Take<ROOT::RVec<double>>("Phis");
  auto phis_v = phis.GetValue();
  auto histo = TH1F("test", "test", 10., 0., 2*M_PI);
  for (auto phi : phis_v[0]) {
    histo.Fill(phi);
  }
  for (int i = 1; i < histo.GetNbinsX()+1; ++i) {
    std::cout << histo.GetBinContent(i) << std::endl;
  }
  auto mean = df_phi.Mean("PhiMean");
  EXPECT_NEAR(mean.GetValue(), M_PI, M_PI*0.01);
}

