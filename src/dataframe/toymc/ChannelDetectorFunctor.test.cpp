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

#include <ROOT/RDataFrame.hxx>
#include "ParticleGeneratorFunctor.hpp"
#include "ChannelDetectorFunctor.hpp"
#include <gtest/gtest.h>

TEST(ChannelDetectorFunctor, Basic) {
  Qn::AxisD phi_bins("phi",4,0.,2*M_PI);
  std::vector<double> efficiencies = {1.0,1.0,1.0,1.0};
  Qn::ToyMC::ChannelDetectorFunctor channel_detector(phi_bins, efficiencies);

  ROOT::RVec<double> phi;
  for (int i = 0; i < 100; ++i) {
    if (i < 50) {
      phi.emplace_back(0.);
    } else {
      phi.emplace_back(M_PI+0.1);
    }
  }
  auto weights = channel_detector(phi);
  auto phis = channel_detector.GetPhiBins()();
}

TEST(ChannelDetectorFunctor, RDataFrame) {
  ROOT::RDataFrame df(1000);
  Qn::ToyMC::ParticleGenerator phi_generator({0.},100);
  auto df_phi = df.Define("Psi",[](){return 0.;})
      .Define("Mult",[]()->int{return 1000;})
      .DefineSlot("Phis",phi_generator, {"Psi","Mult"});

  Qn::AxisD phi_bins("phi",4,0.,2*M_PI);
  std::vector<double> efficiencies = {1.0,1.0,1.0,1.0};
  Qn::ToyMC::ChannelDetectorFunctor channel_detector(phi_bins, efficiencies);

  auto df_det = df_phi.Define("MultiplicityPerChannel", channel_detector, {"Phis"})
                      .Define("PhiPerChannel", channel_detector.GetPhiBins());

  auto multiplicites = df_det.Take<ROOT::RVec<double>>("MultiplicityPerChannel");
  auto phis = df_det.Take<ROOT::RVec<double>>("PhiPerChannel");

  EXPECT_NEAR(multiplicites.GetValue()[0][0],250,10);
  EXPECT_NEAR(phis.GetValue()[0][0],M_PI*1/4., M_PI*1./20);
}
