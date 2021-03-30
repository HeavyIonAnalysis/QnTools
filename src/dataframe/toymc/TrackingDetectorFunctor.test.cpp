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
#include "TrackingDetectorFunctor.hpp"
#include <gtest/gtest.h>

TEST(TrackingDetectorFunctor, Basic) {
  auto efficiency_function = std::function<double(double)>([](double phi){
    double is_detected = 0.;
    if (phi < M_PI) is_detected = 1.;
    return is_detected;
  });
  auto tracking_detector = Qn::ToyMC::TrackingDetectorFunctor(efficiency_function);
  ROOT::RVec<double> phi;
  for (int i = 0; i < 100; ++i) {
    if (i < 50) {
      phi.emplace_back(0.);
    } else {
      phi.emplace_back(M_PI+0.1);
    }
  }
  auto detected = tracking_detector(0, phi);
  EXPECT_EQ(detected.size(), 50);
}

TEST(TrackingDetectorFunctor, OneEvent) {
  auto efficiency_function = std::function<double(
      double)>([](double phi) {
    bool is_detected = 0.;
    if (phi < M_PI) is_detected = 1.;
    return is_detected;
  });
  auto tracking_detector = Qn::ToyMC::TrackingDetectorFunctor(efficiency_function);
  ROOT::RVec<double> phi;
  for (int i = 0; i < 100; ++i) {
    if (i < 50) {
      phi.emplace_back(0.);
    } else {
      phi.emplace_back(M_PI+0.1);
    }
  }
  auto detected = tracking_detector(0, phi);
  auto weights = Qn::ToyMC::MakeUnitWeights(phi);
  auto selected_weights = ROOT::VecOps::Take(weights, detected);
  auto selected_phi = ROOT::VecOps::Take(phi, detected);
  EXPECT_EQ(selected_weights.size(), 50);
  EXPECT_EQ(selected_phi.size(), 50);
  EXPECT_EQ(detected.size(), 50);
}

TEST(TrackingDetectorFunctor, RDataFrame) {
  ROOT::RDataFrame df(1000);
  Qn::ToyMC::ParticleGenerator phi_generator({0.},100);
  auto df_phi = df.Define("Psi",[](){return 0.;})
      .Define("Mult",[]()->int{return 1000;})
      .DefineSlot("Phis",phi_generator, {"Psi","Mult"});

  Qn::AxisD phi_bins("phi",4,0.,2*M_PI);
  std::vector<double> efficiencies = {1.0,1.0,1.0,1.0};
  auto efficiency_function = std::function<double(double)>([](double phi) {
    bool is_detected = 0.;
    if (phi < M_PI) is_detected = 1.;
    return is_detected;
  });
  Qn::ToyMC::TrackingDetectorFunctor tracking_detector(efficiency_function);

  auto filter = [](ROOT::RVec<double> &phis, ROOT::RVec<std::size_t> &filter){
    return ROOT::VecOps::Take(phis, filter);
  };
  auto df_det = df_phi.DefineSlot("DetectionFilter", tracking_detector, {"Phis"})
                      .Define("Phis_Detected",filter, {"Phis","DetectionFilter"})
                      .Define("Weights", Qn::ToyMC::MakeUnitWeightsLambda(), {"Phis"})
                      .Define("Weights_Detected", filter, {"Weights","DetectionFilter"});

  auto weights = df_det.Take<ROOT::RVec<double>>("Weights_Detected");
  auto phis = df_det.Take<ROOT::RVec<double>>("Phis_Detected");

  EXPECT_NEAR(weights.GetValue()[0].size(),500,50);
  EXPECT_NEAR(phis.GetValue()[0].size(),500, 50);
}
