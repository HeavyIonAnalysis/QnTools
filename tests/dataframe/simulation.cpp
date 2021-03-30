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
#include <random>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "QnTools/QnDataFrame.hpp"
//#include "QnDataFrame.hpp"
#include <fmt/core.h>


void simulation() {
  using namespace std::string_literals;
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame df(1e4);

  // randomize the reaction plane
  auto random_dev = std::random_device();
  std::vector<std::mt19937_64> random_engines(ROOT::GetThreadPoolSize()+1, std::mt19937_64{random_dev()});
  auto reaction_plane = [&random_engines](unsigned int slot) {
    auto psi = std::uniform_real_distribution<>(0., 2*M_PI);
    return psi(random_engines.at(slot));
  };

  // generate all particles
  auto particle_generator = Qn::ToyMC::ParticleGenerator({0.2, 0.},100);
  auto df_generated =
      df.Define("int", "return 0.5")
        .DefineSlot("psi", reaction_plane)
        .Define("psi_0", "return 0.")
        .Define("multiplicity", "return 2500")
        .DefineSlot("phis", particle_generator, {"psi", "multiplicity"})
        .Define("phisa", [](const ROOT::RVec<double> &phi){ROOT::RVec<double>phia(phi.size()/2);
          std::copy(phi.begin(), phi.begin()+phi.size()/2, phia.begin());
          return phia;}, {"phis"})
        .Define("phisb", [](const ROOT::RVec<double> &phi){ROOT::RVec<double>phib(phi.size()/2);
          std::copy(phi.begin()+phi.size()/2, phi.end(), phib.begin());
          return phib;}, {"phis"})
        .DefineSlot("phis_0", particle_generator, {"psi_0", "multiplicity"})
        .Define("weights", Qn::ToyMC::MakeUnitWeightsLambda(), {"phis"})
        .Alias("phi_gen", "phis")
        .Alias("mult_gen", "weights");

  std::vector<double> norm_channels(32, 1.);
  std::vector<double> shifted_channels(norm_channels.size(), 1.);
  std::vector<double> twisted_channels(norm_channels.size(), 1.);
  std::vector<double> squashed_channels(norm_channels.size(), 1.);
  float bad_channels = 0.3;
  for (auto i= 0u; i< norm_channels.size(); ++i) {
    norm_channels[i] = 1.;
    shifted_channels[i] =  (i < (norm_channels.size()/2)) ? 0.5 : 1.0;
    int i_quarter = i/(norm_channels.size() / 4);
    if (i_quarter == 1 || i_quarter == 3) {
      twisted_channels[i] = bad_channels;
    }
    int i_eigth = i/(norm_channels.size() /8);
    if (i_eigth == 0 || i_eigth ==3 || i_eigth ==4 || i_eigth == 7) {
      squashed_channels[i] = bad_channels;
    }
  }

  // channelized detectors
  int n_ch = norm_channels.size();
  auto ch_default = Qn::ToyMC::ChannelDetectorFunctor({"phi",n_ch,0.,2*M_PI},norm_channels);
  auto ch_shift = Qn::ToyMC::ChannelDetectorFunctor({"phi",n_ch,0.,2*M_PI},shifted_channels);
  auto ch_twist = Qn::ToyMC::ChannelDetectorFunctor({"phi",n_ch,0.,2*M_PI},twisted_channels);
  auto ch_rescale = Qn::ToyMC::ChannelDetectorFunctor({"phi",n_ch,0.,2*M_PI},squashed_channels);
  auto all_effects_channels = std::vector<double>(norm_channels);
  for (auto i = 0u; i < n_ch; ++i) {
    all_effects_channels[i] = shifted_channels[i] * twisted_channels[i] * squashed_channels[i];
  }
  auto ch_shift_twist = Qn::ToyMC::ChannelDetectorFunctor({"phi",n_ch,0.,2*M_PI}, all_effects_channels);

  // tracking detectors
  auto shift = [](double phi){return (1.0 + 0.2 *cos(phi)) / 1.2;};
  auto tr_default = Qn::ToyMC::TrackingDetectorFunctor([](double phi){return true;});
  auto tr_shift = Qn::ToyMC::TrackingDetectorFunctor(shift);
  auto tr_shift_twist = Qn::ToyMC::TrackingDetectorFunctor([](double phi){return 1.0 + 0.2* cos(phi) + 0.5 * cos(2*phi+0.4) / 1.7;});

  auto filter = [](ROOT::RVec<double> &phis, ROOT::RVec<std::size_t> &filter){return ROOT::VecOps::Take(phis, filter);};

  auto add_tracking_detector = [filter](auto &df, auto name, auto tracking_detector){
    df = df.DefineSlot(fmt::format("filter_{}",name), tracking_detector, {"phis"})
      .Define(fmt::format("mult_{}",name), filter, {"weights", fmt::format("filter_{}",name)})
      .Define(fmt::format("phi_{}",name), filter, {"phis", fmt::format("filter_{}",name)});
  };

  auto add_channel_detector = [](auto &df, auto name, auto channel_detector){
    df = df.Define(fmt::format("mult_{}",name), channel_detector, {"phis"})
           .Define(fmt::format("phi_{}",name), channel_detector.GetPhiBins())
        .Alias(fmt::format("phi_{}_diag",name), fmt::format("phi_{}",name))
        .Alias(fmt::format("mult_{}_diag",name), fmt::format("mult_{}",name))
        .Alias(fmt::format("phi_{}_rescale",name), fmt::format("phi_{}",name))
        .Alias(fmt::format("mult_{}_rescale",name), fmt::format("mult_{}",name))
        .Alias(fmt::format("phi_{}_recenter",name), fmt::format("phi_{}",name))
        .Alias(fmt::format("mult_{}_recenter",name), fmt::format("mult_{}",name));
  };

  //add detectors to dataframe
  auto df_detected = df_generated;
  add_tracking_detector(df_detected, "tr_default", tr_default);
//  add_tracking_detector(df_detected, "tr_shift", tr_shift);
//  add_tracking_detector(df_detected, "tr_shift_twist", tr_shift_twist);
  add_channel_detector(df_detected, "ch_default", ch_default);
  add_channel_detector(df_detected, "ch_shift", ch_shift);
  add_channel_detector(df_detected, "ch_shift_twist_squash", ch_shift_twist);
  add_channel_detector(df_detected, "ch_twist", ch_twist);
  add_channel_detector(df_detected, "ch_squash", ch_rescale);


  //create flow vectors and add corrections
  auto normalizer = Qn::QVectorNormalizationFunctor(Qn::QVector::Normalization::M);
  auto add_qvector = [&normalizer](auto &df, auto name){
    auto q = Qn::MakeQVectorFunctor(fmt::format("q_{}",name), {1, 2});
    auto name_q = q.GetName();
    auto name_norm = normalizer.GetName(name_q);
    df = df.Define(q.GetName(), q, {fmt::format("phi_{}", name), fmt::format("mult_{}", name)})
           .Define(name_norm, normalizer, {name_q});
    return q;
  };
  std::vector<std::string> detectors = {
      "gen",
      "ch_default","ch_shift",
      "ch_shift_twist_squash_recenter",
      "ch_shift_twist_squash_diag",
      "ch_shift_twist_squash_rescale",
      "ch_twist_diag", "ch_squash_rescale",
//      "ch_default_a","ch_shift_a","ch_shift_twist_rescale_a",
//      "ch_twist_a", "ch_rescale_a",
//      "ch_default_b","ch_shift_b","ch_shift_twist_rescale_b",
//      "ch_twist_b", "ch_rescale_b",
      "tr_default"
  };
  auto df_qvector = df_detected;
  auto correction_file = TFile::Open("correction_file.root", "RECREATE");
  auto axis = Qn::MakeAxes(Qn::AxisD{"int", 1, 0., 1.});
  auto recenter = Qn::Correction::RecenterBuilder(axis);
  recenter.SetMinimumEntries(0);
  std::vector<std::string> qvectors;
  Qn::DataContainerQVector q_proto;
  for (const auto &detector : detectors) {
    auto q = add_qvector(df_qvector, detector);
    auto name = normalizer.GetName(q.GetName());
    recenter.AddCorrection(name, "corrected", "", *q.GetPrototype());
    if (name.find("diag") != std::string::npos) {
      recenter.EnableDiagonalization(name);
//      recenter.EnableRescaling(name);
    }
    if (name.find("rescale")!= std::string::npos) {
      recenter.EnableDiagonalization(name);
      recenter.EnableRescaling(name);
    }
    qvectors.emplace_back(name);
    qvectors.emplace_back(name+"_corrected");
    q_proto = *q.GetPrototype();
  }
  //apply corrections
  auto df_corrected = recenter.Apply(df_qvector, correction_file, nullptr);
  recenter.Write(correction_file);

  auto add_xy = [](auto &df, std::string name) {
    df = df.Define(fmt::format("{}_x", name), [](const Qn::DataContainerQVector &q){return q.At(0).x(1);}, {name})
           .Define(fmt::format("{}_y", name), [](const Qn::DataContainerQVector &q){return q.At(0).y(1);}, {name})
           .Define(fmt::format("{}_x2", name), [](const Qn::DataContainerQVector &q){return q.At(0).x(2);}, {name})
           .Define(fmt::format("{}_y2", name), [](const Qn::DataContainerQVector &q){return q.At(0).y(2);}, {name});
  };
  for (const auto &q : qvectors) { add_xy(df_corrected, q); }
  auto df_histo = df_corrected;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histos1d;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> histos2d;
  for (const auto &q : qvectors) {
    histos1d.emplace_back(
      df_histo.Histo1D({fmt::format("h_{}_x", q).data(), ";x;N",100,-1.,1.}, fmt::format("{}_x", q))
    );
    histos1d.emplace_back(
      df_histo.Histo1D({fmt::format("h_{}_y", q).data(), ";y;N",100,-1.,1.}, fmt::format("{}_y", q))
    );
    histos2d.emplace_back(
      df_histo.Histo2D({fmt::format("h_{}_x_y", q).data(), ";x;y;N",100,-1.,1.,100,-1.,1.}, fmt::format("{}_x", q), fmt::format("{}_y", q))
    );
    histos2d.emplace_back(
        df_histo.Histo2D({fmt::format("h_{}_x_y_2", q).data(), ";x2;y2;N",100,-1.,1.,100,-1.,1.}, fmt::format("{}_x2", q), fmt::format("{}_y2", q))
    );
  }

  int n = 10;
  auto df_samples = Qn::Correlation::Resample(df_histo, n);
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wn = Qn::Correlation::UseWeights::No;
  auto correlations = Qn::Correlation::CorrelationBuilder{&df_samples, n, axis};

  auto xx = [](const Qn::QVector &q1, const Qn::QVector &q2) {
    return q1.x(1) * q2.x(1);
  };

  auto c2 = [](const Qn::QVector &q) {
    auto ret = std::numeric_limits<float>::quiet_NaN();
    const auto S = q.sumweights();
    if (S >= 2.) {
      const auto Q = q.DeNormal();
      ret = (Qn::ScalarProduct(Q, Q, 1) - S) / (S * S - S);
    }
    return ret;
  };

  auto n2 = [](const Qn::QVector &q) {
    auto ret = 0.;
    const auto S = q.sumweights();
    if (q.n() >= 2) {
      ret = S * S - S;
    }
    return ret;
  };

  auto w1 = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  namespace P2 = Qn::Correlation::TwoParticle;
  std::array<std::string, 1> c2q{"q_gen_norm_M"};

  auto add_correlation = [& correlations, w1, wn, q_proto](const std::string &channel_detector) {
    std::array<std::string, 2> tr_ch{"q_tr_default_norm_M", channel_detector};
    std::array<std::string, 2> ch_gen{channel_detector, channel_detector};
    std::array<std::string, 2> tr_gen{"q_tr_default_norm_M", channel_detector};
    auto name = fmt::format("{}",channel_detector);
    correlations.AddCorrelationWithInitializationObject(name+"_tr_ch", P2::ScalarProduct(1,1), w1, wn, tr_ch, {q_proto,q_proto});
    correlations.AddCorrelationWithInitializationObject(name+"_ch_gen", P2::ScalarProduct(1,1), w1, wn, ch_gen, {q_proto,q_proto});
    correlations.AddCorrelationWithInitializationObject(name+"_tr_gen", P2::ScalarProduct(1,1), w1, wn, tr_gen, {q_proto,q_proto});
  };

//  auto add_correlation = [& correlations, w1, wn, q_proto](const std::string &channel_detector) {
//    auto first_part = channel_detector.substr(0, channel_detector.find("_norm_M"));
//    auto second_part = channel_detector.substr(channel_detector.find("_norm_M"), std::string::npos);
//    auto name_a = first_part + "_a" + second_part;
//    auto name_b = first_part + "_b" + second_part;
//    std::cout << name_a << " " << name_b << std::endl;
//    std::array<std::string, 2> tr_ch{"q_tr_default_norm_M", channel_detector};
//    std::array<std::string, 2> res{name_a, name_b};
//    auto name = fmt::format("{}",channel_detector);
//    correlations.AddCorrelationWithInitializationObject(name+"xx", P2::xx(1,1), w1, wn, tr_ch, {q_proto,q_proto});
//    correlations.AddCorrelationWithInitializationObject(name+"yy", P2::yy(1,1), w1, wn, tr_ch, {q_proto,q_proto});
//    correlations.AddCorrelationWithInitializationObject(name+"res_xx", P2::xx(1,1), w1, wn, res, {q_proto,q_proto});
//    correlations.AddCorrelationWithInitializationObject(name+"res_yy", P2::yy(1,1), w1, wn, res, {q_proto,q_proto});
//  };

  std::vector<std::string> correlations_v = {
    "q_ch_shift_twist_squash_recenter_norm_M",
    "q_ch_shift_twist_squash_recenter_norm_M_corrected",
    "q_ch_shift_twist_squash_diag_norm_M",
    "q_ch_shift_twist_squash_diag_norm_M_corrected",
    "q_ch_shift_twist_squash_rescale_norm_M",
    "q_ch_shift_twist_squash_rescale_norm_M_corrected",
    "q_ch_shift_norm_M",
    "q_ch_shift_norm_M_corrected",
    "q_ch_twist_diag_norm_M",
    "q_ch_twist_diag_norm_M_corrected",
    "q_ch_squash_rescale_norm_M",
    "q_ch_squash_rescale_norm_M_corrected"
  };

  for (auto &name : correlations_v) { add_correlation(name); }

  correlations.AddCorrelationWithInitializationObject("c2", c2, n2, wn, c2q, {q_proto});

  auto correlation_results = correlations.GetResults();

  auto output = TFile::Open("histos.root", "RECREATE");


//  for (int i =0; i < correlation_results.size()-1; i=i+4) {
//    std::cout << correlation_results[i]->GetName() << " "
//              << correlation_results[i+1]->GetName() << " "
//        << correlation_results[i+2]->GetName() << " "
//    << correlation_results[i+3]->GetName() << std::endl;
//    auto txx = correlation_results[i]->GetDataContainer();
//    auto tyy = correlation_results[i+1]->GetDataContainer();
//    auto tresxx = correlation_results[i+2]->GetDataContainer();
//    auto tresyy = correlation_results[i+3]->GetDataContainer();
//    auto xx = Qn::DataContainerStatCalculate{txx};
//    auto yy = Qn::DataContainerStatCalculate{tyy};
//    auto resxx = Qn::DataContainerStatCalculate{tresxx};
//    auto resyy = Qn::DataContainerStatCalculate{tresyy};
//    auto v1xx = xx / Qn::Sqrt((resxx)) * sqrt(2);
//    auto v1yy = yy / Qn::Sqrt((resyy)) * sqrt(2);
//    auto namexx = fmt::format("vxx_{}", correlations_v[(i)/4]);
//    auto nameyy = fmt::format("vyy_{}", correlations_v[(i)/4]);
//    v1xx.Write(namexx.data());
//    v1yy.Write(nameyy.data());
//    xx.Write(fmt::format("xx_{}", correlations_v[i/4]).data());
//    resxx.Write(fmt::format("resxx_{}", correlations_v[i/4]).data());
//  }

  for (int i =0; i < correlation_results.size()-1; i=i+3) {
    std::cout << correlation_results[i]->GetName() << " " << correlation_results[i+1]->GetName() << " "  << correlation_results[i+2]->GetName() << std::endl;
    auto txxab = correlation_results[i]->GetDataContainer();
    auto txxbc = correlation_results[i+1]->GetDataContainer();
    auto txxac = correlation_results[i+2]->GetDataContainer();
    auto xxab = Qn::DataContainerStatCalculate{txxab};
    auto xxac = Qn::DataContainerStatCalculate{txxac};
    auto xxbc = Qn::DataContainerStatCalculate{txxbc};
    auto v1xx = xxab / Qn::Sqrt((xxab*xxbc) / xxac);
    auto name = fmt::format("v_{}", correlations_v[(i)/3]);
    v1xx.Write(name.data());
    auto graph = Qn::ToTGraph(v1xx);
    graph->Write(fmt::format("g_{}",name).data());
  }

  for (auto histo : histos1d) histo->Write();
  for (auto histo : histos2d) histo->Write();

  output->Close();
}
