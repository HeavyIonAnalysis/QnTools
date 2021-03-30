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
#include <TChain.h>
#include <gtest/gtest.h>

#include "AverageHelper.hpp"
#include "AxesConfiguration.hpp"
#include "CorrelationAction.hpp"
#include "CorrelationBuilder.hpp"
#include "ReSampleFunctor.hpp"
#include "CorrelationFunctions.hpp"
#include "QnDataFrame.hpp"

TEST(CorrelationBuilderUnitTest, RDataFrame) {
  ROOT::RDataFrame df(100);
  std::size_t phi_size = 4;
  Qn::DataContainerQVector qvec_proto;
  auto event_axis = Qn::AxisD("event", 2, 0, 2);
  for (auto &q : qvec_proto) {
    q.ActivateHarmonic(1);
    q.InitializeHarmonics();
  }
  auto df0 = df.DefineSlotEntry("event",
                                [](unsigned int, ULong64_t entry) {
                                  double a = 0;
                                  if (entry > 90) a = 1.;
                                  return a;
                                },
                                {});
  auto df1 = df0.Define("q",
                        [qvec_proto](double event) -> Qn::DataContainerQVector {
                          auto ret = qvec_proto;
                          for (int i = 0; i < 100; ++i) {
                            ret.At(0).Add(0, 1.);
                          }
                          ret.At(0).CheckQuality();
                          ret.At(0).Normal(Qn::QVector::Normalization::M);
                          return ret;
                        },
                        {"event"});
  auto axes = Qn::MakeAxes(event_axis);
  auto c1 = [](const Qn::QVector &a) { return a.x(1); };
  std::vector<Qn::DataContainerQVector> initializations{qvec_proto};
  int n = 5;
  auto df2 = Qn::Correlation::Resample(df1, n);
  auto y = Qn::Correlation::UseWeights::Yes;
  std::array<std::string, 1> qs = {"q"};
  auto correlationvector = Qn::Correlation::CorrelationBuilder(&df2, n, axes);
  correlationvector.AddCorrelationWithInitializationObject("test", c1, c1, y,
                                                           qs, initializations);

  auto results = correlationvector.GetResults();
  auto x = results[0]->GetDataContainer();
  auto xx = Qn::DataContainerStatCalculate(x);
}

TEST(CorrelationBuilderUnitTest, FullChain) {
  const int n_slots = 40;
  ROOT::EnableImplicitMT(n_slots);
  std::string line =
      "/Users/lukas/phd/HeavyIonAnalysis/kreisl/QnTools/runtree.root";
  bool firstfile = true;
  TChain chain("tree");
  std::map<long long, double> runid_runnumber;
  double irunid = 0;
  std::vector<double> vx_mean;
  std::vector<double> vx_stddev;
  std::vector<double> vy_mean;
  std::vector<double> vy_stddev;
  std::vector<double> vz_mean;
  std::vector<double> vz_stddev;

  auto run = 137161;
  runid_runnumber.emplace(run, irunid);
  ROOT::RDataFrame df_vtx("tree", line, {"RunNumber", "VtxX", "VtxY", "VtxZ"});
  auto vxm = df_vtx.Mean("VtxX");
  auto vxs = df_vtx.StdDev("VtxX");
  auto vym = df_vtx.Mean("VtxY");
  auto vys = df_vtx.StdDev("VtxY");
  auto vzm = df_vtx.Mean("VtxZ");
  auto vzs = df_vtx.StdDev("VtxZ");
  vx_mean.push_back(vxm.GetValue());
  vx_stddev.push_back(vxs.GetValue());
  vy_mean.push_back(vym.GetValue());
  vy_stddev.push_back(vys.GetValue());
  vz_mean.push_back(vzm.GetValue());
  vz_stddev.push_back(vzs.GetValue());
  chain.Add(line.data());
  TTreeReader r(&chain);
  ROOT::RDataFrame df(chain);
  auto dfm = df.Define("Multiplicity", [](const Qn::DataContainerQVector &a){return a.At(0).n();}, {"TPC_PLAIN"})
                 .Define("ITSLayer0", "NClsITSLayers[0]")
                 .Define("ITSLayer1", "NClsITSLayers[1]")
                 .Define("ITSLayer2", "NClsITSLayers[2]")
                 .Define("ITSLayer3", "NClsITSLayers[3]")
                 .Define("ITSLayer4", "NClsITSLayers[4]")
                 .Define("ITSLayer5", "NClsITSLayers[5]")
                 .Define("ITSLayerSum",
                         "NClsITSLayers[0]+NClsITSLayers[1]+NClsITSLayers[2]+"
                         "NClsITSLayers[3]+NClsITSLayers[4]+NClsITSLayers[5]");

  std::vector<ROOT::RDF::RResultPtr<TProfile>> pprofiles;
  for (int i = 0; i < 6; ++i) {
    auto name = TString::Format("nTPCTracksHybrid_NClsITSLayer_%d_cut", i);
    auto name_axes = TString::Format(";n cls ITS Layer {};N Tracks Hybrid");
    auto its_layer_variable = TString::Format("ITSLayer%d", i);
    pprofiles.push_back(
        dfm.Profile1D({name.Data(), name_axes.Data(), 250, 0., 8000., "s"},
                      its_layer_variable.Data(), "NTPCTracksHybrid"));
  }
  std::vector<TProfile> profiles;
  for (auto pr : pprofiles) profiles.push_back(pr.GetValue());

  auto get_vertex_mean_sigma = [](const std::vector<double> &mean,
                                  const std::vector<double> &sigma) {
    return [&mean, sigma](double id, double x) {
      auto xm = mean[(int)id];
      auto xs = sigma[(int)id];
      return (x - xm) / xs;
    };
  };

  auto dfid =
      df.Define("RunID",
                [runid_runnumber](long long run) {
                  try {
                    return runid_runnumber.at(run);
                    /// catch if invalid run number was written
                  } catch (std::out_of_range &e) {
                    return -1.;
                  }
                },
                {"RunNumber"})
          .Define("VtxXnSigma", get_vertex_mean_sigma(vx_mean, vx_stddev),
                  {"RunID", "VtxX"})
          .Define("VtxYnSigma", get_vertex_mean_sigma(vy_mean, vy_stddev),
                  {"RunID", "VtxY"})
          .Define("VtxZnSigma", get_vertex_mean_sigma(vz_mean, vz_stddev),
                  {"RunID", "VtxZ"})
          .Filter("CentralityV0M < 80. && CentralityV0M > 0.")
          .Filter("VtxZ > -10 && VtxZ < 10")
          .Filter(
              [&profiles](ROOT::RVec<double> nits, double ntrackshybrid) {
                bool flag = true;
                int i = 0;
                for (const auto &prof : profiles) {
                  auto axis = prof.GetXaxis();
                  auto bin = axis->FindBin(nits[i]);
                  auto mid = prof.GetBinContent(bin);
                  auto err = prof.GetBinError(bin);
                  if ((ntrackshybrid > (mid + 3. * err)) ||
                      (ntrackshybrid < (mid - 3. * err)))
                    flag = false;
                  ++i;
                }
                return flag;
              },
              {"NClsITSLayers", "NTPCTracksHybrid"})
          .Filter("RunID >= 0.");

  Qn::AxisD run_axis ("RunID", (int) irunid,  0., irunid);
  Qn::AxisD cent_axis ("CentralityV0M", 80,   0., 80.);
  Qn::AxisD vtxx_axis("VtxXnSigma", {-5., -0.43072, 0.43072, 5.});
  Qn::AxisD vtxy_axis("VtxYnSigma", {-5., -0.43072, 0.43072, 5.});
  Qn::AxisD vtxz_axis("VtxZnSigma", {-5., -0.43072, 0.43072, 5.});
  auto axes_correlation = Qn::MakeAxes(cent_axis);

  auto correction_file = new TFile("correction.root","UPDATE");
  auto recenter = Qn::Correction::CorrectionBuilder(axes_correlation);
  recenter.SetMinimumEntries(0);
  recenter.AddCorrection("ZNA_PLAIN", "4d", "");
  recenter.AddCorrection("ZNC_PLAIN", "4d", "");
  auto df_after = recenter.Apply(df, correction_file, &r);


  // ---------------- //
  // correlation step //
  // ---------------- //
  int n = 50;
  auto df_samples = Qn::Correlation::Resample(df_after, n);

  namespace MH = Qn::Correlation::MixedHarmonics;
  namespace P2 = Qn::Correlation::TwoParticle;
  auto wn = Qn::Correlation::UseWeights::No;
  auto w1 = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  std::string cs = "_4d";
  std::string cs_tpc = "";
  std::string tpc = "TPCn_PLAIN";
  std::string zna = "ZNA_PLAIN";
  std::string znc = "ZNC_PLAIN";

  auto correlation_builder =
      Qn::Correlation::CorrelationBuilder{&df_samples, n, axes_correlation};
  std::array<std::string, 2> zn{zna, znc};
  std::array<std::string, 2> zncorrected{zna+cs, znc+cs};
  correlation_builder.AddCorrelationWithInternalReader("XX", P2::xx(1, 1), w1,
                                                       wn, zncorrected, zn);
  correlation_builder.AddCorrelationWithInternalReader("YY", P2::yy(1, 1), w1,
                                                       wn, zncorrected, zn);
  correlation_builder.AddCorrelationWithInternalReader("XY", P2::xy(1, 1), w1,
                                                       wn, zncorrected, zn);
  correlation_builder.AddCorrelationWithInternalReader("YX", P2::yx(1, 1), w1,
                                                       wn, zncorrected, zn);

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto correlation_file = TFile::Open("test_correlation.root", "RECREATE");
  correlation_file->cd();
  auto results = correlation_builder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  correlation_file->Close();
}
