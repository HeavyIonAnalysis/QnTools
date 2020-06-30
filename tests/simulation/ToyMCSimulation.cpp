// Qn Tools
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

#include <gtest/gtest.h>

#include "QnDataFrame.hpp"

#include "ChannelDetector.hpp"
#include "CorrectionManager.hpp"
#include "ParticleGenerator.hpp"
#include "TrackingDetector.hpp"

inline std::vector<std::string> CreateDetectorNames(unsigned long n_detectors,
                                                    const std::string& base_name,
                                                    const std::vector<std::string>& correction_steps) {
  std::vector<std::string> detector_names;
  std::vector<std::string> names_temp;
  for (std::size_t i = 0; i < n_detectors; ++i) {
    names_temp.push_back(base_name + std::to_string(i) + "_");
  }
  for (const auto &name : names_temp) {
    for (const auto &correction_step : correction_steps) {
      detector_names.push_back(name + correction_step);
    }
  }
  return detector_names;
}

TEST(ToyMCSimulationTest, BasicSimulation) {

std::mt19937_64 engine1(1);
std::mt19937_64 engine2(2);
std::mt19937_64 engine3(3);

std::vector<TrackingDetector<std::mt19937_64>> detectors_trk;

enum variables {
  kEvent = 0,
  kPsi,
  kPsiWeight,
  kPhiNoRotation,
  kPhiNoShift,
  kPhi,
  kPt,
  kRapidity,
  kWeight,
};

Qn::AxisD pT_axis("pT", 5, 0., 1.);
Qn::AxisD rapidity_axis("rapidity", 5, -1., 1.);

Qn::CorrectionManager man;
man.SetFillOutputTree(true);
man.SetFillCalibrationQA(true);
man.SetFillValidationQA(true);

auto displaced_phi = kPhi;
auto displaced_weights = kWeight;

man.AddVariable("Event", kEvent, 1);
man.AddVariable("phi", displaced_phi, 1);
man.AddVariable("psi", kPsi, 1);
man.AddVariable("phinorotation", kPhiNoRotation, 1);
man.AddVariable("phinoshift", kPhiNoShift, 1);

//Add Tracking
const int kGranularity = 1000;
std::vector<double> modulation1(kGranularity);
std::vector<double> modulation2(kGranularity);
std::vector<double> modulation3(kGranularity);
std::vector<double> modulation4(kGranularity);
std::vector<double> baseline(kGranularity);
float impact = 0.05;
float ax = 0.05;
float ay = 0.03;
for (int i = 0; i < kGranularity; ++i) {
auto phi_step = (2. * M_PI / kGranularity) * i;
auto efficiency_bin1 = 1 / (1 + 2 * (impact + impact + ax + ay))
    * (1 + 2 * impact * std::cos(2 * phi_step) + 2 * impact * std::sin(2 * phi_step) + 2 * ax * std::cos(4 * phi_step)
        + 2 * ay * std::sin(4 * phi_step));
auto efficiency_bin2 = 1 / (1 + 2 * impact) * (1 + 2 * impact * std::cos(2 * phi_step));
auto efficiency_bin3 = 1 / (1 + 2 * (impact + impact + ax + ay))
    * (1 + 2 * impact * std::cos(2 * phi_step) + impact * std::sin(2 * phi_step) + 2 * ax * std::cos(4 * phi_step)
        + 2 * ay * std::sin(4 * phi_step));
auto efficiency_bin4 = 1 / (1 + 2 * impact) * (1 + 2 * impact * std::cos(4 * phi_step));
modulation1[i] = efficiency_bin1;
modulation2[i] = efficiency_bin2;
modulation3[i] = efficiency_bin3;
modulation4[i] = efficiency_bin4;
baseline[i] = 1.0;
}

std::vector<std::vector<double>> efficiencies{
    baseline,
    modulation1,
    modulation2,
    modulation3,
    modulation4};

for (std::size_t i = 0; i < efficiencies.size(); ++i) {
TrackingDetector<std::mt19937_64>
    det({"DetTrk" + std::to_string(i), kGranularity, 0, 2 * M_PI}, []([[maybe_unused]] const double phi) {return true;});
det.SetChannelEfficencies(efficiencies.at(i));
detectors_trk.push_back(det);
man.AddVariable("weight_trk_" + std::to_string(i), displaced_weights + i, 1);
}

TrackingDetector<std::mt19937_64>
    detpsi({"DetPsi", kGranularity, 0, 2 * M_PI}, []([[maybe_unused]] const double phi) { return true; });
detpsi.SetChannelEfficencies(efficiencies.at(0));
man.AddVariable("weight_psi", kPsiWeight, 1);
man.AddVariable("pT", kPt, 1);
man.AddVariable("rapidity", kRapidity, 1);

ParticleGenerator<std::mt19937_64, 2, 10000> gen({0., 0.05});
std::uniform_real_distribution<> psi_distribution(0, 2 * M_PI);
std::uniform_int_distribution<> n_distribution(2000, 2000);

auto outfile = TFile::Open("correctionout.root", "RECREATE");
auto tree = new TTree();
man.ConnectOutputTree(tree);
man.SetCalibrationInputFileName("correctionin.root");

man.AddEventVariable("Event");
man.AddCorrectionAxis({"Event", 1, 0, 1});
man.AddEventHisto1D({"psi", 100, 0, 2 * M_PI});

//Tracking
for (unsigned int i = 0; i < detectors_trk.size(); ++i) {
std::string name = detectors_trk.at(i).Name();
std::string weight_name = std::string("weight_trk_") + std::to_string(i);
man.AddDetector(name, Qn::DetectorType::TRACK, "phi", weight_name, {pT_axis, rapidity_axis}, {1, 2});
Qn::Recentering recentering;
recentering.SetApplyWidthEqualization(false);
man.AddCorrectionOnQnVector(name, recentering);
Qn::TwistAndRescale twist_and_rescale;
twist_and_rescale.SetApplyRescale(true);
twist_and_rescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
man.AddCorrectionOnQnVector(name, twist_and_rescale);
man.SetOutputQVectors(name, {Qn::QVector::CorrectionStep::PLAIN, Qn::QVector::CorrectionStep::RECENTERED, Qn::QVector::CorrectionStep::TWIST, Qn::QVector::CorrectionStep::RESCALED});
man.AddHisto1D(name, {"phi", 10000, 0, 2 * M_PI}, weight_name);
man.AddHisto1D(name, {"phinorotation", 10000, 0, 2 * M_PI}, weight_name);
man.AddHisto1D(name, {"phinoshift", 10000, 0, 4 * M_PI}, weight_name);
man.AddHisto1D(name, {"pT", 1000, 0., 1.}, "Ones");
man.AddHisto1D(name, {"rapidity", 1000, -1., 1.}, "Ones");
}

std::string weight_name("weight_psi");
man.AddDetector("DetPsi", Qn::DetectorType::TRACK, "psi", weight_name, {}, {1, 2});
man.SetOutputQVectors("DetPsi", {Qn::QVector::CorrectionStep::PLAIN});

man.InitializeOnNode();
auto correction_list = man.GetCorrectionList();
auto correction_qa_list = man.GetCorrectionQAList();
man.SetCurrentRunName("test");

double *values = man.GetVariableContainer();
for (int iev = 0; iev < 50; ++iev) {
man.Reset();
values[kEvent] = 0.5;
if (man.ProcessEvent()) {
auto psi = psi_distribution(engine1);
//      values[kPsi] = psi;
//Tracking
for (int n = 0; n < 100; ++n) {
auto phi = gen.GetPhi(engine1, psi);
auto phino = gen.GetPhi(engine2, 0);
values[kPhiNoShift] = phi;
phi = std::atan2(sin(phi), cos(phi)) + M_PI;
values[kPhiNoRotation] = phino;
values[kPt] = (n + 0.5) / 100.;      // from 0 to 1
values[kRapidity] = (n - 49.5) / 50.;// from -1 to 1
//        std::cout << values[kPt] << " " << values[kRapidity] << std::endl;
int ndettrk = 0;
for (std::size_t j = 0; j < detectors_trk.size(); ++j) {
detectors_trk.at(j).Detect(phi);
detectors_trk.at(j).FillDataRec(engine3, values, displaced_phi, displaced_weights + j);
ndettrk++;
}
detpsi.Detect(psi);
detpsi.FillDataRec(engine3, values, kPsi, kPsiWeight);
man.FillTrackingDetectors();
}
}
man.ProcessCorrections();
}
man.Finalize();
outfile->cd();
tree->Write("tree");
correction_list->Write("CorrectionHistograms", TObject::kSingleKey);
correction_qa_list->Write("CorrectionQAHistograms", TObject::kSingleKey);
outfile->Close();

auto begin = std::chrono::high_resolution_clock::now();// start of timing

std::string tree_file_name = "correctionout.root";
ROOT::EnableImplicitMT();
auto file = TFile::Open(tree_file_name.data());
if (!file) return;
TTreeReader reader("tree", file);

std::vector<std::string> correction_steps_t{"PLAIN", "RECENTERED", "TWIST", "RESCALED"};
std::vector<std::string> correction_steps;
for (const auto &step : correction_steps_t) {
auto det = detectors_trk.at(0);
std::string name(det.Name() + "_" + step);
TTreeReaderValue<Qn::DataContainerQVector> value_test(reader, name.data());
reader.SetLocalEntry(1);
if (value_test.GetSetupStatus() < 0) {
std::cout << "Branch " << value_test.GetBranchName() << std::endl;
std::cout << "correction_step " + step + " not available." << std::endl;
} else {
correction_steps.push_back(step);
}
reader.Restart();
}

const unsigned int n_samples = 1000;

std::cout << reader.GetEntries(true) << std::endl;
Qn::AxisD event("Event", 1, 0, 1);
ROOT::RDataFrame df("tree", tree_file_name);
auto dfs = Qn::Correlation::Resample(df, n_samples);

auto detector_names_trk = CreateDetectorNames(detectors_trk.size(), "DetTrk", correction_steps);

constexpr auto kObs = Qn::Stats::Weights::OBSERVABLE;
auto event_axes = Qn::MakeAxes(event);

std::map<std::string, ROOT::RDF::RResultPtr<Qn::Correlation::CorrelationActionBase>> correlations;
for (const auto &detector : detector_names_trk) {
std::array<std::string, 1> inputs{detector};
std::string correlation_name{"v22"};
auto c22 = Qn::MakeAverageHelper(Qn::Correlation::MakeCorrelationAction(correlation_name + detector,
                                                                        Qn::Correlation::TwoParticle::Cumulant(2),
                                                                        inputs, {kObs}, event_axes, n_samples))
    .BookMe(dfs);
correlations.emplace(detector + "_" + correlation_name, c22);
}

auto out_file = new TFile("correlationout.root", "RECREATE");
out_file->cd();
std::vector<Qn::DataContainerStats> stats;
for (auto &correlation : correlations) {
//    correlation.second.GetValue().GetDataContainer().Write(correlation.first.data());
stats.push_back(correlation.second.GetValue().GetDataContainer());
}
out_file->Close();

auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> fp_time = end - begin;
auto seconds = std::chrono::duration_cast<std::chrono::seconds>(fp_time);
auto minutes = std::chrono::duration_cast<std::chrono::minutes>(fp_time);
std::cout << minutes.count() << " minutes " << seconds.count() - minutes.count() * 60 << " seconds" << std::endl;
}