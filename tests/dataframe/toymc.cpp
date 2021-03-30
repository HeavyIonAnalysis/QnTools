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
#include <ROOT/RVec.hxx>
#include "QnTools/QnDataFrame.hpp"
#include <random>
void toymc() {
//  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame df(1e4);

  auto particle_generator = Qn::ToyMC::ParticleGenerator({0.2, 0.},100);

  auto random_dev = std::random_device();
  std::vector<std::mt19937_64>
      random_engines(ROOT::GetThreadPoolSize()+1, std::mt19937_64{random_dev()});

  auto reaction_plane = [&random_engines](unsigned int slot) {
    auto psi = std::uniform_real_distribution<>(0., 2*M_PI);
    return psi(random_engines.at(slot));
  };
  std::vector<std::mt19937_64>
      random_engines_mult(ROOT::GetThreadPoolSize()+1, std::mt19937_64{random_dev()});
  auto multiplicity = [&random_engines_mult](unsigned int slot) {
    auto mult = std::uniform_int_distribution<>(100, 1000);
    return mult(random_engines_mult.at(slot));
  };


  auto df_generated = df.DefineSlot("psi", reaction_plane)
                        .Define("psi0", "return 0.")
                        .DefineSlot("multiplicity", multiplicity)
                        .Define("int", "return 0.5")
      .DefineSlot("phis", particle_generator, {"psi", "multiplicity"})
      .DefineSlot("phis0", particle_generator, {"psi0", "multiplicity"})
                      .Define("weights", Qn::ToyMC::MakeUnitWeightsLambda(), {"phis"});

  auto ch_default = Qn::ToyMC::ChannelDetectorFunctor({"phi",8,0.,2*M_PI},{1.,1.,1.,1.,1.,1.,1.,1.});
  auto ch_shift = Qn::ToyMC::ChannelDetectorFunctor({"phi",8,0.,2*M_PI},{0.5,0.5,1.,1.,1.,1.,0.5,0.5});
  auto ch_twist = Qn::ToyMC::ChannelDetectorFunctor({"phi",8,0.,2*M_PI},{1.,1.,0.5,0.5,1.,1.,0.5,0.5});
  auto ch_shift_twist = Qn::ToyMC::ChannelDetectorFunctor({"phi",8,0.,2*M_PI},{0.7,0.7,0.5,0.5,1.,1.,0.3,0.3});

  auto tr_default = Qn::ToyMC::TrackingDetectorFunctor([](double phi){return true;});
  auto tr_shift = Qn::ToyMC::TrackingDetectorFunctor(
            [](double phi){
              return 1.0;
            }
          );

  auto filter = [](ROOT::RVec<double> &phis, ROOT::RVec<std::size_t> &filter){
    return ROOT::VecOps::Take(phis, filter);
  };

  auto df_detected =
      df_generated.Define("mult_ch_default",ch_default, {"phis"})
          .Define("mult_ch_shift",ch_shift, {"phis"})
          .Define("mult_ch_twist",ch_twist, {"phis"})
          .Define("mult_ch_shift_twist",ch_shift_twist, {"phis"})
                  .Define("phi_ch_default", ch_default.GetPhiBins())
          .Define("phi_ch_shift", ch_shift.GetPhiBins())
          .Define("phi_ch_twist", ch_twist.GetPhiBins())
          .Define("phi_ch_shift_twist", ch_shift_twist.GetPhiBins())
                  .DefineSlot("filter_tr_default", tr_default, {"phis"})
                  .DefineSlot("filter_tr_shift", tr_shift, {"phis"})
                  .Define("mult_tr_default", filter, {"weights", "filter_tr_default"})
                  .Define("mult_tr_shift", filter, {"weights", "filter_tr_shift"})
                  .Define("phi_tr_default", filter, {"phis", "filter_tr_default"})
                  .Define("phi_tr_shift", filter, {"phis", "filter_tr_shift"});

  auto normalizer = Qn::QVectorNormalizationFunctor(Qn::QVector::Normalization::M);

  auto add_to_df = [&normalizer](auto dataframe, auto q, auto phi, auto weight){
    auto ret_df = dataframe
      .Define(q.GetName(), q, {phi, weight})
      .Define(normalizer.GetName(q.GetName()), normalizer, {q.GetName()})
      .Define(normalizer.GetName(q.GetName())+"_x",[](const Qn::DataContainerQVector &q){return q.At(0).x(1);},
                              {normalizer.GetName(q.GetName())})
      .Define(normalizer.GetName(q.GetName())+"_y",[](const Qn::DataContainerQVector &q){return q.At(0).y(1);},
                              {normalizer.GetName(q.GetName())});
    return ret_df;
  };

  auto q_gen = Qn::MakeQVectorFunctor("q_gen", {1, 2});
  auto q_ch_default = Qn::MakeQVectorFunctor("q_ch_default", {1, 2});
  auto q_ch_shift = Qn::MakeQVectorFunctor("q_ch_shift", {1, 2});
  auto q_ch_twist = Qn::MakeQVectorFunctor("q_ch_twist", {1, 2});
  auto q_ch_shift_twist = Qn::MakeQVectorFunctor("q_ch_shift_twist", {1, 2});

  auto df_1 = add_to_df(df_detected, q_gen, "phis", "weights");
  auto df_2 = add_to_df(df_1, q_ch_default, "phi_ch_default", "mult_ch_default");
  auto df_3 = add_to_df(df_2, q_ch_shift, "phi_ch_shift", "mult_ch_shift");
  auto df_4 = add_to_df(df_3, q_ch_twist, "phi_ch_twist", "mult_ch_twist");
  auto df_5 = add_to_df(df_4, q_ch_shift_twist, "phi_ch_shift_twist", "mult_ch_shift_twist");

  auto names = df_4.GetColumnNames();
  for (auto name : names) {
    std::cout << name << std::endl;
  }

  auto correction_file = TFile::Open("correction_file.root", "RECREATE");
  auto integrated_axis = Qn::AxisD("int",1,0.,1.);
  auto correction_axis = Qn::MakeAxes(integrated_axis);
  auto recenter = Qn::Correction::RecenterBuilder(correction_axis);
  recenter.SetMinimumEntries(0);
  recenter.AddCorrection("q_ch_shift_norm_M", "rec", "", *q_ch_shift.GetPrototype());
  recenter.AddCorrection("q_ch_default_norm_M", "rec","", *q_ch_default.GetPrototype());
  recenter.AddCorrection("q_gen_norm_M", "rec", "", *q_gen.GetPrototype());
  recenter.AddCorrection("q_ch_twist_norm_M", "rectwistrescale", "", *q_ch_twist.GetPrototype());
  recenter.AddCorrection("q_ch_shift_twist_norm_M", "rectwistrescale", "", *q_ch_shift_twist.GetPrototype());
  recenter.EnableDiagonalization("q_ch_twist_norm_M_rectwistrescale");
  recenter.EnableRescaling("q_ch_twist_norm_M_rectwistrescale");
  recenter.EnableDiagonalization("q_ch_shift_twist_norm_M_rectwistrescale");
  recenter.EnableRescaling("q_ch_shift_twist_norm_M_rectwistrescale");
  auto df_after = recenter.Apply(df_5, correction_file, nullptr);
  recenter.Write(correction_file);


  auto df_final = df_after.Define("q_ch_shift_norm_M_rec_x", [](const Qn::DataContainerQVector &q){return q.At(0).x(1);}, {"q_ch_shift_norm_M_rec"})
                          .Define("q_ch_shift_norm_M_rec_y", [](const Qn::DataContainerQVector &q){return q.At(0).y(1);}, {"q_ch_shift_norm_M_rec"})
      .Define("q_ch_twist_norm_M_rectwistrescale_x", [](const Qn::DataContainerQVector &q){return q.At(0).x(1);}, {"q_ch_twist_norm_M_rectwistrescale"})
      .Define("q_ch_twist_norm_M_rectwistrescale_y", [](const Qn::DataContainerQVector &q){return q.At(0).y(1);}, {"q_ch_twist_norm_M_rectwistrescale"})
      .Define("q_ch_shift_twist_norm_M_rectwistrescale_x", [](const Qn::DataContainerQVector &q){return q.At(0).x(1);}, {"q_ch_shift_twist_norm_M_rectwistrescale"})
      .Define("q_ch_shift_twist_norm_M_rectwistrescale_y", [](const Qn::DataContainerQVector &q){return q.At(0).y(1);}, {"q_ch_shift_twist_norm_M_rectwistrescale"});
  auto histo1 = df_final.Histo2D({"gen","h",100,-1.,1.,100, -1., 1.},"q_gen_norm_M_x","q_gen_norm_M_y");
  auto histo2 = df_final.Histo2D({"default","h",100,-1.,1.,100, -1., 1.},"q_ch_default_norm_M_x","q_ch_default_norm_M_y");
  auto histo3 = df_final.Histo2D({"shift","h",100,-1.,1.,100, -1., 1.},"q_ch_shift_norm_M_x","q_ch_shift_norm_M_y");
  auto histo4 = df_final.Histo2D({"shift_rec","h",100,-1.,1.,100, -1., 1.},"q_ch_shift_norm_M_rec_x","q_ch_shift_norm_M_rec_y");
  auto histo5 = df_final.Histo2D({"twist","h",100,-1.,1.,100, -1., 1.},"q_ch_twist_norm_M_x","q_ch_twist_norm_M_y");
  auto histo6 = df_final.Histo2D({"twist_rec_twist","h",100,-1.,1.,100, -1., 1.},"q_ch_twist_norm_M_rectwistrescale_x","q_ch_twist_norm_M_rectwistrescale_y");
  auto histo7 = df_final.Histo2D({"twist_rec_shift_twist","h",100,-1.,1.,100, -1., 1.},"q_ch_shift_twist_norm_M_rectwistrescale_x","q_ch_shift_twist_norm_M_rectwistrescale_y");
  auto histo8 = df_final.Histo2D({"before_rec_shift_twist","h",100,-1.,1.,100, -1., 1.},"q_ch_shift_twist_norm_M_x","q_ch_shift_twist_norm_M_y");

  auto c1 = new TCanvas("x","g", 800 ,600);
  histo1->DrawClone("colz");
  auto c2 = new TCanvas("c2","c2",800, 600);
  histo2->DrawClone("colz");
  auto c3 = new TCanvas("c3","c2",800, 600);
  histo3->DrawClone("colz");
  auto c4 = new TCanvas("c4","c2",800, 600);
  histo4->DrawClone("colz");
  auto c5 = new TCanvas("c5","c2",800, 600);
  histo5->DrawClone("colz");
  auto c6 = new TCanvas("c6","c2",800, 600);
  histo6->DrawClone("colz");
  auto c7 = new TCanvas("c7","c2",800, 600);
  histo7->DrawClone("colz");
  auto c8 = new TCanvas("c8","c2",800, 600);
  histo8->DrawClone("colz");
//  auto histopsi = df_generated.Histo1D({"psi","psi", 100, 0., 2*M_PI}, "psi");
//  histopsi->DrawClone("");
//

}
