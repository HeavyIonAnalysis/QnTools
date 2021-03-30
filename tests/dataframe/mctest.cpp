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
#include <TH1F.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <random>
void mctest() {
  auto histo = new TH1F("h","h", 100, 0., 2*M_PI);
  std::mt19937_64 random_engine(std::random_device{}());
  TRandom3 random;
  std::vector<double> phis(10000);
  for (auto &phi : phis) {
    auto phig = std::normal_distribution<>(5.,0.5);
    phi = phig(random_engine);
//    phi = random.Gaus(5.,0.5);
  }
  for (int iev = 0; iev < 10000; ++iev) {
    auto psig = std::uniform_real_distribution(-4.,4.);
    auto psi = psig(random_engine);
//    auto psi = random.Uniform(-4.,4.);
    std::vector<double> tphis(phis.size());
    std::copy(phis.begin(),phis.end(), tphis.begin());
    for (auto &phi : tphis) {
      histo->Fill(psi);
    }
  }
//  auto c1 = TCanvas("c1","C1",800,600);
  histo->Draw();
//  c1.SaveAs("test.pdf");
}
