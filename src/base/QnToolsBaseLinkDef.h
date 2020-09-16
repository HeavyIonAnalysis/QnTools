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

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace Qn;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class Qn::Axis < double> + ;
#pragma link C++ class Qn::QVec + ;
#pragma link C++ class Qn::QVector + ;
#pragma link C++ class Qn::Statistics + ;
#pragma link C++ class Qn::BootStrap + ;
#pragma link C++ class Qn::Stat + ;
#pragma link C++ class Qn::StatCollect + ;
#pragma link C++ class Qn::StatCalculate + ;
#pragma link C++ class Qn::DataContainer < Qn::Statistics, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainer < Qn::StatCollect, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainer < Qn::StatCalculate, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainer < Qn::QVector, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainer < double, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainer < TH1F, Qn::Axis < double>> + ;
#pragma link C++ class Qn::DataContainerHelper + ;

#pragma link C++ typedef Qn::AxisF;
#pragma link C++ typedef Qn::AxisD;
#pragma link C++ typedef Qn::BinnedStatistics;
#pragma link C++ typedef Qn::DataContainerStatistics;
#pragma link C++ typedef Qn::DataContainerStatCollect;
#pragma link C++ typedef Qn::DataContainerStatCalculate;
#pragma link C++ typedef Qn::DataContainerQVector;

#pragma link C++ function Qn::ToTGraph;


#endif
