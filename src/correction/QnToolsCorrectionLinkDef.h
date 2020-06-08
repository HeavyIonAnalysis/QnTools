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

#pragma link C++ class Qn::CorrectionBase + ;
#pragma link C++ class Qn::CorrectionOnInputData + ;
#pragma link C++ class Qn::CorrectionOnQnVector + ;
#pragma link C++ class Qn::CorrectionsSetOnInputData + ;
#pragma link C++ class Qn::CorrectionsSetOnQvector + ;
#pragma link C++ class Qn::SubEvent + ;
#pragma link C++ class Qn::SubEventChannels + ;
#pragma link C++ class Qn::SubEventTracks + ;
#pragma link C++ class Qn::CorrectionAxis + ;
#pragma link C++ class Qn::CorrectionAxisSet + ;
#pragma link C++ class Qn::CorrectionHistogram + ;
#pragma link C++ class Qn::CorrectionHistogramBase + ;
#pragma link C++ class Qn::CorrectionHistogramChannelized + ;
#pragma link C++ class Qn::CorrectionHistogramChannelizedSparse + ;
#pragma link C++ class Qn::CorrectionHistogramSparse + ;
#pragma link C++ class Qn::GainEqualization + ;
#pragma link C++ class Qn::CorrectionCalculator + ;
#pragma link C++ class Qn::CorrectionProfile + ;
#pragma link C++ class Qn::CorrectionProfile3DCorrelations + ;
#pragma link C++ class Qn::CorrectionProfileChannelized + ;
#pragma link C++ class Qn::CorrectionProfileChannelizedIngress + ;
#pragma link C++ class Qn::CorrectionProfileComponents + ;
#pragma link C++ class Qn::CorrectionProfileCorrelationComponents + ;
#pragma link C++ class Qn::CorrectionProfileCorrelationComponentsHarmonics + ;
#pragma link C++ class Qn::Alignment + ;
#pragma link C++ class Qn::Recentering + ;
#pragma link C++ class Qn::TwistAndRescale + ;
#pragma link C++ class Qn::DataContainer < std::unique_ptr < Qn::SubEvent>, Qn::Axis < double>> + ;
#pragma link C++ class Qn::InputVariable + ;
#pragma link C++ class std::map < std::string, InputVariable> + ;
#pragma link C++ class Qn::CorrectionBase + ;
#pragma link C++ class Qn::InputVariableManager + ;
#pragma link C++ class Qn::Detector + ;
#pragma link C++ class Qn::DetectorList + ;
#pragma link C++ class Qn::RunList + ;
#pragma link C++ class Qn::QAHistogram + ;
#pragma link C++ class Qn::QAHistograms + ;
#pragma link C++ class Qn::CorrectionManager + ;

#endif
