// Qn Tools
//
// Copyright (C) 2019  Lukas Kreis, Ilya Selyuzhenkov
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

#ifndef FLOW_RUNLIST_H
#define FLOW_RUNLIST_H

#include <string>
#include <utility>
#include <vector>
#include <Rtypes.h>

namespace Qn {
class RunList {
 public:
  explicit RunList(std::vector<std::string> list) : run_list_(std::move(list)) {}
  RunList() = default;
  virtual ~RunList() = default;
  std::vector<std::string>::const_iterator begin() const { return run_list_.begin(); }
  std::vector<std::string>::const_iterator end() const { return run_list_.end(); }
  void SetCurrentRun(std::string name) {
    current_run_name_ = std::move(name);
    if (std::find(run_list_.begin(), run_list_.end(), current_run_name_)==run_list_.end()) {
      run_list_.emplace_back(current_run_name_);
    }
  }
  std::string GetCurrent() const { return current_run_name_; }
  bool empty() const { return run_list_.empty(); }
 private:
  std::string current_run_name_;
  std::vector<std::string> run_list_;

  /// \cond CLASSIMP
 ClassDef(RunList, 1);
  /// \endcond
};
}
#endif