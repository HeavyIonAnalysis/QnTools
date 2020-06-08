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

#ifndef FLOW_INPUTVARIABLEMANAGER_H
#define FLOW_INPUTVARIABLEMANAGER_H

#include <string>
#include <map>
#include <utility>
#include <cmath>
#include "TTree.h"

#include "InputVariable.hpp"

/**
 * @brief Attaches a variable to a chosen tree.
 * Type is converted to the chosen type.
 * @tparam T type of variable saved in the tree.
 */
template<typename T>
class OutputValue {
 public:
  explicit OutputValue(Qn::InputVariable *var) : var_(var) {}
  /**
   * @brief updates the value. To be called every event.
   */
  void UpdateValue() { value_ = (T) (*var_->begin()); }
  /**
   * @brief Creates a new branch in the tree.
   * @param tree output tree
   */
  void SetToTree(TTree *tree) {
    tree->Branch(var_->GetName().data(), &value_);
  }
 private:
  T value_; /// value which is written
  Qn::InputVariable *var_; /// Variable to be written to the tree
};

namespace Qn {
/**
 * @brief Manages the input variables for the correction step.
 * A variable consist of a name and a unsigned integer position in the value array and an unsigned integer length.
 * Different variables all access the same variable container.
 */
class InputVariableManager {
 public:

  using f_type = Double32_t;
  using i_type = Long64_t;

  InputVariableManager() {
    CreateVariableOnes();
  }

  virtual ~InputVariableManager() {
    for (auto &var : channel_variables_) {
      delete[] var;
    }
  }

  void Initialize() {
    variable_values_float_ = new f_type[kMaxSize]; /// non-owning pointer to variables
    variable_values_ones_ = new f_type[kMaxSize]; /// values container of ones.

    for (auto &var : variable_map_) {
      var.second.values_container_ = variable_values_float_;
    }
    variable_map_["Ones"].values_container_ = variable_values_ones_;
    for (int i = 0; i < kMaxSize; ++i) { variable_values_float_[i] = NAN; }
    for (unsigned int i = 0; i < kMaxSize; ++i) { variable_values_ones_[i] = 1.0; }
  }

  void InitVariable(InputVariable &var) {
    if (var.name_=="Ones") {
      var.values_container_ = variable_values_ones_;
    } else {
      var.values_container_ = variable_values_float_;
    }
  }

  /**
   * @brief Creates a new variable in the variable manager.
   * @param name Name of the new variable.
   * @param id position of the variable in the values container.
   * @param length length of the variable in the values container.
   */
  void CreateVariable(std::string name, const int id, const int length) {
    InputVariable var(id, length, name);
    variable_map_.emplace(name, var);
  }
  /**
   * @brief Initializes the variable container for ones.
   */
  void CreateVariableOnes() {
    InputVariable var(0, kMaxSize, "Ones");
    variable_map_.emplace("Ones", var);
  }
  /**
   * @brief Creates a channel variable (variables which counts from 0 to the size-1).
   * @param name name of the variable
   * @param size number of channels
   */
  void CreateChannelVariable(const std::string &name, const int size) {
    InputVariable var(0, size, name);
    auto *arr = new f_type[size];
    for (int i = 0; i < size; ++i) { arr[i] = i; }
    var.values_container_ = arr;
    variable_map_.emplace(name, var);
    channel_variables_.push_back(arr);
  }
  /**
   * @brief Finds the variable in the variable manager.
   * @param name Name of the variable.
   * @return Variable of the given name.
   */
  InputVariable FindVariable(const std::string &name) const {
    try{
      return variable_map_.at(name);
    } catch (std::out_of_range&){
      throw std::out_of_range("Requested variable " + name + " cannot be found.");
    }
  }
  /**
 * @brief Finds the variable in the variable manager.
 * @param name Name of the variable.
 * @return Variable of the given name.
 */
  InputVariable *FindVariablePtr(const std::string &name) { return &variable_map_.at(name); }
  /**
   * @brief Find the position in the values container of a variable with a given name.
   * @param name Name of the variable.
   * @return position in the values container.
   */
  int FindNum(const std::string &name) const { return variable_map_.at(name).id_; }
  /**
   * @brief Get the values container.
   * @return a pointer to the values container.
   */
  f_type *GetVariableContainer() { return variable_values_float_; }
  /**
   * @brief Register the variable to be saved in the output tree as float.
   * Beware of conversion from double.
   * @param name Name of the variable.
   */
  void RegisterOutputF(const std::string &name) {
    OutputValue<f_type> val(FindVariablePtr(name));
    variable_output_float_.push_back(val);
  }

  /**
   * @brief Register the variable to be saved in the output tree as long.
   * Beware of conversion from double.
   * @param name Name of the variable.
   */
  void RegisterOutputL(const std::string &name) {
    OutputValue<i_type> val(FindVariablePtr(name));
    variable_output_integer_.push_back(val);
  }

  /**
   * @brief Creates Branches in tree for saving the event information.
   * @param tree output tree to contain the event information
   */
  void SetOutputTree(TTree *tree) {
    for (auto &element : variable_output_float_) { element.SetToTree(tree); }
    for (auto &element : variable_output_integer_) { element.SetToTree(tree); }
  }

  /**
   * @brief Updates the output variables.
   */
  void UpdateOutVariables() {
    for (auto &element : variable_output_float_) { element.UpdateValue(); }
    for (auto &element : variable_output_integer_) { element.UpdateValue(); }
  }

 private:
  static constexpr int kMaxSize = 11000; /// Maximum number of variables.
  f_type *variable_values_float_ = nullptr; //!<! non-owning pointer to variables
  f_type *variable_values_ones_ = nullptr; //!<! values container of ones.
  std::vector<double *> channel_variables_; //!<!
  std::map<std::string, InputVariable> variable_map_; /// name to variable map
  std::vector<OutputValue<f_type>> variable_output_float_; //!<! variables registered for output as float
  std::vector<OutputValue<i_type>> variable_output_integer_; //!<! variables registered for output as long
  /// \cond CLASSIMP
 ClassDef(InputVariableManager, 1);
/// \endcond
};
}

#endif //FLOW_INPUTVARIABLEMANAGER_H
