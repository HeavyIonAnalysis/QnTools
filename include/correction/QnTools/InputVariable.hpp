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

#ifndef QN_INPUTVARIABLE
#define QN_INPUTVARIABLE

namespace Qn {
/**
 * Variable
 * Constructor is private so it is only created in the variable manager
 * with the correct pointer to the values container.
 */
class InputVariable {
 private:
  InputVariable(const unsigned int id, const unsigned int length, std::string name) :
      id_(id),
      size_(length),
      name_(std::move(name)) {}
  InputVariable(const unsigned int id, const unsigned int length) : id_(id), size_(length) {}
  unsigned int id_{0}; /// position in the values container
  unsigned int size_{0}; /// length in the  values container
  double *values_container_ = nullptr; //!<! pointer to the values container
  std::string name_; /// name of the variable
  friend class InputVariableManager;
  friend class CorrectionCuts;
 public:
  InputVariable() = default;
  virtual ~InputVariable() = default;
  inline double operator[](std::size_t i) const noexcept { return values_container_[id_ + i];}
  inline double *begin() noexcept { return &values_container_[id_]; } /// implements begin for iteration
  inline double *end() noexcept { return &values_container_[id_ + size_]; } /// implements end for iteration
  const double *begin() const noexcept { return &values_container_[id_]; }  /// implements begin for iteration
  inline double *end() const noexcept { return &values_container_[id_ + size_]; }  /// implements end for iteration
  inline double *at(std::size_t i) noexcept { return &values_container_[id_ + i]; }
  const inline double *Get() const noexcept { return &values_container_[id_]; }
  inline unsigned int size() const noexcept { return size_; }
  inline unsigned int GetID() const noexcept { return id_; }
  std::string GetName() const { return name_; }

  /// \cond CLASSIMP
 ClassDef(InputVariable, 1);
/// \endcond
};
}

#endif