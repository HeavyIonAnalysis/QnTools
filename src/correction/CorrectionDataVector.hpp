#ifndef QN_DATAVECTOR_H
#define QN_DATAVECTOR_H

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

namespace Qn {
/**
 * @class CorrectionDataVector
 * @brief Class models a data vector.
 * It allows to model data vector of different detector types.
 */
class CorrectionDataVector {
 public:

  /**
   * Default constructor
   */
  CorrectionDataVector() = default;

  /**
   * Constructor
   * @param id id of the channel
   * @param phi azimuthal angle of the channel or track
   * @param weight weight applied to the channel or track
   * @param radial_offset radial offset of the channel is only used for certain detector geometries.
   */
  CorrectionDataVector(int id, float phi, float weight, float radial_offset) :
      id_(id),
      phi_(phi),
      radial_offset_(radial_offset),
      weight_(weight),
      equalized_weight_(weight) {}

  ~CorrectionDataVector() = default;

  /**
   * Sets the equalized weight
   * @param weight equalized weight after channel equalization
   */
  inline void SetEqualizedWeight(const float weight) { equalized_weight_ = weight; }
  /**
   * Gets the channel id associated with the data vector
   * @return the channel id
   */
  constexpr int GetId() const { return id_; }

  /**
   * Gets the azimuthal angle for the data vector
   * @return phi
   */
  constexpr float Phi() const { return phi_; }
  /**
   * Gets the radial offset of the data vector
   * @return radial offset
   */
  constexpr float RadialOffset() const {return radial_offset_; }
  /**
   * Gets the weight for the data vector
   * @return defaults to 1.0
   */
  constexpr float Weight() const { return weight_; }
  /**
   * Gets the equalized weight for the data vector
   * @return defaults to weights
   */
  constexpr float EqualizedWeight() const { return equalized_weight_; }
 private:
  int   id_;                    //!<! the id associated with the data vector
  float phi_;                   //!<! the azimuthal angle of the data vector
  float radial_offset_; //!<! radial offset of the channel represented by the data vector
  float weight_;                //!<! raw weight assigned to the data vector
  float equalized_weight_;      //!<! eq weight assigned to the data vector
};
}
#endif /* QNCORRECTIONS_DATAVECTORS_H */
