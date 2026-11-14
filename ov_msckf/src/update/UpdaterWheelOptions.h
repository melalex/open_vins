/*
 * OpenVINS: An Open Platform for Visual-Inertial Research
 * Copyright (C) 2018-2023 Patrick Geneva
 * Copyright (C) 2018-2023 Guoquan Huang
 * Copyright (C) 2018-2023 OpenVINS Contributors
 * Copyright (C) 2018-2019 Kevin Eckenhoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef OV_MSCKF_UPDATER_WHEEL_OPTIONS_H
#define OV_MSCKF_UPDATER_WHEEL_OPTIONS_H

#include <string>

#include "utils/print.h"

namespace ov_msckf {

/**
 * @brief Struct which stores the wheel odometry updater options
 */
struct UpdaterWheelOptions {

  /**
   * @brief What the two values of a wheel measurement are
   *
   * - WHEEL2D_ANG: left and right wheel angular velocities (rad/s). Intrinsics: left radius, right radius, baseline.
   * - WHEEL2D_LIN: left and right wheel linear velocities (m/s). Intrinsics: left scale, right scale, baseline.
   * - WHEEL2D_CEN: angular (rad/s) and linear (m/s) velocity of the odometry frame. No intrinsics.
   */
  enum WheelType { WHEEL2D_ANG, WHEEL2D_LIN, WHEEL2D_CEN };

  /// Type of the wheel measurements
  WheelType type = WheelType::WHEEL2D_ANG;

  /// White noise sigma of an angular velocity reading (rad/s)
  double noise_w = 0.01;

  /// White noise sigma of a linear velocity reading (m/s)
  double noise_v = 0.01;

  /// What chi-squared multipler we should apply
  double chi2_multipler = 1;

  /// If the relative yaw of the preintegrated measurement is used in the update.
  /// Disable for vehicles whose wheel speed difference is a poor yaw source (e.g. skid-steering).
  bool use_yaw = true;

  /// Sigma added to the preintegrated relative yaw (rad), models what the kinematic model does not
  double noise_floor_yaw = 0.0;

  /// Sigma added to the preintegrated forward displacement (m)
  double noise_floor_x = 0.0;

  /// Sigma added to the preintegrated lateral displacement (m), models side slip and non-planar motion
  double noise_floor_y = 0.0;

  /// Sigma added to the preintegrated displacement per radian turned over the interval (m/rad), models slip in turns
  double noise_turn = 0.0;

  /// Max time between two clones we will preintegrate over (s)
  double max_dt = 1.0;

  /// Max time between two wheel readings (or a reading and an interval end) we will integrate over (s)
  double max_gap = 0.5;

  /// Parses the wheel type from a string
  static bool type_from_string(const std::string &str, WheelType &type) {
    if (str == "Wheel2DAng") {
      type = WheelType::WHEEL2D_ANG;
    } else if (str == "Wheel2DLin") {
      type = WheelType::WHEEL2D_LIN;
    } else if (str == "Wheel2DCen") {
      type = WheelType::WHEEL2D_CEN;
    } else {
      return false;
    }
    return true;
  }

  /// Nice print function of what parameters we have loaded
  void print() {
    PRINT_DEBUG("    - wheel_type: %d\n", (int)type);
    PRINT_DEBUG("    - wheel_noise_w: %.4f\n", noise_w);
    PRINT_DEBUG("    - wheel_noise_v: %.4f\n", noise_v);
    PRINT_DEBUG("    - wheel_chi2_multipler: %.1f\n", chi2_multipler);
    PRINT_DEBUG("    - wheel_use_yaw: %d\n", use_yaw);
    PRINT_DEBUG("    - wheel_noise_floor (yaw, x, y): %.4f, %.4f, %.4f\n", noise_floor_yaw, noise_floor_x, noise_floor_y);
    PRINT_DEBUG("    - wheel_noise_turn: %.4f\n", noise_turn);
    PRINT_DEBUG("    - wheel_max_dt: %.3f\n", max_dt);
    PRINT_DEBUG("    - wheel_max_gap: %.3f\n", max_gap);
  }
};

} // namespace ov_msckf

#endif // OV_MSCKF_UPDATER_WHEEL_OPTIONS_H
