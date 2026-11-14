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

#ifndef OV_MSCKF_UPDATER_WHEEL_H
#define OV_MSCKF_UPDATER_WHEEL_H

#include <Eigen/Eigen>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

#include "utils/sensor_data.h"

#include "UpdaterWheelOptions.h"

namespace ov_msckf {

class State;

/**
 * @brief Will update the two newest clones with preintegrated 2D wheel odometry.
 *
 * The wheel readings between the two newest clones are integrated into a single 2D relative pose
 * (yaw and planar displacement of the odometry frame), which is then used as an EKF update.
 * The wheel intrinsics, the wheel-to-IMU extrinsics and the time offset can be calibrated online.
 * This follows:
 *
 * > W. Lee, K. Eckenhoff, Y. Yang, P. Geneva and G. Huang,
 * > "Visual-Inertial-Wheel Odometry with Online Calibration", IROS 2020.
 *
 * The odometry frame {O} has x forward, y left and z up, and the wheels move it along its x axis only.
 */
class UpdaterWheel {

public:
  /// A preintegrated wheel measurement
  struct Preintegration {
    /// Relative yaw and position of the end frame in the start frame [theta, x, y]
    Eigen::Vector3d z = Eigen::Vector3d::Zero();
    /// Covariance of the measurement due to the wheel reading noise
    Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
    /// Jacobian of the measurement in respect to the intrinsics
    Eigen::Matrix3d dz_dintr = Eigen::Matrix3d::Zero();
  };

  /// Jacobians of the predicted measurement
  struct Jacobians {
    /// In respect to the newer IMU pose [theta, p]
    Eigen::Matrix<double, 3, 6> H_pose1;
    /// In respect to the older IMU pose [theta, p]
    Eigen::Matrix<double, 3, 6> H_pose0;
    /// In respect to the extrinsics [theta_ItoO, p_IinO]
    Eigen::Matrix<double, 3, 6> H_ext;
  };

  /**
   * @brief Default constructor for our wheel updater
   * @param options Updater options (noises, measurement type, etc)
   */
  UpdaterWheel(UpdaterWheelOptions &options);

  /**
   * @brief Feed function for wheel data
   * @param message Contains our timestamp and wheel information
   */
  void feed_wheel(const ov_core::WheelData &message);

  /**
   * @brief Will try to update the two newest clones with the wheel readings between them
   * @param state State of the filter (it should have just been cloned at its current time)
   * @param w_clone Angular velocity in the IMU frame at the state's current time
   * @return True if the state was updated
   */
  bool update(std::shared_ptr<State> state, const Eigen::Vector3d &w_clone);

  /**
   * @brief Integrates the wheel readings from time0 to time1 (in the wheel clock)
   * @param data Wheel readings sorted by time
   * @param time0 Start of the interval
   * @param time1 End of the interval
   * @param intr Wheel intrinsics (see UpdaterWheelOptions::WheelType)
   * @param options Updater options
   * @param[out] result Preintegrated measurement
   * @return False if the readings do not cover the interval
   */
  static bool preintegrate(const std::vector<ov_core::WheelData> &data, double time0, double time1, const Eigen::Vector3d &intr,
                           const UpdaterWheelOptions &options, Preintegration &result);

  /**
   * @brief Relative 2D pose of the odometry frame the state predicts
   * @param R_GtoI0 Orientation of the older IMU pose
   * @param p_I0inG Position of the older IMU pose
   * @param R_GtoI1 Orientation of the newer IMU pose
   * @param p_I1inG Position of the newer IMU pose
   * @param R_ItoO Rotation from the IMU to the odometry frame
   * @param p_IinO Position of the IMU in the odometry frame
   * @return Relative yaw and position of the newer odometry frame in the older [theta, x, y]
   */
  static Eigen::Vector3d predict(const Eigen::Matrix3d &R_GtoI0, const Eigen::Vector3d &p_I0inG, const Eigen::Matrix3d &R_GtoI1,
                                 const Eigen::Vector3d &p_I1inG, const Eigen::Matrix3d &R_ItoO, const Eigen::Vector3d &p_IinO);

  /**
   * @brief Jacobians of predict() in respect to the error states of its arguments
   */
  static Jacobians jacobians(const Eigen::Matrix3d &R_GtoI0, const Eigen::Vector3d &p_I0inG, const Eigen::Matrix3d &R_GtoI1,
                             const Eigen::Vector3d &p_I1inG, const Eigen::Matrix3d &R_ItoO, const Eigen::Vector3d &p_IinO);

  /// Number of updates applied
  int num_updates() const { return _num_updates; }

  /// Number of updates rejected by the chi-squared test
  int num_rejected() const { return _num_rejected; }

  /// Number of clone pairs skipped (no wheel readings over the interval, or too long an interval)
  int num_skipped() const { return _num_skipped; }

protected:
  /// Angular and linear velocity of the IMU at a clone
  struct CloneKinematics {
    Eigen::Vector3d w;
    Eigen::Vector3d v;
  };

  /// Options used during update
  UpdaterWheelOptions _options;

  /// Our history of wheel readings (sorted by time)
  std::vector<ov_core::WheelData> wheel_data;
  std::mutex wheel_data_mtx;

  /// Velocities at the clone times (used for the time offset jacobian)
  std::map<double, CloneKinematics> clone_kinematics;

  /// Chi squared 95th percentile table (lookup would be size of residual)
  std::map<int, double> chi_squared_table;

  int _num_updates = 0;
  int _num_rejected = 0;
  int _num_skipped = 0;
};

} // namespace ov_msckf

#endif // OV_MSCKF_UPDATER_WHEEL_H
