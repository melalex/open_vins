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

#include "UpdaterWheel.h"

#include "state/State.h"
#include "state/StateHelper.h"

#include "utils/print.h"
#include "utils/quat_ops.h"

#include <boost/math/distributions/chi_squared.hpp>

using namespace ov_core;
using namespace ov_type;
using namespace ov_msckf;

UpdaterWheel::UpdaterWheel(UpdaterWheelOptions &options) : _options(options) {
  for (int i = 1; i <= 3; i++) {
    boost::math::chi_squared chi_squared_dist(i);
    chi_squared_table[i] = boost::math::quantile(chi_squared_dist, 0.95);
  }
}

void UpdaterWheel::feed_wheel(const ov_core::WheelData &message) {
  std::lock_guard<std::mutex> lck(wheel_data_mtx);

  // Append it to our vector, and keep it sorted if it came out of order
  wheel_data.emplace_back(message);
  if (wheel_data.size() > 1 && wheel_data.at(wheel_data.size() - 2).timestamp > message.timestamp) {
    std::sort(wheel_data.begin(), wheel_data.end());
  }

  // Bound the history in case no update is consuming it (e.g. the filter is not initialized yet)
  double oldest_time = wheel_data.back().timestamp - 10.0;
  auto it0 = wheel_data.begin();
  while (it0 != wheel_data.end() && it0->timestamp < oldest_time) {
    it0 = wheel_data.erase(it0);
  }
}

bool UpdaterWheel::preintegrate(const std::vector<ov_core::WheelData> &data, double time0, double time1, const Eigen::Vector3d &intr,
                                const UpdaterWheelOptions &options, Preintegration &result) {

  result = Preintegration();
  if (data.empty() || time1 <= time0) {
    return false;
  }

  // Reading at a time: interpolated between its neighbours, or the nearest one if it is close enough
  auto reading_at = [&](double time, Eigen::Vector2d &m) {
    auto next = std::lower_bound(data.begin(), data.end(), time,
                                 [](const ov_core::WheelData &d, double t) { return d.timestamp < t; });
    if (next == data.begin()) {
      m << next->m1, next->m2;
      return next->timestamp - time <= options.max_gap;
    }
    auto prev = std::prev(next);
    if (next == data.end()) {
      m << prev->m1, prev->m2;
      return time - prev->timestamp <= options.max_gap;
    }
    if (next->timestamp - prev->timestamp > options.max_gap) {
      return false;
    }
    double lambda = (time - prev->timestamp) / (next->timestamp - prev->timestamp);
    m << (1 - lambda) * prev->m1 + lambda * next->m1, (1 - lambda) * prev->m2 + lambda * next->m2;
    return true;
  };

  // The readings we will integrate: the interval ends, and every reading in between
  std::vector<std::pair<double, Eigen::Vector2d>> readings;
  Eigen::Vector2d m_end;
  readings.emplace_back(time0, Eigen::Vector2d::Zero());
  if (!reading_at(time0, readings.back().second) || !reading_at(time1, m_end)) {
    return false;
  }
  for (const auto &d : data) {
    if (d.timestamp > time0 && d.timestamp < time1) {
      readings.emplace_back(d.timestamp, Eigen::Vector2d(d.m1, d.m2));
    }
  }
  readings.emplace_back(time1, m_end);

  // Discrete noise of a reading
  Eigen::Matrix2d Q = Eigen::Matrix2d::Zero();
  if (options.type == UpdaterWheelOptions::WheelType::WHEEL2D_ANG) {
    Q.diagonal() << std::pow(options.noise_w, 2), std::pow(options.noise_w, 2);
  } else if (options.type == UpdaterWheelOptions::WheelType::WHEEL2D_LIN) {
    Q.diagonal() << std::pow(options.noise_v, 2), std::pow(options.noise_v, 2);
  } else {
    Q.diagonal() << std::pow(options.noise_w, 2), std::pow(options.noise_v, 2);
  }

  // Integrate one interval at a time
  // The measurement is [theta, x, y], and each interval has constant velocities (the mean of its two readings)
  for (size_t i = 0; i + 1 < readings.size(); i++) {
    double dt = readings.at(i + 1).first - readings.at(i).first;
    if (dt > options.max_gap) {
      return false;
    }
    Eigen::Vector2d m = 0.5 * (readings.at(i).second + readings.at(i + 1).second);

    // Angular and linear velocity of the odometry frame, and their jacobians
    // J_m is in respect to the readings, J_i is in respect to the intrinsics
    double w, v;
    Eigen::Matrix2d J_m = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 3> J_i = Eigen::Matrix<double, 2, 3>::Zero();
    if (options.type == UpdaterWheelOptions::WheelType::WHEEL2D_CEN) {
      w = m(0);
      v = m(1);
    } else {
      double rl = intr(0), rr = intr(1), b = intr(2);
      w = (m(1) * rr - m(0) * rl) / b;
      v = (m(1) * rr + m(0) * rl) / 2.0;
      J_m << -rl / b, rr / b, rl / 2.0, rr / 2.0;
      J_i << -m(0) / b, m(1) / b, -w / b, m(0) / 2.0, m(1) / 2.0, 0.0;
    }

    // Midpoint integration of the 2D kinematics
    double th_mid = result.z(0) + 0.5 * w * dt;
    double c = std::cos(th_mid), s = std::sin(th_mid);
    Eigen::Matrix3d Phi_tr = Eigen::Matrix3d::Identity();
    Phi_tr(1, 0) = -v * s * dt;
    Phi_tr(2, 0) = v * c * dt;
    Eigen::Matrix<double, 3, 2> G;
    G << dt, 0.0, -0.5 * v * s * dt * dt, c * dt, 0.5 * v * c * dt * dt, s * dt;

    result.z(0) += w * dt;
    result.z(1) += v * c * dt;
    result.z(2) += v * s * dt;
    result.P = Phi_tr * result.P * Phi_tr.transpose() + G * J_m * Q * J_m.transpose() * G.transpose();
    result.dz_dintr = Phi_tr * result.dz_dintr + G * J_i;
  }
  return true;
}

Eigen::Vector3d UpdaterWheel::predict(const Eigen::Matrix3d &R_GtoI0, const Eigen::Vector3d &p_I0inG, const Eigen::Matrix3d &R_GtoI1,
                                      const Eigen::Vector3d &p_I1inG, const Eigen::Matrix3d &R_ItoO, const Eigen::Vector3d &p_IinO) {

  // Rotation from the older to the newer odometry frame, a positive yaw of the frame is a negative rotation of it
  Eigen::Matrix3d R_O0toO1 = R_ItoO * R_GtoI1 * R_GtoI0.transpose() * R_ItoO.transpose();

  // Position of the newer odometry frame in the older one
  Eigen::Vector3d p_O1inO0 = R_ItoO * R_GtoI0 * (p_I1inG - p_I0inG) - R_O0toO1.transpose() * p_IinO + p_IinO;

  Eigen::Vector3d z;
  z << -log_so3(R_O0toO1)(2), p_O1inO0(0), p_O1inO0(1);
  return z;
}

UpdaterWheel::Jacobians UpdaterWheel::jacobians(const Eigen::Matrix3d &R_GtoI0, const Eigen::Vector3d &p_I0inG,
                                                const Eigen::Matrix3d &R_GtoI1, const Eigen::Vector3d &p_I1inG,
                                                const Eigen::Matrix3d &R_ItoO, const Eigen::Vector3d &p_IinO) {

  Eigen::Matrix3d R_O0toO1 = R_ItoO * R_GtoI1 * R_GtoI0.transpose() * R_ItoO.transpose();
  Eigen::Matrix3d R_O1toO0 = R_O0toO1.transpose();
  Eigen::Vector3d dp_inO0 = R_ItoO * R_GtoI0 * (p_I1inG - p_I0inG);
  Eigen::Matrix<double, 1, 3> e3 = Eigen::Vector3d::UnitZ().transpose();
  Eigen::Matrix<double, 2, 3> Lambda = Eigen::Matrix3d::Identity().topRows(2);

  Jacobians J;
  J.H_pose1.setZero();
  J.H_pose0.setZero();
  J.H_ext.setZero();

  // Newer IMU pose
  J.H_pose1.block(0, 0, 1, 3) = e3 * R_ItoO;
  J.H_pose1.block(1, 0, 2, 3) = Lambda * R_O1toO0 * skew_x(p_IinO) * R_ItoO;
  J.H_pose1.block(1, 3, 2, 3) = Lambda * R_ItoO * R_GtoI0;

  // Older IMU pose
  J.H_pose0.block(0, 0, 1, 3) = -e3 * R_O0toO1 * R_ItoO;
  J.H_pose0.block(1, 0, 2, 3) =
      Lambda * (R_ItoO * skew_x(R_GtoI0 * (p_I1inG - p_I0inG)) - R_O1toO0 * skew_x(p_IinO) * R_O0toO1 * R_ItoO);
  J.H_pose0.block(1, 3, 2, 3) = -Lambda * R_ItoO * R_GtoI0;

  // Extrinsics
  J.H_ext.block(0, 0, 1, 3) = e3 * (Eigen::Matrix3d::Identity() - R_O0toO1);
  J.H_ext.block(1, 0, 2, 3) = Lambda * (skew_x(dp_inO0) + R_O1toO0 * skew_x(p_IinO) * (Eigen::Matrix3d::Identity() - R_O0toO1));
  J.H_ext.block(1, 3, 2, 3) = Lambda * (Eigen::Matrix3d::Identity() - R_O1toO0);
  return J;
}

bool UpdaterWheel::update(std::shared_ptr<State> state, const Eigen::Vector3d &w_clone) {

  // Remember the velocities at this clone, and forget the ones of the clones that are gone
  double time1 = state->_timestamp;
  if (state->_clones_IMU.find(time1) == state->_clones_IMU.end()) {
    return false;
  }
  clone_kinematics[time1] = {w_clone, state->_imu->vel()};
  for (auto it = clone_kinematics.begin(); it != clone_kinematics.end();) {
    if (state->_clones_IMU.find(it->first) == state->_clones_IMU.end()) {
      it = clone_kinematics.erase(it);
    } else {
      it++;
    }
  }

  // Find the clone before the newest one
  double time0 = -1;
  for (const auto &clone : state->_clones_IMU) {
    if (clone.first < time1 && clone.first > time0) {
      time0 = clone.first;
    }
  }
  if (time0 < 0 || clone_kinematics.find(time0) == clone_kinematics.end()) {
    return false;
  }
  if (time1 - time0 > _options.max_dt) {
    _num_skipped++;
    return false;
  }

  // Integrate the wheel readings between the two clones
  // The clone times are in the camera clock, the wheel readings are in the wheel clock: t_imu = t_wheel + dt_wheel
  Eigen::Vector3d intr = state->_calib_wheel_intrinsics->value();
  double t_off_cam = state->_calib_dt_CAMtoIMU->value()(0);
  double t_off_wheel = state->_calib_dt_WHEELtoIMU->value()(0);
  double time0_wheel = time0 + t_off_cam - t_off_wheel;
  double time1_wheel = time1 + t_off_cam - t_off_wheel;
  Preintegration preint;
  {
    std::lock_guard<std::mutex> lck(wheel_data_mtx);
    bool success = preintegrate(wheel_data, time0_wheel, time1_wheel, intr, _options, preint);

    // The next interval starts where this one ends, so we only need the readings around its start
    auto it0 = wheel_data.begin();
    while (it0 != wheel_data.end() && std::next(it0) != wheel_data.end() &&
           std::next(it0)->timestamp < time1_wheel - _options.max_gap) {
      it0 = wheel_data.erase(it0);
    }
    if (!success) {
      _num_skipped++;
      return false;
    }
  }

  // The two clones and the extrinsics
  std::shared_ptr<PoseJPL> clone0 = state->_clones_IMU.at(time0);
  std::shared_ptr<PoseJPL> clone1 = state->_clones_IMU.at(time1);
  std::shared_ptr<PoseJPL> calib = state->_calib_IMUtoWHEEL;

  // Residual with the current estimates, jacobians with the first estimates if we are doing FEJ
  Eigen::Vector3d z_pred = predict(clone0->Rot(), clone0->pos(), clone1->Rot(), clone1->pos(), calib->Rot(), calib->pos());
  Eigen::Vector3d res = preint.z - z_pred;
  res(0) = std::atan2(std::sin(res(0)), std::cos(res(0)));
  bool fej = state->_options.do_fej;
  Jacobians J = jacobians(fej ? clone0->Rot_fej() : clone0->Rot(), fej ? clone0->pos_fej() : clone0->pos(),
                          fej ? clone1->Rot_fej() : clone1->Rot(), fej ? clone1->pos_fej() : clone1->pos(),
                          fej ? calib->Rot_fej() : calib->Rot(), fej ? calib->pos_fej() : calib->pos());

  // Stack the jacobian of every variable that is in the state
  std::vector<std::shared_ptr<Type>> Hx_order = {clone1, clone0};
  int H_size = 12;
  if (state->_options.do_calib_wheel_extrinsics) {
    Hx_order.push_back(calib);
    H_size += calib->size();
  }
  if (state->_options.do_calib_wheel_timeoffset) {
    Hx_order.push_back(state->_calib_dt_WHEELtoIMU);
    H_size += 1;
  }
  if (state->_options.do_calib_wheel_intrinsics) {
    Hx_order.push_back(state->_calib_wheel_intrinsics);
    H_size += 3;
  }
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, H_size);
  H.block(0, 0, 3, 6) = J.H_pose1;
  H.block(0, 6, 3, 6) = J.H_pose0;
  int H_id = 12;
  if (state->_options.do_calib_wheel_extrinsics) {
    H.block(0, H_id, 3, 6) = J.H_ext;
    H_id += 6;
  }
  if (state->_options.do_calib_wheel_timeoffset) {
    // An error in the time offset moves both poses along their velocities
    const CloneKinematics &kin0 = clone_kinematics.at(time0);
    const CloneKinematics &kin1 = clone_kinematics.at(time1);
    H.block(0, H_id, 3, 1) = J.H_pose1.block(0, 0, 3, 3) * kin1.w + J.H_pose1.block(0, 3, 3, 3) * kin1.v +
                             J.H_pose0.block(0, 0, 3, 3) * kin0.w + J.H_pose0.block(0, 3, 3, 3) * kin0.v;
    H_id += 1;
  }
  if (state->_options.do_calib_wheel_intrinsics) {
    H.block(0, H_id, 3, 3) = -preint.dz_dintr;
    H_id += 3;
  }

  // Measurement noise: the integrated reading noise, plus what the kinematic model does not capture
  Eigen::Matrix3d R = preint.P;
  double sigma_turn = _options.noise_turn * std::abs(z_pred(0));
  R(0, 0) += std::pow(_options.noise_floor_yaw, 2);
  R(1, 1) += std::pow(_options.noise_floor_x, 2) + std::pow(sigma_turn, 2);
  R(2, 2) += std::pow(_options.noise_floor_y, 2) + std::pow(sigma_turn, 2);

  // Drop the yaw if we should not use it
  Eigen::MatrixXd H_used = H;
  Eigen::VectorXd res_used = res;
  Eigen::MatrixXd R_used = R;
  if (!_options.use_yaw) {
    H_used = H.bottomRows(2);
    res_used = res.tail(2);
    R_used = R.bottomRightCorner(2, 2);
  }

  // Chi2 distance check
  Eigen::MatrixXd P_marg = StateHelper::get_marginal_covariance(state, Hx_order);
  Eigen::MatrixXd S = H_used * P_marg * H_used.transpose() + R_used;
  double chi2 = res_used.dot(S.llt().solve(res_used));
  if (chi2 > _options.chi2_multipler * chi_squared_table[(int)res_used.rows()]) {
    _num_rejected++;
    return false;
  }

  // We are good!!! Update the state
  StateHelper::EKFUpdate(state, Hx_order, H_used, res_used, R_used);
  _num_updates++;
  return true;
}
