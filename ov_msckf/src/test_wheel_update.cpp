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

#include <cmath>
#include <iostream>
#include <random>

#include "update/UpdaterWheel.h"
#include "utils/quat_ops.h"

using namespace ov_core;
using namespace ov_msckf;

// Checks the wheel odometry measurement model:
// 1. the preintegration against the closed form of a constant-velocity arc
// 2. the preintegration jacobian in respect to the intrinsics against finite differences
// 3. the predicted measurement against the preintegrated one, for a 3D IMU trajectory of a planar vehicle
// 4. the predicted measurement jacobians against finite differences
int main() {

  int failed = 0;
  auto check = [&](const std::string &name, double error, double tolerance) {
    bool ok = error < tolerance;
    std::cout << (ok ? "[ OK ] " : "[FAIL] ") << name << " (error " << error << ", tolerance " << tolerance << ")" << std::endl;
    if (!ok)
      failed++;
  };

  std::mt19937 gen(0);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  auto random_vec = [&](double scale) { return Eigen::Vector3d(scale * uni(gen), scale * uni(gen), scale * uni(gen)); };

  //===================================================================================
  // 1. Constant wheel speeds drive a circular arc
  //===================================================================================
  UpdaterWheelOptions options;
  options.type = UpdaterWheelOptions::WheelType::WHEEL2D_ANG;
  Eigen::Vector3d intr(0.30, 0.31, 1.5);
  double wl = 3.0, wr = 4.0;
  std::vector<WheelData> data;
  for (int i = -2; i < 60; i++) {
    data.push_back({0.01 * i, wl, wr});
  }
  double time0 = 0.003, time1 = 0.487;
  UpdaterWheel::Preintegration preint;
  bool success = UpdaterWheel::preintegrate(data, time0, time1, intr, options, preint);
  double w = (wr * intr(1) - wl * intr(0)) / intr(2);
  double v = (wr * intr(1) + wl * intr(0)) / 2.0;
  double th = w * (time1 - time0);
  Eigen::Vector3d z_arc(th, v / w * std::sin(th), v / w * (1 - std::cos(th)));
  check("preintegration covers the interval", success ? 0.0 : 1.0, 0.5);
  check("preintegration of an arc", (preint.z - z_arc).norm(), 1e-6);

  //===================================================================================
  // 2. Jacobian in respect to the intrinsics (varying wheel speeds)
  //===================================================================================
  data.clear();
  for (int i = -2; i < 60; i++) {
    data.push_back({0.01 * i, 3.0 + std::sin(0.3 * i), 4.0 + std::cos(0.2 * i)});
  }
  UpdaterWheel::preintegrate(data, time0, time1, intr, options, preint);
  Eigen::Matrix3d dz_dintr_num;
  for (int i = 0; i < 3; i++) {
    Eigen::Vector3d intr_p = intr, intr_m = intr;
    intr_p(i) += 1e-6;
    intr_m(i) -= 1e-6;
    UpdaterWheel::Preintegration preint_p, preint_m;
    UpdaterWheel::preintegrate(data, time0, time1, intr_p, options, preint_p);
    UpdaterWheel::preintegrate(data, time0, time1, intr_m, options, preint_m);
    dz_dintr_num.col(i) = (preint_p.z - preint_m.z) / 2e-6;
  }
  check("jacobian in respect to the intrinsics", (preint.dz_dintr - dz_dintr_num).norm(), 1e-6);

  // Readings that do not cover the interval should fail
  options.max_gap = 0.05;
  check("missing readings are refused", UpdaterWheel::preintegrate(data, 0.5, 0.9, intr, options, preint) ? 1.0 : 0.0, 0.5);
  options.max_gap = 0.5;

  //===================================================================================
  // 3. A planar vehicle with a tilted, offset IMU: the state predicts what the wheels measure
  //===================================================================================
  Eigen::Matrix3d R_ItoO = exp_so3(Eigen::Vector3d(0.3, -0.2, 1.2));
  Eigen::Vector3d p_IinO(0.4, -0.1, 0.8);
  // Odometry frame poses in a tilted global plane: R_GtoO = Rz(-yaw) * R_GtoP, p_OinG = R_GtoP^T * [x y 0]
  Eigen::Matrix3d R_GtoP = exp_so3(Eigen::Vector3d(0.1, 0.05, 0.4));
  auto imu_pose = [&](double yaw, double x, double y, Eigen::Matrix3d &R_GtoI, Eigen::Vector3d &p_IinG) {
    Eigen::Matrix3d R_GtoO = exp_so3(Eigen::Vector3d(0, 0, -yaw)) * R_GtoP;
    Eigen::Vector3d p_OinG = R_GtoP.transpose() * Eigen::Vector3d(x, y, 0);
    R_GtoI = R_ItoO.transpose() * R_GtoO;
    p_IinG = p_OinG + R_GtoO.transpose() * p_IinO;
  };
  Eigen::Matrix3d R_GtoI0, R_GtoI1;
  Eigen::Vector3d p_I0inG, p_I1inG;
  double yaw0 = 0.7, x0 = 2.0, y0 = -1.0;
  imu_pose(yaw0, x0, y0, R_GtoI0, p_I0inG);
  // The second pose is the first one moved along the arc of test 1
  double x1 = x0 + std::cos(yaw0) * z_arc(1) - std::sin(yaw0) * z_arc(2);
  double y1 = y0 + std::sin(yaw0) * z_arc(1) + std::cos(yaw0) * z_arc(2);
  imu_pose(yaw0 + z_arc(0), x1, y1, R_GtoI1, p_I1inG);
  Eigen::Vector3d z_pred = UpdaterWheel::predict(R_GtoI0, p_I0inG, R_GtoI1, p_I1inG, R_ItoO, p_IinO);
  check("predicted measurement of a planar motion", (z_pred - z_arc).norm(), 1e-9);

  //===================================================================================
  // 4. Jacobians of the predicted measurement (general 3D poses, small relative rotation)
  //===================================================================================
  // Error states: R = exp(-dtheta) * R_hat, p = p_hat + dp
  double max_error_yaw = 0, max_error_pos = 0;
  for (int n = 0; n < 100; n++) {
    Eigen::Matrix3d R0 = exp_so3(random_vec(3.0));
    Eigen::Vector3d p0 = random_vec(10.0);
    Eigen::Matrix3d R1 = exp_so3(random_vec(0.02)) * R0;
    Eigen::Vector3d p1 = p0 + random_vec(0.5);
    Eigen::Matrix3d Ro = exp_so3(random_vec(3.0));
    Eigen::Vector3d po = random_vec(1.0);
    UpdaterWheel::Jacobians J = UpdaterWheel::jacobians(R0, p0, R1, p1, Ro, po);
    Eigen::Matrix<double, 3, 6> H1_num, H0_num, He_num;
    double eps = 1e-6;
    for (int i = 0; i < 6; i++) {
      Eigen::Vector3d dth = Eigen::Vector3d::Zero(), dp = Eigen::Vector3d::Zero();
      if (i < 3)
        dth(i) = eps;
      else
        dp(i - 3) = eps;
      Eigen::Matrix3d dRp = exp_so3(-dth), dRm = exp_so3(dth);
      H1_num.col(i) = (UpdaterWheel::predict(R0, p0, dRp * R1, p1 + dp, Ro, po) -
                       UpdaterWheel::predict(R0, p0, dRm * R1, p1 - dp, Ro, po)) / (2 * eps);
      H0_num.col(i) = (UpdaterWheel::predict(dRp * R0, p0 + dp, R1, p1, Ro, po) -
                       UpdaterWheel::predict(dRm * R0, p0 - dp, R1, p1, Ro, po)) / (2 * eps);
      He_num.col(i) = (UpdaterWheel::predict(R0, p0, R1, p1, dRp * Ro, po + dp) -
                       UpdaterWheel::predict(R0, p0, R1, p1, dRm * Ro, po - dp)) / (2 * eps);
    }
    for (const auto &diff : {J.H_pose1 - H1_num, J.H_pose0 - H0_num, J.H_ext - He_num}) {
      max_error_yaw = std::max(max_error_yaw, diff.topRows(1).cwiseAbs().maxCoeff());
      max_error_pos = std::max(max_error_pos, diff.bottomRows(2).cwiseAbs().maxCoeff());
    }
  }
  // The yaw jacobian is first order in the relative rotation (here up to 0.035 rad), the position one is exact
  check("yaw jacobians of the predicted measurement", max_error_yaw, 0.03);
  check("position jacobians of the predicted measurement", max_error_pos, 1e-6);

  std::cout << (failed == 0 ? "all wheel tests passed" : "some wheel tests FAILED") << std::endl;
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
