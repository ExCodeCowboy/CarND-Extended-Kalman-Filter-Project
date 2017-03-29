#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new state estimate
  x_ = x_ + (K * y);
  //new probability distribution
  P_ = (I_ - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd Hx = GetPolarState();

  //Calculate
  VectorXd y = z - Hx;
  MatrixXd Ht = Hj_radar_.transpose();
  MatrixXd S = Hj_radar_ * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new state estimate
  x_ = x_ + (K * y);
  //new state probability distribution
  P_ = (I_ - K * Hj_radar_) * P_;
}

VectorXd KalmanFilter::GetPolarState() const {//Convert current state from cartesian
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  float range = sqrt(pow(px, 2) + pow(py, 2));
  float angle = atan(py/px);
  float velocity = (px*vx + py*vy) / range;

  //Create current polar state
  VectorXd Hx = VectorXd(3);
  Hx << range, angle, velocity;
  return Hx;
}
