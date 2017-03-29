#include "kalman_filter.h"

using Eigen::MatrixXf;
using Eigen::VectorXf;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXf Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXf &z) {
  VectorXf z_pred = H_laser_ * x_;
  VectorXf y = z - z_pred;
  MatrixXf Ht = H_laser_.transpose();
  MatrixXf S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXf Si = S.inverse();
  MatrixXf PHt = P_ * Ht;
  MatrixXf K = PHt * Si;

  //new state estimate
  x_ = x_ + (K * y);
  //new probability distribution
  P_ = (I_ - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXf &z) {
  VectorXf Hx = GetPolarState();

  //Calculate
  VectorXf y = z - Hx;
  MatrixXf Ht = Hj_radar_.transpose();
  MatrixXf S = Hj_radar_ * P_ * Ht + R_radar_;
  MatrixXf Si = S.inverse();
  MatrixXf PHt = P_ * Ht;
  MatrixXf K = PHt * Si;

  //new state estimate
  x_ = x_ + (K * y);
  //new state probability distribution
  P_ = (I_ - K * Hj_radar_) * P_;
}

VectorXf KalmanFilter::GetPolarState() const {//Convert current state from cartesian
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  float range = (float) sqrt(pow(px, 2) + pow(py, 2));
  float angle = (float) atan(py / px);
  float velocity = (px*vx + py*vy) / range;

  //Create current polar state
  VectorXf Hx;
  Hx = VectorXf(3);
  Hx << range, angle, velocity;
  return Hx;
}
