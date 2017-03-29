#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transistion matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_laser_;

  // measurement matrix
  Eigen::MatrixXd Hj_radar_;

  // measurement covariance matrix
  Eigen::MatrixXd R_laser_;

  // measurement covariance matrix
  Eigen::MatrixXd R_radar_;

  // state identity matrix
  Eigen::MatrixXd I_;


    /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  /**
   * Calculate the current state in polar coordinates
   * @return Hx The current state in polar coordinates
   */
  Eigen::VectorXd GetPolarState() const;
};

#endif /* KALMAN_FILTER_H_ */
