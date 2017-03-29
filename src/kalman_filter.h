#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXf x_;

  // state covariance matrix
  Eigen::MatrixXf P_;

  // state transistion matrix
  Eigen::MatrixXf F_;

  // process covariance matrix
  Eigen::MatrixXf Q_;

  // measurement matrix
  Eigen::MatrixXf H_laser_;

  // measurement matrix
  Eigen::MatrixXf Hj_radar_;

  // measurement covariance matrix
  Eigen::MatrixXf R_laser_;

  // measurement covariance matrix
  Eigen::MatrixXf R_radar_;

  // state identity matrix
  Eigen::MatrixXf I_;


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
  void Update(const Eigen::VectorXf &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXf &z);

  /**
   * Calculate the current state in polar coordinates
   * @return Hx The current state in polar coordinates
   */
  Eigen::VectorXf GetPolarState() const;
};

#endif /* KALMAN_FILTER_H_ */
