#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  noise_ax = 9;
  noise_ay = 9;

  // initializing kalman matrices
  ekf_.R_laser_ = MatrixXd(2, 2);
  ekf_.R_radar_ = MatrixXd(3, 3);
  ekf_.H_laser_ = MatrixXd(2, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.x_ = VectorXd(4);
  ekf_.I_ = MatrixXd::Identity(4, 4);

  /******************************
  /Setup static kalman matrices
  /*****************************/
  //measurement covariance matrix - laser
  ekf_.R_laser_ <<  0.0225, 0,
                    0, 0.0225;

  //measurement covariance matrix - radar
  ekf_.R_radar_ <<  0.09, 0, 0,
                    0, 0.0009, 0,
                    0, 0, 0.09;

  //measurement matrix - laser
  ekf_.H_laser_ <<  1, 0, 0, 0,
                    0, 1, 0, 0;

  /******************************
  /Setup initial kalman state
  /*****************************/
  //kalmen filter initial state transition matrix (without time component)
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  //kalmen filter initial state covariance
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  //kalman filter initial state
  ekf_.x_ << 0, 0, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float v = measurement_pack.raw_measurements_(2);
      float px = ro * cos(phi);
      float py = ro * sin(phi);
      float vx = v * cos(phi);
      float vy = v * sin(phi);
      ekf_.x_ << px, py, vx, vy;

      return;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time component is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //2. Set the process covariance matrix Q
  float t2 = dt * dt;
  float t3 = t2 * dt;
  float t4 = t3 * dt;
  float t4_4 = t4 / 4.0;
  float t3_2 = t3 / 2.0;

  ekf_.Q_ <<  t4_4 * noise_ax, 0, t3_2 * noise_ax, 0,
              0, t4_4 * noise_ay, 0, t3_2 * noise_ay,
              t3_2 * noise_ax, 0, t2 * noise_ax, 0,
              0, t3_2 * noise_ay, 0, t2 * noise_ay;

  //3. Predict new position and probability distribution.
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //skip zero/initialization measurements
    if (measurement_pack.raw_measurements_[0] == 0) {
      return;
    }
    // Generate new jacobian matrix
    ekf_.Hj_radar_ = tools.CalculateJacobian(ekf_.x_);

    // Process radar update
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Process laser update
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
