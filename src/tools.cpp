#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  unsigned long eSize = estimations.size();
  unsigned long gSize = ground_truth.size();
  if (eSize == 0 || gSize != eSize) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
  }

  for (unsigned int i = 0; i < eSize; i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array().pow(2);
    rmse += residual;
  }

  rmse = rmse / eSize;
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  //Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float px_py = pow(px, 2) + pow(py, 2);

  MatrixXd Hj(3, 4);

  //Check division by zero
  if (px_py == 0) {

    //Input placeholder Hj
    Hj << 0, 0, 0, 0,
        1e+9, 1e+9, 0, 0,
        0, 0, 0, 0;

    //Output a warning
    std::cout << "CalculateJacobian() - error - Division by Zero" << std::endl;
    return Hj;
  }

  //Calculate our component terms
  float root_px_py = sqrt(px_py);
  float px_py_1_5 = pow(px_py, 1.5);
  float pyv = py * ((vx * py) - (vy * px));
  float pxv = px * ((vy * px) - (vx * py));

  //Assemble the Jacobian matrix
  Hj << px / root_px_py, py / root_px_py, 0, 0,
        -py / px_py, px / px_py, 0, 0,
        pyv / px_py_1_5, pxv / px_py_1_5, px / root_px_py, py / root_px_py;

  return Hj;
}
