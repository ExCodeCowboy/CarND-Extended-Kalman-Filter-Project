#include <iostream>
#include "tools.h"

using Eigen::VectorXf;
using Eigen::MatrixXf;
using std::vector;


Tools::Tools() {}

Tools::~Tools() {}

VectorXf Tools::CalculateRMSE(const vector<VectorXf> &estimations,
                              const vector<VectorXf> &ground_truth) {
  VectorXf rmse(4);
  rmse << 0, 0, 0, 0;

  unsigned long eSize = estimations.size();
  unsigned long gSize = ground_truth.size();
  if (eSize == 0 || gSize != eSize) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
  }

  for (unsigned int i = 0; i < eSize; i++) {
    VectorXf residual = estimations[i] - ground_truth[i];
    residual = residual.array().pow(2);
    rmse += residual;
  }

  rmse = rmse / eSize;
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXf Tools::CalculateJacobian(const VectorXf &x_state) {
  //Recover state parameters
  float px = (float) x_state(0);
  float py = (float) x_state(1);
  float vx = (float) x_state(2);
  float vy = (float) x_state(3);

  float px_py = (float) (pow(px, 2) + pow(py, 2));

  MatrixXf Hj(3, 4);

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
  float root_px_py = (float) sqrt(px_py);
  float px_py_1_5 = (float) pow(px_py, 1.5);
  float pyv = py * ((vx * py) - (vy * px));
  float pxv = px * ((vy * px) - (vx * py));

  //Assemble the Jacobian matrix
  Hj << px / root_px_py, py / root_px_py, 0, 0,
        -py / px_py, px / px_py, 0, 0,
        pyv / px_py_1_5, pxv / px_py_1_5, px / root_px_py, py / root_px_py;

  return Hj;
}
