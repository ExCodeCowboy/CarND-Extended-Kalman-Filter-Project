#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXf CalculateRMSE(const std::vector<Eigen::VectorXf> &estimations, const std::vector<Eigen::VectorXf> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXf CalculateJacobian(const Eigen::VectorXf& x_state);

};

#endif /* TOOLS_H_ */
