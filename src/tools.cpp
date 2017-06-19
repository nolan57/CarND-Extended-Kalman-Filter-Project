#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
MatrixXd Tools::CalculateQ(const float ax2, const float ay2, const float deltaT) {
  /**
   * TODO
   * Calculate process noise Q here
   */
    MatrixXd Q = MatrixXd(4 ,4);
    double deltaT2 = pow(deltaT, 2);
    double deltaT4 = deltaT2*deltaT2;
    double deltaT3 = deltaT2*deltaT;

    Q << deltaT4*ax2/4, 0, deltaT3*ax2/2, 0,
         0, deltaT4*ay2/4, 0, deltaT3*ay2/2,
         deltaT3*ax2/2, 0, deltaT2*ax2, 0,
         0, deltaT3*ay2/2, 0, deltaT2*ay2;

    return Q;
}
