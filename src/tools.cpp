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
    MatrixXd Hj (3, 4);
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    double c1 = px*px + py*py;
    double c2 = sqrt(c1);
    double c3 = c1*c2;

    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    Hj <<  (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    return  Hj;
}
MatrixXd Tools::CalculateQ(const float ax2, const float ay2, const float deltaT) {
  /**
   * TODO
   * Calculate process noise Q here
   */
    MatrixXd Q(4 ,4);
    double deltaT2 = pow(deltaT, 2);
    double deltaT4 = deltaT2*deltaT2;
    double deltaT3 = deltaT2*deltaT;

    Q << deltaT4*ax2/4, 0, deltaT3*ax2/2, 0,
         0, deltaT4*ay2/4, 0, deltaT3*ay2/2,
         deltaT3*ax2/2, 0, deltaT2*ax2, 0,
         0, deltaT3*ay2/2, 0, deltaT2*ay2;

    return Q;
}
