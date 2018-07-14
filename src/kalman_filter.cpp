#include "kalman_filter.h"
#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_ ;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S = H_ * P_ * Ht_ + R_ ;
  MatrixXd K = P_ * Ht_ * S.inverse();

  x_ = x_ + K * y;
  float x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_ ;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd y;

  float c = x_(0) * x_(0) + x_(1)*x_(1);
  float ro = sqrt(c);
  float ro_dot;
  if (fabs(ro) < 0.0001){
    ro_dot = 0;
  }else{
    ro_dot = (x_(0) * x_(2) + x_(1) * x_(3)) / ro;
  }
  float theta = atan2(x_(1), x_(0) );
  VectorXd zp(3);
  zp << ro, theta, ro_dot;
  y = z - zp;

  MatrixXd Ht_ = H_.transpose();
  MatrixXd S = H_ * P_ * Ht_ + R_ ;
  MatrixXd K = P_ * Ht_ * S.inverse();

  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_ ;
}
