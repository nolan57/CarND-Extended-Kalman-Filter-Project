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

  // initializing matrices
  //R_laser_ = MatrixXd(2, 2);
  MatrixXd r_laser_(2, 2);
  //R_radar_ = MatrixXd(3, 3);
  MatrixXd r_radar_(3, 3);
  //H_laser_ = MatrixXd(2, 4);
  MatrixXd h_laser_(2, 4);
  //Hj_ = MatrixXd(3, 4);
  MatrixXd hj_(3, 4);

  //measurement covariance matrix - laser
  r_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  r_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  h_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

  R_laser_ = r_laser_;
  R_radar_ = r_radar_;
  H_laser_ = h_laser_;
  Hj_      = hj_;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  VectorXd x_s(4);
  //Tools tools;

  if (is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro;
      double theta;
      double ro_dot;

      ro = measurement_pack.raw_measurements_(0);
      theta = measurement_pack.raw_measurements_(1);
      ro_dot = measurement_pack.raw_measurements_(2);

      double px = cos(theta) * ro;
      double py = - sin(theta) * ro;
      double vx = cos(theta) * ro_dot;
      double vy = - sin(theta) * ro_dot;
      x_s << px, py, vx, vy;

      ekf_.H_ = tools.CalculateJacobian(x_s);
      ekf_.R_ = R_radar_;

      //ekf_.x_ = x_s;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double px = measurement_pack.raw_measurements_(0);
      double py = measurement_pack.raw_measurements_(1);

      double dpx = px - ekf_.x_(0);
      double dpy = py - ekf_.x_(1);
      double vx = dpx / dt;
      double vy = dpy / dt;
      x_s << px, py, vx, vy;

      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      //ekf_.x_ = x_s;

    }}else{

    // done initializing, no need to predict or update
    VectorXd x_in(4);
    MatrixXd P_in(4 ,4);
    MatrixXd F_in(4, 4);
    MatrixXd Q_in(4, 4);

    x_in << 0, 0, 0, 0;
    Q_in << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;

    P_in << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    F_in << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      ekf_.Init(x_in, P_in, F_in, Hj_, R_radar_, Q_in);
      cout << "initialized_ for RADAR" << endl;
    } else {
      // Laser updates
      ekf_.Init(x_in, P_in, F_in, H_laser_, R_laser_, Q_in);
      cout << "initialized_ for LIDAR" << endl;
    }

    is_initialized_ = true;

    cout << "initialized_" << endl;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  ekf_.Q_ = tools.CalculateQ(9.0, 9.0, dt);
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(x_s);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
