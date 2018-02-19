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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  // measurement matrix Radar: calculate for every step the jacobian Hj_

  // acceleration noise
  float noise_ax = 6;//9;
  float noise_ay = 6;// 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
  *  Ellapsed time and rest
  ****************************************************************************/
  if (!is_initialized_) // only the very first measurement
	previous_timestamp_ = measurement_pack.timestamp_;
  //compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

  //Check if the delta time is > 3 secs and reset
  // This is needed in case of reset within the simulator
  if ((dt > 3) || (dt < 0)) {
		is_initialized_ = false;
		cout << "Re-initialize" << endl;
  }
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0.25, 0.25;
	double p_x;
	double p_y;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  // Polar coordinates
	  double rho = measurement_pack.raw_measurements_[0];
	  double theta = measurement_pack.raw_measurements_[1];
	  double vrho = measurement_pack.raw_measurements_[2];
	  // Cartesian coordinates
	  p_x = rho * cos(theta);
	  p_y = rho * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
	  p_x = measurement_pack.raw_measurements_[0];
	  p_y = measurement_pack.raw_measurements_[1];
    }
	// Initialize state
	ekf_.x_(0) = p_x;
	ekf_.x_(1) = p_y;

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;

	// state covariance matrix (from lesson)
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	// process noise covariance matrix Q_
	//ekf_.Q_

    // done initializing, no need to predict or update
	//previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   //Update the state transition matrix F according to the elapsed time
   ekf_.F_(0, 2) = dt;
   ekf_.F_(1, 3) = dt;

  // Update the process noise covariance matrix Q
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
	  0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
	  dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
	  0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  cout << "skip radar" << endl;
	  //ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  //ekf_.R_ = R_radar_;
	  //ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
	
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
