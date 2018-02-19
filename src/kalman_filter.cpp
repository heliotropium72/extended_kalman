#include "kalman_filter.h"
#include "math.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //Predicted state in cartesian coordinates
  float p_x = x_(0);
  float p_y = x_(1);
  float v_x = x_(2);
  float v_y = x_(3);
  // Map the predicted state to radar measurements (polar coordinated)
  VectorXd z_pred(3);
  z_pred(0) = sqrt(p_x*p_x + p_y*p_y);
  z_pred(1) = atan2(p_y, p_x);
  z_pred(2) = (p_x*v_x + p_y*v_y) / z_pred(0);
  // Difference between measurement and prediction
  VectorXd y = z - z_pred;
  // Second entry of y (theta) should be between -Pi and Pi  
  if (y(1) > M_PI || y(1) < -M_PI) {
	while (y(1) > M_PI)
		y(1) -= 2.0 * M_PI;
	while (y(1) < -M_PI)
		y(1) += 2.0 * M_PI;
  }

  // H is the Jacobian of the current predicted state now
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
