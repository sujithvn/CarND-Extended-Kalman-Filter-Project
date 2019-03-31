#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

VectorXd c2p(const VectorXd &x_state) {
	/*
	 * Converts radar measurement from px,py,vx,vy (Cartesian) to rho,phi,rho-dot (polar)
	 */
	float px = x_state[0];
	float py = x_state[1];
	float vx = x_state[2];
	float vy = x_state[3];
	float rho, phi, rho_dot;
	VectorXd z_pred = VectorXd(3);
  z_pred << 0,0,0;

  if(px == 0 && py == 0) {
    return z_pred;
  }

	rho = sqrt((px * px) + (py * py));
	phi = atan2(py, px);
	// When rho is very small, reset to 0.0001
  // This avoids division by 0 error while computing rho_dot
	if (rho < 0.000001)
		rho = 0.000001;
  rho_dot = ((px * vx) + (py * vy)) / rho;

	z_pred << rho, phi, rho_dot;
	return z_pred;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
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
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
	VectorXd z_pred = c2p(x_);
	VectorXd y = z - z_pred;

	// normalize ϕ in the y vector so that it is an angle between −pi and pi
	while (y(1) > M_PI) {
		y(1) -= M_PI;
	}

	while (y(1) < -M_PI) {
		y(1) += M_PI;
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
