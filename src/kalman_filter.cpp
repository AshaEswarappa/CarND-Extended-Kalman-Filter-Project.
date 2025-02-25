#include "kalman_filter.h"

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

void KalmanFilter::Predict() 
{
  /**
   * prediction of state
   */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) 
{
  /**
   * KF Measurement update step
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  UpdateNewState(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  float px = x_(0); 
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho     = sqrt(px*px + py*py); // sqrt(px^2+py^2)
  float theta   = atan2(py, px); 
  
  float rho_dot;

 // checks for the negative number to avoid complex number 
  if (fabs(rho) < 0.0001) 
  {
    rho_dot = 0; // assign 0 if its negative
  } 
  else 
  {
    rho_dot = (px*vx + py*vy) / rho;  // (px*vx/py*vy/sqrt(px^2+py^2))
  }

  VectorXd h = VectorXd(3);

  h << rho, theta, rho_dot;
  
  VectorXd y = z-h;

  //Normalize the angle between -pi to pi
  while (y[1] < -M_PI)
  {
     y[1] += 2 * M_PI;
    
  }   
  while (y[1] > M_PI)
  {
    y[1] -= 2 * M_PI;
    
  }
    
  UpdateNewState(y);

}

void KalmanFilter::UpdateNewState(const VectorXd &y)
{

  MatrixXd Ht = H_.transpose(); // Transpose of measurement matrix
  MatrixXd S = H_ * P_ * Ht + R_; 
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}