#include "tools.h"
#include <iostream>

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
   
	VectorXd rmse(4);
	 rmse << 0,0,0,0;

	  // check the validity of the following inputs:
	  //  * the estimation vector size should not be zero
	  //  * the estimation vector size should equal ground truth vector size
	 if (estimations.size() != ground_truth.size() || estimations.size() == 0) 
		  
		{
			cout << "Invalid estimation or ground_truth data" << endl;
			return rmse;
		}

	  // accumulate squared residuals
	  for (unsigned int i=0; i < estimations.size(); ++i) 
		{

			VectorXd residual = estimations[i] - ground_truth[i];

			// coefficient-wise multiplication
			residual = residual.array()*residual.array();
			rmse += residual;
		}

	  // calculate the mean
	  rmse = rmse/estimations.size();

	  // calculate the squared root
	  rmse = rmse.array().sqrt();
	  
	  // return the result
	  return rmse;
	
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
   
	   MatrixXd Hj(3,4);
	  // recover state parameters
	  float px = x_state(0);
	  float py = x_state(1);
	  float vx = x_state(2);
	  float vy = x_state(3);

	  // pre-compute a set of terms used in denominator to avoid repeated calculation
	  float d1 = px*px+py*py;
	  float d2 = sqrt(d1);

	  // check division by zero
	  if (fabs(d1) < 0.0001)
	  {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	  }

	  // compute the Jacobian matrix
	  Hj << (px/d2), (py/d2), 0, 0,
		  	-(py/d1), (px/d1), 0, 0,
		  	py*(vx*py - vy*px)/(d1*d2), px*(px*vy - py*vx)/(d1*d2), px/d2, py/d2;

	  return Hj;
}
