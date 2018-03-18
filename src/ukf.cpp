#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  //Measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Initialization Boolean
  is_initialized_ = false;
  
  // No. of states
  n_x_ = 5;

  // No. of augmented states
  n_aug_ = 7;
  lambda_ = 3 - n_aug_ ;

  // Initalizing Measurement Noise Covariance for lidar as well as radar
  R_radar_ = MatrixXd(3,3);
  R_lidar_ = MatrixXd(2,2);
  R_radar_ << std_radr_*std_radr_,0,0,
		  0,std_radphi_*std_radphi_,0,
		  0,0,std_radrd_*std_radrd_;

  R_lidar_ << std_laspx_*std_laspx_,0,
		  0,std_laspy_*std_laspy_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // initial predicted sigma matrix
  Xsig_pred_ =  MatrixXd(n_x_,2*n_aug_ +1);

  //Defining weights size
  weights_ = VectorXd(2*n_aug_ +1);
  time_us_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	if (is_initialized_ == false)
	{

	    time_us_ =meas_package.timestamp_;
	    // Intilializing Process Covariance Matrix
	    P_ << 1,0,0,0,0,
	    		0,1,0,0,0,
				0,0,10,0,0,
				0,0,0,1,0,
				0,0,0,0,4;
	    // Defining  weights
		weights_(0) = lambda_/(lambda_ +n_aug_);
		for (int i =1; i < 2*n_aug_ + 1; i++)
		{
			weights_(i) = 0.5/(n_aug_+lambda_);

		}
		//Initializing state depending on first measurement
		 if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		 {
		      /**
		      Converting radar from polar to cartesian coordinates and initializing state.
		      */

		    	float rho, phi,rho_dot;
		    	rho = meas_package.raw_measurements_[0];
		    	phi = meas_package.raw_measurements_[1];
		    	rho_dot= meas_package.raw_measurements_[2];
		    	float px, py, v,yaw,yawd;
		    	px = rho*cos(phi);
		    	py = rho*sin(phi);
		    	v =0;
		    	yaw =0;
		    	yawd =0;
		    	x_ << px,py,v,yaw,yawd;


		 }
		 else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		 {
			 x_<< meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;

		 }
		   // done initializing
		   is_initialized_ = true;

		   return;

	}
	  // time calculation between current and previous measurements
	  double delta_t = (meas_package.timestamp_ -time_us_) / 1000000.0;

	  time_us_ = meas_package.timestamp_;

	  // Calling Predict function
	  Prediction(delta_t);

	  // Measurement Update

	  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)&&(use_radar_ ==true))
	  {
		  UpdateRadar(meas_package);

	  }

	  else if((meas_package.sensor_type_ == MeasurementPackage::LASER)&&(use_laser_==true))
	  {

		  UpdateLidar(meas_package);

	  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
	int col_aug = 2*n_aug_ +1 ;
	MatrixXd Xsig_aug = MatrixXd(n_aug_,col_aug);

	// Augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	// Augmented Covariance matrix
	P_aug.fill(0);
	P_aug.topLeftCorner(n_x_,n_x_) = P_;
	double var_a = std_a_ *std_a_ ;
	double var_yawdd =std_yawdd_ *std_yawdd_;

	P_aug(5,5) = var_a;
	P_aug(6,6) = var_yawdd;

	// Square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//Augmented Sigma points
	Xsig_aug.col(0) = x_aug ;
	for (int i =0;i < n_aug_;i++)
	{
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ +n_aug_)*L.col(i) ;
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ +n_aug_)*L.col(i);

	}

	x_.fill(0);
	//Calculating state sigma prediction matrix
	for(int i =0;i< col_aug;i++)
	{
		double p_x,p_y,v,yaw,yawd,a,yawdd;

		p_x = Xsig_aug(0,i);
		p_y = Xsig_aug(1,i);
		v = Xsig_aug(2,i);
		yaw = Xsig_aug(3,i);
		yawd = Xsig_aug(4,i);
		a = Xsig_aug(5,i);
		yawdd = Xsig_aug(6,i);

		// predicted state values
		double px_p,py_p,v_p,yaw_p,yawd_p;
		double dt2 = delta_t*delta_t;

		if (fabs(yawd)<0.001)
		{
			px_p = p_x + v*delta_t*cos(yaw) +0.5*a*dt2*cos(yaw);
			py_p = p_y +v*delta_t*sin(yaw) + 0.5*a*dt2*sin(yaw);

		}
		else
		{
			px_p = p_x +v/yawd*(sin(yaw +yawd*delta_t)-sin(yaw)) + 0.5*a*dt2*cos(yaw);
			py_p = p_y +v/yawd *(cos(yaw)- cos(yaw+yawd*delta_t)) + 0.5*a*dt2*sin(yaw);
		}
		v_p = v + a*delta_t;
		yaw_p = yaw + yawd*delta_t + 0.5*yawdd*dt2;
		yawd_p = yawd + yawdd*delta_t;

		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;

		// State Mean Prediction

		x_ = x_ + weights_(i)*Xsig_pred_.col(i);

	}
	// State Covariance Matrix Prediction
	P_.fill(0);
	for (int i =0;i <col_aug;i++)
	{
		VectorXd x_del = Xsig_pred_.col(i) - x_;
		// Normalization
		x_del(3) = atan2(sin(double(x_del(3))),cos(double(x_del(3))));
		P_ = P_ + weights_(i)*x_del *x_del.transpose();
	}


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	VectorXd z  = VectorXd(2);
	z<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
	MatrixXd H_lidar = MatrixXd(2,5);
	H_lidar << 1,0,0,0,0,
			0,1,0,0,0;
	VectorXd z_pred =H_lidar*x_;
	VectorXd y = z - z_pred ;
	MatrixXd H_t = H_lidar.transpose();
	MatrixXd PH_t = P_*H_t;
	MatrixXd S = H_lidar*PH_t + R_lidar_ ;
	MatrixXd S_inv = S.inverse();


	// Kalman Filter Gain
	MatrixXd K = PH_t*S_inv;

	// New Estimation of State and State Covariance Matrix
	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I - K*H_lidar)*P_;

	// Calcuation of Lidar NIS.It checks consistency of filter
	double lidar_NIS;
	lidar_NIS = y.transpose()*S_inv*y;
	cout<<"lidar_NIS="<<lidar_NIS<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	int n_z =3;
	int col  = 2*n_aug_ +1 ;

	MatrixXd Zsig = MatrixXd(n_z,col);
	// Calculating Measurement Sigma matrix
	for(int i =0; i< col;i++)
	{
		double p_x,p_y,v,yaw;
		p_x = Xsig_pred_(0,i);
		p_y = Xsig_pred_(1,i);
		v = Xsig_pred_(2,i);
		yaw = Xsig_pred_(3,i);

		double v_x ,v_y;
		v_x = v*cos(yaw);
		v_y = v*sin(yaw);

		double root;
		root = sqrt(p_x*p_x + p_y*p_y);

		// divide by zero safety
		if (fabs(root) < 0.001)
			{
			root = 1;
			cout << "Division by zero";
			}
		// Measurement Model
		double rho,phi,rho_dot;
		rho = root;
		phi = atan2(p_y,p_x);
		rho_dot = (p_x*v_x  + p_y*v_y)/root;
		Zsig(0,i) = rho;
		Zsig(1,i) = phi;
		Zsig(2,i) = rho_dot;
	}

	// Mean Calculation for measurement model
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0);
	for (int i =0;i<col;i++)
	{
		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}
	// innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0);

	// Cross Correlation Matrix
	MatrixXd Tc = MatrixXd(n_x_,n_z);
	Tc.fill(0);
	//Calculating innovation covariance and Cross corelation matrix
	for (int i =0;i <col;i++)
	{
		VectorXd z_del = Zsig.col(i) - z_pred;
		// Normalization
		z_del(1) = atan2(sin(double(z_del(1))),cos(double(z_del(1))));
		S = S + weights_(i)*z_del *z_del.transpose();

		VectorXd x_del = Xsig_pred_.col(i) - x_;
		// Normalization
		x_del(3) = atan2(sin(double(x_del(3))),cos(double(x_del(3))));

		Tc = Tc + weights_(i)*x_del*z_del.transpose();
	}
	// Adding Measurement noise
	S = S + R_radar_;

	// Measurement from sensor
	VectorXd z= VectorXd(n_z);
	z = meas_package.raw_measurements_;

	MatrixXd S_inv = S.inverse();
	//Kalman gain
	MatrixXd K = Tc*S_inv;
	VectorXd zdiff = z-z_pred;
	// Update State Mean  and CoVariance Matrix
	x_ = x_ + K*zdiff;
	P_ = P_ - K*S*K.transpose();

	//Calculating Radar NIS,Checks for consistency of filter
	double radar_NIS;
	radar_NIS = zdiff.transpose()*S_inv*zdiff;
	cout<<"radar_NIS="<<radar_NIS<<endl;


}
