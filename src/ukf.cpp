#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // Set this to false so that the initialization happens for the first measurement.
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.35;

  // State dimension
  n_x_ = 5;

 // Augmented dimension
 n_aug_ = 7;

 // spreading parameter
 lambda_ = 3 - n_aug_;

 NIS_radar_ = 0.0;

 NIS_laser_ = 0.0;

 R_radar_ = MatrixXd(3, 3);

 R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_ * std_radrd_;

 R_laser_ = MatrixXd(2, 2);
 R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;
} 

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
        P_ = MatrixXd::Identity(5, 5);
	P_(0,0) = 0.1;
	P_(1,1) = 0.1;
  	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        /**
      	  Convert radar from polar to cartesian coordinates and initialize state.
         */

      	    double rho = meas_package.raw_measurements_[0];
      	    double phi = meas_package.raw_measurements_[1];
            double rho_dot = meas_package.raw_measurements_[2];

            x_(0) = rho * cos(phi);
            x_(1) = rho * sin(phi);
            x_(2) = sqrt(((rho_dot * cos(phi)) * (rho_dot * cos(phi))) + ((rho_dot * sin(phi)) * (rho_dot * sin(phi))));
            x_(3) = 0;
            x_(4) = 0;

            // Avoid divide by zero
            for (int i = 0; i < 3; i++)
               if (x_(i) < EPSILON)
                     x_(i) = EPSILON;
       }
       else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
         /**
          Initialize state.
         */
         x_(0) = meas_package.raw_measurements_[0];
         x_(1) = meas_package.raw_measurements_[1];
         x_(2) = 0;
         x_(3) = 0;
         x_(4) = 0;

         // Avoid divide by zero
         for (int i = 0; i < 2; i++)
             if (x_(i) < EPSILON)
                   x_(i) = EPSILON;
       }


       // Nothing to predict or update for the first measurement.
       previous_timestamp_ = meas_package.timestamp_;
       is_initialized_ = true;
       //cout << "Initialization completed. " << endl;
       return;
   }

   float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = meas_package.timestamp_;

   // Predict
   Prediction(delta_t);

   // Update
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    UpdateLidar(meas_package);
  }
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  //write result
  *Xsig_out = Xsig_aug;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Step # 1: Predict sigma points

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); 
  AugmentedSigmaPoints(&Xsig_aug);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Step # 2: Predict the state mean

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i < 2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    x_ = x_ + weights(i) * Xsig_pred_.col(i);
  }

  // Step # 3: Predict state covariance
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //cout << "Entered Lidar" << endl;

  //cout << "LIDAR: Entered Lidar Update ... " << endl;
  int n_z = 2;

  VectorXd weights = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i < 2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights(i) = weight;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Zsig(0, i) = Xsig_pred_(0,i);
      Zsig(1, i) = Xsig_pred_(1,i);
  }
  
  //calculate mean predicted measurement 
  for (int i = 0; i < n_z; i++) {
    z_pred(i) = 0;
    for (int j = 0; j < 2 * n_aug_ + 1; j++) {
        z_pred(i) += weights(j) *  Zsig(i, j);
    }
  } 

  //cout << "LIDAR: z_pred = " << z_pred << endl;
  //cout << "LIDAR: Zsig = " << Zsig << endl;

  UpdateUKF(z_pred, Zsig, n_z, meas_package);
  //cout << "LIDAR: After UpdateUKF... " << endl;
}

void UKF::UpdateUKF(VectorXd z_pred, MatrixXd Zsig, int n_z, MeasurementPackage meas_package) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);


  // set weights
  VectorXd weights = VectorXd(2 * n_aug_+1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i < 2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights(i) = weight;
  }


  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //cout << "weights == " << weights << endl;

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //cout << "S matrix == " << S << endl;
  //cout << "S matrix inverse == " << S.inverse() << endl;

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    R = R_laser_;
  }

  S = S + R;

  //cout << "S matrix == " << S << endl;
  
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  //if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  //}

  //update state mean and covariance matrix

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose(); 


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar

        //cout << "S inverse == " << S.inverse() << endl;
        //cout << "z_transpose * S inverse ==" << z.transpose() * S.inverse() << endl;
        //cout << "z.transpose() == " << z.transpose() << endl;
	NIS_radar_ = z.transpose() * S.inverse() * z;
	//cout << "NIS_radar_ = " << NIS_radar_ << endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) { // Lidar
	NIS_laser_ = z.transpose() * S.inverse() * z;
	//cout << "NIS_laser_ = " << NIS_laser_ << endl;
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //cout << "Entered Radar" << endl;

  int n_z = 3;

  //cout << "RADAR: Inside UpdateRadar.." << endl;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  VectorXd weights = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i < 2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights(i) = weight;
  }


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
      double rho = sqrt(p_x*p_x + p_y*p_y);
      double si = atan2(p_y, p_x);
      double rho_dot; 
      if (rho != 0) 
 	rho_dot = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/ rho;
      else
	rho_dot = 0;
      
      Zsig(0, i) = rho;
      Zsig(1, i) = si;
      Zsig(2, i) = rho_dot;
  }

  //calculate mean predicted measurement
  for (int i = 0; i < n_z; i++) {
    z_pred(i) = 0;
    for (int j = 0; j < 2 * n_aug_ + 1; j++) {
        z_pred(i) += weights(j) *  Zsig(i, j);
    }
  }

  UpdateUKF(z_pred, Zsig, n_z, meas_package);
  //cout << "RADAR: After completing the UKF update .." << endl;
}
