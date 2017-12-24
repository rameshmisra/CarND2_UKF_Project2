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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  // Parameters above this line are scaffolding, do not modify
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  weights = VectorXd(2*n_aug_+1);
  //set weights for sigma points used for Radar
  weights(0) = lambda_/(lambda_ + n_aug_);
  otherweights = 1/(2*(lambda_ + n_aug_));
  for(int i=1; i<(2*n_aug_ +1); i++){
      weights(i) = otherweights;
  }


  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1,0,0,0,0,
              0,1,0,0,0;
  

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ *std_radr_ , 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_ *std_radrd_ ;

  P_ = MatrixXd::Identity(n_x_, n_x_);
  Xsig_pred_ = MatrixXd(n_x_, n_x_);


  //is_initialized_ = false;

  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    

    time_us_= meas_package.timestamp_;
    x_.fill(0.0);


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho *cos(theta);
      x_(1) = rho *sin(theta);
      
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_= meas_package.timestamp_;

  
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    Prediction(delta_t);
    UpdateLidar(meas_package);
    //cout<< "update lidar done" <<endl;
  } 
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    Prediction(delta_t);
    UpdateRadar(meas_package);
    //cout<< "update radar done" <<endl;
  }

  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  

  /*MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  MatrixXd P_sqrt = P_.llt().matrixL();

  Xsig.col(0) = x_;
  MatrixXd B = sqrt(lambda_ + n_x_)*P_sqrt;
  MatrixXd C = -B;
  B.colwise() += x_;
  C.colwise() += x_;
  for(int i = 0; i< n_x_; i++){
      Xsig.col(i+1) = B.col(i);
      Xsig.col(i+1+n_x_) = C.col(i);
  }*/

  //Creation of augmented sigma points

  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_*std_yawdd_;
  
  
  //create square root matrix
  MatrixXd P_augsqrt = P_aug.llt().matrixL();
  
  //create augmented sigma points
  
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_; i++){
      
      Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_)*P_augsqrt.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*P_augsqrt.col(i);
  }

  //cout << "UKF Prediction of augmented sigma points done " << endl;

  //Prediction of sigma points for next timestep:

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
  //cout << "predicted  Xsig_aug : " << Xsig_aug(0,0) <<endl;

  //Calculate State Mean & Covariance from above sigma points

  //predict state mean
  x_.fill(0.0);
  P_.fill(0.0);
  x_ = Xsig_pred_ * weights;

  /*for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights(i) * Xsig_pred_.col(i);
  }*/

  //cout << "predicted  X_ : " << x_(0) <<endl;
  //predict state covariance matrix
  
  MatrixXd Diff = Xsig_pred_.colwise() - x_;
  for(int i=0; i<2*n_aug_ +1; i++){
      while (Diff(3,i) > M_PI) Diff(3,i) -= 2.0 *M_PI;
      while (Diff(3,i) < -M_PI) Diff(3,i) += 2.0 *M_PI;

  }
  

  MatrixXd W = MatrixXd(n_x_, 2*n_aug_ +1);
  W = otherweights * Diff;
  
  //Make correction for zeroth column of W
  for(int i=0; i<n_x_; i++){
      W(i,0) = lambda_/(lambda_ + n_aug_) * Diff(i, 0);
  }
  P_ = W * Diff.transpose();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  
  
  
  
  //Actual measurement data
  VectorXd z_ = VectorXd(2);
  z_(0) = meas_package.raw_measurements_(0);
  z_(1) = meas_package.raw_measurements_(1);
  //cout<< z_(0) << z_(1) <<endl;

  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z_ - z_pred;
  
  MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
  
  
  MatrixXd K = P_ * H_laser_.transpose() * S.inverse();
  //cout<< "P_:  "<<P_<<endl;

  x_ = x_ + (K * y);
  //cout<< x_(0) << endl;
  
  MatrixXd I = MatrixXd::Identity(5, 5);
  P_ = (I - K * H_laser_) * P_;

  //NIS calculation

  double NIS_laser = y.transpose() * S.inverse() * y;


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;
  
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  Zsig.fill(0.0);
  for(int l=0; l<(2*n_aug_ + 1); l++){
      Zsig(0, l) = sqrt(Xsig_pred_(0, l)*Xsig_pred_(0, l) + Xsig_pred_(1, l)*Xsig_pred_(1, l));
      Zsig(1, l) = atan2(Xsig_pred_(1, l), Xsig_pred_(0, l));
      Zsig(2, l) = Xsig_pred_(2, l)*(Xsig_pred_(0, l)*cos(Xsig_pred_(3, l)) + Xsig_pred_(1, l)*sin(Xsig_pred_(3, l)))/Zsig(0, l);
  }
  
  //calculate mean predicted measurement
  z_pred = Zsig * weights;
  
  //calculate predicted measurement covariance matrix S
  
  MatrixXd ZDiff = Zsig.colwise() - z_pred;

  for(int i=0; i<2*n_aug_ +1; i++){
      while (ZDiff(1,i) > M_PI) ZDiff(1,i) -= 2.0 *M_PI;
      while (ZDiff(1,i) < -M_PI) ZDiff(1,i) += 2.0 *M_PI;
  }

  MatrixXd wxZDiff = MatrixXd(n_z, 2*n_aug_ + 1);
  wxZDiff = otherweights * ZDiff;
  for(int i=0; i<n_z; i++){
      wxZDiff(i,0) = lambda_/(lambda_ + n_aug_) * ZDiff(i, 0);
  }
  S = wxZDiff * ZDiff.transpose() + R_radar_;


  //Now calculate updated state mean & covariance, using Predicted Measurement, Predicted State, & Actual Measurement

  //Actual measurement data
  VectorXd z_ = VectorXd(3);
  z_(0) = meas_package.raw_measurements_(0);
  z_(1) = meas_package.raw_measurements_(1);
  z_(2) = meas_package.raw_measurements_(2); 
  

  //create matrix for cross correlation Txc
  MatrixXd Txc = MatrixXd(n_x_, 3);
  //calculate cross correlation matrix
  MatrixXd XDiff = Xsig_pred_.colwise() - x_;

  for(int i=0; i<2*n_aug_ +1; i++){
      while (XDiff(3,i) > M_PI) XDiff(3,i) -= 2.0 *M_PI;
      while (XDiff(3,i) < -M_PI) XDiff(3,i) += 2.0 *M_PI;
  }
  
  MatrixXd wxXDiff = MatrixXd(n_x_, 2*n_aug_ +1);
  wxXDiff = otherweights * XDiff;
  for(int i=0; i<n_x_; i++){
      wxXDiff(i,0) = lambda_/(lambda_ + n_aug_) * XDiff(i, 0);
  }
  
  
  Txc = wxXDiff * ZDiff.transpose();
  //std::cout<<"Tc \n"<<Tc<<std::endl;

  //calculate Kalman gain K;
  MatrixXd K = Txc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd error_y = z_ - z_pred;
  x_ += K*error_y;
  //cout<< x_(0) << x_(1) << endl;
  
  P_ -= K*S*K.transpose();


  //NIS calculation

  double NIS_radar = error_y.transpose() * S.inverse() * error_y;




}
