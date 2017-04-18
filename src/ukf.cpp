#include "ukf.h"
#include "tools.h"
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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 0.5, 0, 0, 0, 0,
        0, 0.5, 0, 0, 0,
        0, 0, 0.6, 0, 0,
        0, 0, 0, 0.6, 0,
        0, 0, 0, 0, 0.6;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.85;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.85;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //previous timestamp
  previous_timestamp_ = 0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Lidar measurement dimension
  n_z_L_ = 2;

  // Radar measurement dimension
  n_z_R_ = 3;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;
  // Augmented Sigma point spreading parameter
   aug_lambda_ = 3 - n_aug_;
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  //Lidar measurement covariance matrix
  R_laser_ =  MatrixXd(n_z_L_, n_z_L_);

  R_laser_<< std_laspx_, 0,
             0, std_laspy_;

  R_radar_ = MatrixXd(n_z_R_,n_z_R_);
  R_radar_ <<   std_radr_*std_radr_, 0, 0,
               0, std_radphi_*std_radphi_, 0,
               0, 0,std_radrd_*std_radrd_;

  // the current NIS for radar
  NIS_radar_=0.001;

  // the current NIS for laser
  NIS_laser_=0.002;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /*****************************************************************************
       *  Initialization
       ****************************************************************************/
      if (!is_initialized_) {

        // first measurement
        cout << "UKF: " << endl;


        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          /**
          Convert radar from polar to cartesian coordinates and initialize state.
          */
            float r = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            float ro_dot = meas_package.raw_measurements_[2];
           // cout << "initialing position RADAR Measurement" << endl;

            x_ << r * cos(phi), r * sin(phi), ro_dot,  phi,  0;
            cout << "initial position RADAR Measurement: x=" << x_[0]  << " y=" << x_[1] << endl;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          /**
          Initialize state.
          */
           //cout << "initialing position LIDAR Measurement" << endl;

           // px, py, v, yaw, yaw rate
           x_ <<  meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 1, 0.1, 0;


           cout << "initial position LIDAR Measurement: x=" << x_[0]  << " y=" <<  x_[1] << endl;
        }

       GenerateSigmaPoints();

        // set weights
        double weight_0 =  aug_lambda_/( aug_lambda_+n_aug_);
        weights_(0) = weight_0;
        for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
          double weight = 0.5/(n_aug_+ aug_lambda_);
          weights_(i) = weight;
         }
        // done initializing, no need to predict or update
        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
      }


      //compute the time elapsed between the current and previous measurements
      double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
      cout << dt << endl;
     // if(dt > 0.001){ //if two measurements are are approximately coincident in the same time, just use previous state, no predittion is needed

          previous_timestamp_ = meas_package.timestamp_;

          cout << "Do Prediction" << endl;
          /*****************************************************************************
           *  Prediction
           ****************************************************************************/


          Prediction(dt);
          cout << "Done Prediction" << endl;


          // cout << "x_ = " << ekf_.x_ << endl;
     // }else{
     //     cout << "if two measurements are are approximately coincident in the same time, just use previous state, no predittion is needed" << endl;
    //  }

      /*****************************************************************************
       *  Update
       ****************************************************************************/


      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates

          //cout << "Do Radar udpate" << endl;
          UpdateRadar(meas_package);
          cout << "Done Radar udpate" << endl;


      } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
        // Laser updates
         // cout << "Do Lidar udpate" << endl;
          UpdateLidar(meas_package);
          cout << "Done Lidar udpate" << endl;
      }

      // print the output
      cout << "x_ = " << x_ << endl;
      cout << "P_ = " << P_ << endl;
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

      //Generate augented Sigma points
      MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

      //Craete Augmented Sigma Points
      AugmentedSigmaPoints(&Xsig_aug);

      //Predict Sigma Points
      SigmaPointPrediction(Xsig_aug, delta_t);

      //predicted state mean
      x_.fill(0.0);
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
      }

      //predicted state covariance matrix
      P_.fill(0.0);
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) =   NormalizeYaw(x_diff(3));

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
      }


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
    //Measurement vector
    VectorXd z = meas_package.raw_measurements_;
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_L_, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
         // extract values for better readibility
         double px = Xsig_pred_(0,i);
         double py = Xsig_pred_(1,i);
         // measurement model
         Zsig(0,i) = px;
         Zsig(1,i) = py;
   }

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z_L_);
   z_pred.fill(0.0);
   for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

    //measurement covariance matrix
   MatrixXd S = MatrixXd(n_z_L_,n_z_L_);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred; //residual      
      S = S + weights_(i) * z_diff * z_diff.transpose();
   }


   S = S + R_laser_;

   //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_L_);

    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // state difference
         VectorXd x_diff = Xsig_pred_.col(i) - x_;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
      }
      //calculate Kalman gain K;

      MatrixXd SI =  S.inverse();

      MatrixXd K = Tc * SI;

      //residual
      VectorXd z_diff = z - z_pred;

      //update state mean and covariance matrix
      x_ = x_ + K * z_diff;
     // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
      P_ = P_ - K * S * K.transpose();

      //angle normalization
      x_(3) =   NormalizeYaw(x_(3));

      //Calculate Z(K+1|k)
      VectorXd z_new = VectorXd(2);
      z_new(0) = x_(0);
      z_new(1) = x_(1);

      VectorXd z_diff_new = z - z_new;

      NIS_laser_ = z_diff_new.transpose() * SI * z_diff_new;

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
     //create  vector for incoming lidar measurement
    //Measurement vector
    VectorXd z = meas_package.raw_measurements_;


    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_R_ , 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        // extract values for better readibility
            double p_x = Xsig_pred_(0,i);
            double p_y = Xsig_pred_(1,i);
            double v  = Xsig_pred_(2,i);
            double yaw = Xsig_pred_(3,i);

            double v1 = cos(yaw)*v;
            double v2 = sin(yaw)*v;

            // measurement model
            Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);           //rho
            if(p_x < 0.001 && p_y < 0.001) //both px and py are zero
                Zsig(1,i) = 0;                                //phi
            else
               Zsig(1,i) = atan2(p_y,p_x);              //phi
            if( Zsig(0,i) > 0.001)
                Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
            else
               Zsig(2,i) = 0;

   }

   //mean predicted measurement
   VectorXd z_pred = VectorXd(n_z_R_);
   z_pred.fill(0.0);
   for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
   }

    //measurement covariance matrix
   MatrixXd S = MatrixXd(n_z_R_,n_z_R_);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred; //residual

       //angle normalization
      z_diff(1) =   NormalizeYaw(z_diff(1));
      S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //add measurement noise covariance matrix
   S = S + R_radar_;

   //create matrix for cross correlation Tc
     MatrixXd Tc = MatrixXd(n_x_, n_z_R_);


     //calculate cross correlation matrix

     Tc.fill(0.0);
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

       //residual
       VectorXd z_diff = Zsig.col(i) - z_pred;
       //angle normalization

       z_diff(1) =   NormalizeYaw(z_diff(1));

       // state difference
       VectorXd x_diff = Xsig_pred_.col(i) - x_;
       //angle normalization

       x_diff(3) =   NormalizeYaw(x_diff(3));


       Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
     }
     //calculate Kalman gain K;
     MatrixXd SI =  S.inverse();
     MatrixXd K = Tc * SI;

     //residual
     VectorXd z_diff = z - z_pred;

     //angle normalization
     //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
     z_diff(1) =   NormalizeYaw(z_diff(1));
     //NormalizeAngle(&z_diff(1));



     //update state mean and covariance matrix
     x_ = x_ + K * z_diff;

     //angle normalization

     x_(3) =   NormalizeYaw(x_(3));


    // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
     P_ = P_ - K * S * K.transpose();


     // Caculate Z(k+1|k)
     VectorXd z_new = VectorXd(3);

     double p_x = x_(0);
     double p_y = x_(1);
     double v  = x_(2);
     double yaw = x_(3);

     double v1 = cos(yaw)*v;
     double v2 = sin(yaw)*v;

     // measurement model
     z_new(0) = sqrt(p_x*p_x + p_y*p_y);           //rho
     if(p_x < 0.001 && p_y < 0.001) //both px and py are zero
         z_new(1) = 0;                                //phi
     else
        z_new(1) = atan2(p_y,p_x);              //phi
     if( z_new(0) > 0.001)
        z_new(2) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
     else
       z_new(2)  = 0;

     VectorXd z_diff_new = z - z_new;
     z_diff_new(1) =   NormalizeYaw(z_diff_new(1));

     //Caculate radar NIS
     NIS_radar_ = z_diff_new.transpose() * SI * z_diff_new;

}

/**
 * @brief UKF::GenerateSigmaPoints
 * @
 */
void UKF::GenerateSigmaPoints() {

    //create augmented mean vector
   VectorXd x_aug = VectorXd(n_aug_);

   //create augmented state covariance
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

   //create sigma point matrix
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

 /*******************************************************************************
  * Student part begin
  ******************************************************************************/

   //create augmented mean state
   x_aug.head(5) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;

   //create augmented covariance matrix
   P_aug.fill(0.0);
   P_aug.topLeftCorner(n_x_,n_x_) = P_;
   P_aug(5,5) = std_a_*std_a_;
   P_aug(6,6) = std_yawdd_*std_yawdd_;

   //create square root matrix
   MatrixXd L = P_aug.llt().matrixL();

   //create augmented sigma points
   Xsig_aug.col(0)  = x_aug;
   for (int i = 0; i< n_aug_; i++)
   {
     VectorXd temp = sqrt(aug_lambda_+n_aug_) * L.col(i);
     Xsig_aug.col(i+1)       = x_aug + temp;
     Xsig_aug.col(i+1+n_aug_) = x_aug - temp;
   }

   //write result
   Xsig_pred_ = Xsig_aug.topLeftCorner(n_x_,2*n_aug_ + 1);

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

   //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    VectorXd temp = sqrt(aug_lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1)       = x_aug + temp;
    Xsig_aug.col(i+1+n_aug_) = x_aug - temp;
  }

  //write result
  *Xsig_out = Xsig_aug;

}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  // cout << "Extracted valued form aug sigma points" << endl;
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
   //  cout << "Predicted Simga Points: " << Xsig_aug  << endl;

}

double UKF:: NormalizeYaw(double cal_yaw){
    //cout << "Input yaw: " << cal_yaw << endl;
    /*double nor_yaw = 0.;
    double pi_times = cal_yaw/2./M_PI;
    int int_pi_times = int(pi_times);
    double residual = pi_times - int_pi_times;
    cout << "Residual: " << residual << endl;
    if(cal_yaw > M_PI){

      nor_yaw = (residual -1) * 2. * M_PI;
    }
    if(cal_yaw < -M_PI){
     nor_yaw = (residual + 1) *2. *M_PI;
    }
    cout << "Normalized yaw: " << cal_yaw << endl; */

    double nor_yaw = cal_yaw - round(cal_yaw / (2.0d * M_PI)) * (2.0d * M_PI);
    return nor_yaw;

}


