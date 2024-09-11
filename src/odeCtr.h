#ifndef ODE_CTR_H
#define ODE_CTR_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Ordinary differential equations for CTR
 * 
 * @param[in] s arc-length
 * @param[in] y state variables
 * @param[in] Kxy Vector of bending stiffness for each tube
 * @param[in] Kz Vector of torsional stiffness for each tube
 * @param[in] Ux Vector of precurvature for each tube
 * @param[in] f External punctual force applied at the end-effector
 *
 * @param[out] y_s_out Derivative of the state variables
 */

void odeCtr(
  double s,
  Eigen::Vector<double,CTR_CONST::nStateVar>& y,
  const Eigen::Vector<double,CTR_CONST::n>& Kxy,
  const Eigen::Vector<double,CTR_CONST::n>& Kz,
  const Eigen::Vector<double,CTR_CONST::n>& Ux,
  const Eigen::Vector3d& f,

  Eigen::Vector<double,CTR_CONST::nStateVar>& y_s_out);

/*! Maps a vector in R3 to a 3x3 skew-symettric matrix in so3. */
inline Eigen::Matrix3d hat(Eigen::Vector3d y){
    Eigen::Matrix3d y_hat;
    y_hat <<    0, -y(2),  y(1),
             y(2),     0, -y(0),
            -y(1),  y(0),     0;

    return y_hat;
}

#endif //ODE_CTR_H
