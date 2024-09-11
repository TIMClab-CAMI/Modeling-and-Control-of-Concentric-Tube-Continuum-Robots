#ifndef ODE_INT_CTR_RK4_H
#define ODE_INT_CTR_RK4_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Forward integration of ODEs using Runge-Kutta 4 scheme
 * 
 * @param[in] y0 value of state variables at the beginning of the segment
 * @param[in] s0 arc-length at the beginning of the segment
 * @param[in] sf arc-length at the end of the segment
 * @param[in] Kxy Vector of bending stiffness for each tube on the considered segment
 * @param[in] Kz Vector of torsional stiffness for each tube on the considered segment
 * @param[in] Ux Vector of precurvature for each tube on the considered segment
 * @param[in] f external punctual force applied at the end-effector
 * 
 * @return Matrix containing the values of the state variables for each integration node on the considered segment 
 */

Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nIntPoints> odeIntCtrRK4(
  Eigen::Vector<double,CTR_CONST::nStateVar> y0,
  double s0,
  double sf,
  const Eigen::Vector<double,CTR_CONST::n>& Kxy,
  const Eigen::Vector<double,CTR_CONST::n>& Kz,
  const Eigen::Vector<double,CTR_CONST::n>& Ux,
  const Eigen::Vector3d& f);

#endif //ODE_INT_CTR_RK4_H
