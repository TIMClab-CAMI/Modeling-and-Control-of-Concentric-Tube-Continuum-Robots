#ifndef SOLVE_IVP_H
#define SOLVE_IVP_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Solve the IVP (Forward integration of the ODEs on each segment, considering an initial guess yu0)
 * 
 * @param yu0 Guess of unknown initial conditions : [mx0, my0, u1z0, u2z0, u3z0] (sum of bending moments along x and y axis, and torsion in each tube, at arc length s=0)
 * @param[in] q Actuation variables [beta_i, alpha_i] with beta the translations (m) and alpha the rotations (rad) 
 * @param[in] segmented Data of the segmented CTR
 * @param[in] f External punctual force applied at the end-effector (3d-vector expressed in the robot base fraome)
 * 
 * @param[out] yTot_out Output matrix of state variables at each integration node
 */
void solveIVP(
  const Eigen::Vector<double,CTR_CONST::NB_YU0>                &yu0,
  const Eigen::Vector<double,CTR_CONST::NB_Q>                  &q,
  const segmentedData                                          &segmented,
  const Eigen::Vector3d                                        &f,

  Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> &yTot_out);

#endif //SOLVE_IVP_H
