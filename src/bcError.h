#ifndef BC_ERROR_H
#define BC_ERROR_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Error function computing the residuals of the boundary conditions
 * 
 * @param y     Matrix containing the state variable at each integration node
 * @param iEnd  Vector containing the index of the segment where each tube ends
 * @param kxy1  Bending stiffness of inner tube
 * @param Ux1   Precurvature of inner tube
 * @return Vector containing the residuals of the boundary conditions [u1z, u2z, u3z, mx, my].
 */

Eigen::Vector<double, CTR_CONST::NB_BC> bcError(
  const Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> &y,
  const Eigen::Vector<int,CTR_CONST::n> &iEnd,
  double kxy1,
  double Ux1);
  
#endif //BC_ERROR_H
