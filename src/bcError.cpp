#include <iostream>
#include "Eigen/Dense"
#include "ctrConstants.h"

using namespace Eigen;
using namespace CTR_CONST;

Vector<double, NB_BC> bcError( 
  const Matrix<double,nStateVar,nSegMax*nIntPoints> &y,
  const Vector<int,3> &iEnd,
  double kxy1,
  double Ux1){

  // err(0) : u1z should be zero at the end of first tube
  // err(1) : u2z should be zero at the end of second tube
  // err(2) : u3z should be zero at the end of third tube
  // err(3) : sum of bending moment along x axis should be zero at the end of the CTR
  // err(4) : sum of bending moment along y axis should be zero at the end of the CTR

  Vector<double, NB_BC> err;
  // Convert index of segment to column index in y matrix
  int i_end0 = getYtotIndexFromIend(iEnd(0));
  int i_end1 = getYtotIndexFromIend(iEnd(1));
  int i_end2 = getYtotIndexFromIend(iEnd(2));
  err(0) = y(14,i_end0);
  err(1) = y(15,i_end1);
  err(2) = y(16,i_end2);

  // Compute sum of bending moments along x and y axis
  double R11 = y(3,i_end0);
  double R12 = y(4,i_end0);
  double R21 = y(6,i_end0);
  double R22 = y(7,i_end0);

  double u1x = y(12,i_end0);
  double u1y = y(13,i_end0);

  double delta_u1x = u1x - Ux1;
  
  err(3) = R11 * kxy1 * delta_u1x + R12 * kxy1 * u1y;
  err(4) = R21 * kxy1 * delta_u1x + R22 * kxy1 * u1y;
  return err;
}
