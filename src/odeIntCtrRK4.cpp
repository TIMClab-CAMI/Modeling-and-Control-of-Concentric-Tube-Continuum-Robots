#include "Eigen/Dense"
#include "odeCtr.h"
#include "ctrConstants.h"
#include <iostream>

using namespace Eigen;
using namespace CTR_CONST;

Matrix<double,nStateVar,nIntPoints> odeIntCtrRK4(
  Vector<double,nStateVar> y0,
  double s0,
  double sf,
  const Vector<double,n>& Kxy,
  const Vector<double,n>& Kz,
  const Vector<double,n>& Ux,
  const Vector3d& f){
    
  Matrix<double,nStateVar,nIntPoints> y;
  y.col(0) = y0;

  double ds = (sf-s0)/(nIntPoints-1);
  double half_ds = ds/2;
  double sixth_ds = ds/6;
  double L = sf-s0;
  constexpr int Nm1 = nIntPoints-1;

  //Classic 4th-order Runge-Kutta method
  Vector<double,nStateVar> k0, k1, k2, k3;
  double s = s0;
  for(int i = 0; i < Nm1; i++){
    odeCtr(s,y0, Kxy, Kz, Ux, f, k0);
    y0 += k0*half_ds;
    s += half_ds;
    odeCtr(s,y0, Kxy, Kz, Ux, f, k1);
    y0 = y.col(i) + k1*half_ds;
    odeCtr(s,y0, Kxy, Kz, Ux, f, k2);
    s = s0 + (L*(i+1))/Nm1;
    y0 = y.col(i) + k2*ds;
    odeCtr(s,y0, Kxy, Kz, Ux, f, k3);

    y0 = y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
    y.col(i+1) = y0;
  }
  return y;
}
