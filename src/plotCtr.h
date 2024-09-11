#ifndef  PLOT_CTR_H
#define  PLOT_CTR_H

#include "Eigen/Dense"
#include "ctrConstants.h"

void plotCtr(
  const Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints>& y,
  const Eigen::Vector<int,CTR_CONST::n>& iEnd,
  const Eigen::MatrixXd& xScatter);
#endif // PLOT_CTR_H
