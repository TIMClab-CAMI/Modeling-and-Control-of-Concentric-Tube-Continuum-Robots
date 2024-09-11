#include <iostream>
#include "Eigen/Dense"
#define EIGEN_CORE_H
#include <gnuplot-iostream.h>
#include "ctrConstants.h"

using namespace Eigen;
using namespace CTR_CONST;

void plotCtr(const Matrix<double,nStateVar,nSegMax*nIntPoints>& y, const Vector<int,3>& iEnd, const MatrixXd& xScatter){
  static Gnuplot gp;
  static bool init = false;
  if(!init){
    init = true;
    gp << "set terminal qt" << std::endl;
    gp << "set xrange [-0.1:0.1]" << std::endl;
    gp << "set yrange [-0.1:0.1]" << std::endl;
    gp << "set zrange [0:0.25]" << std::endl;
    gp << "set xlabel 'x (m)'" << std::endl;
    gp << "set ylabel 'y (m)'" << std::endl;
    gp << "set zlabel 'z (m)'" << std::endl;
    gp << "set view equal xyz" << std::endl;
    gp << "set grid xtics ytics ztics" << std::endl;
    gp << "set ticslevel 0" << std::endl;
  }
  int sz_end0 = getYtotIndexFromIend(iEnd(0)) + 1;
  int sz_end1 = getYtotIndexFromIend(iEnd(1)) + 1;
  int sz_end2 = getYtotIndexFromIend(iEnd(2)) + 1;

  MatrixXd r1bis = y.topLeftCorner(n,sz_end0);
  MatrixXd r2bis = y.topLeftCorner(n,sz_end1);
  MatrixXd r3bis = y.topLeftCorner(n,sz_end2);

  MatrixXd r1T = r1bis.transpose();
  MatrixXd r2T = r2bis.transpose();
  MatrixXd r3T = r3bis.transpose();

  auto plots = gp.splotGroup();
  plots.add_plot1d(r1T, "with lines lw 3 linecolor rgb 'red' title 'tube 1'");
  plots.add_plot1d(r2T, "with lines lw 4 linecolor rgb 'blue' title 'tube 2'");
  plots.add_plot1d(r3T, "with lines lw 5 linecolor rgb 'green' title 'tube 3'");
  plots.add_plot1d(xScatter, "with lines lw 1 linecolor rgb 'black' title 'trajectory'");
  gp << plots;
}