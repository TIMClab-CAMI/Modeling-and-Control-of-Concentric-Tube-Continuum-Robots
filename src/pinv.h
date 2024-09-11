#ifndef PINV_H
#define PINV_H

#include "Eigen/Dense"

/**
 * @brief Compute pseudo-inverse of matrix
 */

template<typename _Matrix_Type_>
_Matrix_Type_ pinv(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon()){
  // For a square matrix
  //  (Here we use pinv for square matrices (Bu and temporary 6x6 expression in GDLS))
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
  // For a nonsquare matrix
  //Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}
#endif //PINV_H
