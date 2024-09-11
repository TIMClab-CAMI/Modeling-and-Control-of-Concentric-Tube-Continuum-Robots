#ifndef IVP_FD_H
#define IVP_FD_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Compute  jacobian matrices using the IVP finite differences method
 * 
 * @tparam opt Option for chosing which elements to compute
 * @param[in] q Actuation variables [beta_i, alpha_i] with beta the translations [m] and alpha the rotations [rad] 
 * @param[in] yu0 Guess of unknown initial conditions : [mx0, my0, u1z0, u2z0, u3z0] (sum of bending moments along x and y axis, and torsion in each tube, at arc length s=0)
 * @param[in] Kxy Vector of bending stiffness for each tube
 * @param[in] Kz Vector of torsional stiffness for each tube
 * @param[in] Ux Vector of precurvature for each tube
 * @param[in] l Vector of effective length for each tube
 * @param[in] l_k Vector of precurved length for each tube (we assume a tube is composed of a straight part and a precurved part with constant precurvature)
 * @param[in] f External punctual force applied at the end-effector (3d-vector expressed in the robot base fraome)
 * 
 * @param[out] yTot_out Output state variables at each integration node associated with the "nominal" computation solveIVP(yu0,q,...)
 * @param[out] b_out Output residuals of boundary conditions associated with the "nominal" computation solveIVP(yu0,q,...)
 * @param[out] Eq_out Output jacobian matrix Eq (delta_X/delta_q)
 * @param[out] Eu_out Output jacobian matrix Eu (delta_X/delta_yu0)
 * @param[out] Bq_out Output jacobian matrix Bq (delta_b/delta_q)
 * @param[out] Bu_out Output jacobian matrix Bu (delta_b/delta_yu0)
 * @return int (0 indicates success and < 0 indicates failure )
 */
  template <COMPUTATION_OPTION opt> 
  int ComputeIvpJacobianMatrices( 
    const Eigen::Vector<double,CTR_CONST::NB_Q> &q,
    const Eigen::Vector<double,CTR_CONST::NB_YU0> &yu0,
    const Eigen::Vector<double,CTR_CONST::n> &Kxy,
    const Eigen::Vector<double,CTR_CONST::n> &Kz,
    const Eigen::Vector<double,CTR_CONST::n> &Ux,
    const Eigen::Vector<double,CTR_CONST::n> &l,
    const Eigen::Vector<double,CTR_CONST::n> &l_k,
    const Eigen::Vector3d &f,

    Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> &yTot_out,
    Eigen::Vector<double,CTR_CONST::NB_BC> &b_out,
    Eigen::Matrix<double,6,CTR_CONST::NB_Q> &Eq_out,
    Eigen::Matrix<double,6,CTR_CONST::NB_YU0> &Eu_out,
    Eigen::Matrix<double,CTR_CONST::NB_BC,CTR_CONST::NB_Q> &Bq_out,
    Eigen::Matrix<double,CTR_CONST::NB_BC,CTR_CONST::NB_YU0> &Bu_out,
    uint nThread);
#endif //IVP_FD_H
