#ifndef SEGMENTING_H
#define SEGMENTING_H

#include "Eigen/Dense"
#include "ctrConstants.h"

/**
 * @brief Split the CTR in several segments. A new segment start each time a tube starts, end or at each step change or precurvature.
 * 
 * @param[in] q Vector of actuation variables [beta_i, alpha_i] with beta the translations [m] and alpha the rotations [rad] 
 * @param[in] Kxy Vector of bending stiffness for each tube
 * @param[in] Kz Vector of torsional stiffness for each tube
 * @param[in] Ux Vector of precurvature for each tube
 * @param[in] l Vector of effective length for each tube
 * @param[in] l_k Vector of precurved length for each tube (we assume a tube is composed of a straight part and a precurved part with constant precurvature)
 
 * @param[out] segmented_out    Output structure containing informations about the segmented CTCR
 * @return int (0 indicates success and < 0 indicates failure )
 */
int segmenting(
  const Eigen::Vector<double,CTR_CONST::NB_Q>            &q,
  const Eigen::Vector<double,CTR_CONST::n>               &Kxy,
  const Eigen::Vector<double,CTR_CONST::n>               &Kz,
  const Eigen::Vector<double,CTR_CONST::n>               &Ux,
  const Eigen::Vector<double,CTR_CONST::n>               &l,
  const Eigen::Vector<double,CTR_CONST::n>               &l_k,

  segmentedData                                          &segmented_out);
#endif //SEGMENTING_H
