#ifndef CTR_CONSTANTS_H
#define CTR_CONSTANTS_H

namespace CTR_CONST{
    constexpr int n = 3;                        // number of tubes
    constexpr int nIntPoints = 2;               // number of integration nodes per segment
    constexpr int nSegMax = 2 * n;              // max number of segments
    constexpr int nStateVar = 15 + 2 * (n - 1); // number of state variables

    constexpr int NB_Q = 2 * n;                 // number of actuation variables
    constexpr int NB_YU0 = n + 2;               // number of unknown initial state variables
    constexpr int NB_BC = n + 2;                // number of variables in boundary condition residuals 

    constexpr double pi = 3.1415926535897932384626433832795028841971L;
}

struct segmentedData{
    Eigen::Vector<double, CTR_CONST::nSegMax>              S;    // vector containing the arc length at the end of each segment
    Eigen::Vector<int,    CTR_CONST::n>                    iEnd; // vector of the segment index where each tube ends
    Eigen::Matrix<double, CTR_CONST::n,CTR_CONST::nSegMax> KKxy; // matrix of bending stiffness for each tube (row) at each segment (column)
    Eigen::Matrix<double, CTR_CONST::n,CTR_CONST::nSegMax> KKz;  // matrix of torsional stiffness for each tube (row) at each segment (column)
    Eigen::Matrix<double, CTR_CONST::n,CTR_CONST::nSegMax> UUx;  // matrix of precurvature for each tube (row) at each segment (column)
    };

// templated computation options
// NO_LOAD : Compute only the model, unloaded
// NO_LOAD_J : Compute the model along with the robot Jacobian matrix J, unloaded
// LOAD : Compute only the model, with external loads
// LOAD_J : Compute the model along with the robot Jacobian matrix J, with external loads
// LOAD_J_C : Compute the model along with the robot Jacobian matrix J and compliance matrix C, with external loads
//enum COMPUTATION_OPTION {NO_LOAD, NO_LOAD_J, LOAD, LOAD_J, LOAD_J_C};
enum COMPUTATION_OPTION {LOAD, LOAD_J};
// Warning ! This is a temporary version of the code. Computation options NO_LOAD, NO_LOAD_J, and LOAD_J_C will be implemented shortly.

// Convert index of segment (in SegmentedData) to column index in yTot matrix
inline int getYtotIndexFromIend(int iEnd){return iEnd * CTR_CONST::nIntPoints + CTR_CONST::nIntPoints - 1;}
// get position X from yTot
inline Eigen::Vector3d getXFromYtot(
    const Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> &yTot,
    const segmentedData &segmented){
        
    return yTot(Eigen::seqN(0,3), getYtotIndexFromIend(segmented.iEnd(0)));
};
// get "flattened" 9-element vector representing R from yTot
inline Eigen::Vector<double,9> getRvFromYtot(
    const Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> &yTot,
    const segmentedData &segmented){
        
    return yTot(Eigen::seqN(3,9), getYtotIndexFromIend(segmented.iEnd(0)));
};

#endif //CTR_CONSTANTS_H