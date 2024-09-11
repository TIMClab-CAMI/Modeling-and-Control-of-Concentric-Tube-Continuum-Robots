#ifndef CTR_MODEL_H
#define CTR_MODEL_H
#include "Eigen/Dense"
#include "ctrConstants.h"
#include "loadParameters.h"

class CtrModel
{
public:
    CtrModel(parameters p);
    virtual ~CtrModel();

    template <COMPUTATION_OPTION opt> 
    int Compute(Eigen::Vector<double,CTR_CONST::NB_Q> q); // Compute end-effector position

    const Eigen::Vector<double, 2 * CTR_CONST::n> GetOffset(){return offset;};
    const Eigen::Vector<double,CTR_CONST::n> GetUx(){return Ux;};
    const Eigen::Vector<double,CTR_CONST::n> GetL(){return l;};
    const Eigen::Vector<double,CTR_CONST::n> GetL_k(){return l_k;};
    const Eigen::Vector<double,CTR_CONST::n> GetKxy(){return Kxy;};
    const Eigen::Vector<double,CTR_CONST::n> GetKz(){return Kz;};
    const Eigen::Vector3d GetF(){return f;};
    const Eigen::Matrix<double,CTR_CONST::nStateVar,CTR_CONST::nSegMax*CTR_CONST::nIntPoints> GetYTot(){return yTot;}
    const Eigen::Vector<double,CTR_CONST::NB_YU0> GetYu0(){return yu0;};
    const Eigen::Vector3d GetX(){return X;};
    const Eigen::Matrix<double,6,CTR_CONST::NB_Q> GetJ(){return J;};

    segmentedData segmented;

protected:
    Eigen::Vector<double,CTR_CONST::n> Ux;         // Precurvature along x-axis (m^-1) for each tube
    Eigen::Vector<double,CTR_CONST::n> l;          // Effective length for each tube [m]
    Eigen::Vector<double,CTR_CONST::n> l_k;        // Precurved length for each tube [m]
    Eigen::Vector<double,CTR_CONST::n> Kxy;        // Bending stiffness E*I (Pa.m^4) for each tube 
    Eigen::Vector<double,CTR_CONST::n> Kz;         // Torsional stiffnes G*J (Pa.m^4) for each tube
    Eigen::Vector<double,CTR_CONST::NB_Q> q;       // Actuation variables [beta_i, alpha_i] with beta the translations (m) and alpha the rotations (rad) actuation
    Eigen::Vector<double,CTR_CONST::NB_BC> b;      // Boundary condition errors
    Eigen::Vector3d f;                             // External, punctual force applied on the end-effector
    Eigen::Matrix<double,CTR_CONST::nStateVar,
        CTR_CONST::nSegMax*CTR_CONST::nIntPoints> yTot; // Matrix containing state variables along arc-length (at each integration node, as defined in "ctrConstants.h")
    Eigen::Vector<double,CTR_CONST::NB_YU0> yu0;   // Guess for unknown initial state variables
    Eigen::Vector3d X;                             // Position of the end-effector
    Eigen::Matrix3d R;                             // Orientation of the end-effector as a rotation matrix

    Eigen::Matrix<double,6,CTR_CONST::NB_Q> J;    // Robot Jacobian

    Eigen::Vector<double, 2 * CTR_CONST::n> offset; // Translational and rotational offsets (between the "actuator zero position" and the "zero position" defined in the model) [m ; m; m; rad; rad; rad]
private:
    

};

#endif // CTR_MODEL_H
