#ifndef LOAD_PARAMETERS_H
#define LOAD_PARAMETERS_H

#include <string>
#include <vector>
#include "Eigen/Dense"

struct parameters{
    // Meta-data
    int id;                               // identifier (useful when running multiples tests with different parameters)
    std::string name;                     // name of the set of parameters
    std::string comment;                  // description of the set of parameters

    // Material and geometric parameters
    Eigen::Vector3d E;                    // Young's modulus of each tube [Pa]
    Eigen::Vector3d mu;                   // Poisson ratio of each tube [/]
    Eigen::Vector3d Ux;                   // Precurvature along x-axis  for each tube [m^-1]
    Eigen::Vector3d l;                    // Effective length for each tube [m]
    Eigen::Vector3d l_k;                  // Precurved length for each tube [m]
    Eigen::Vector3d rIn;                  // Inner radius of each tube [m]
    Eigen::Vector3d rOut;                 // Outer radius of each tube [m]

    // <UNUSED IN THIS CURRENT VERSION> Additional parameters when using an experimental setup 
    Eigen::Vector<double,6> Tref;         // Homogeneous transformation matrix between the reference sensor frame and the robot bas frame [m for the translation part]
    Eigen::Vector3d Ttip;                 // Translation between the origin of the end-effector sensor and the origin of the end effector frame (expressed in the end-effector sensor frame) [m]
    Eigen::Vector<double,6> offset;       // Translational and rotational offsets (between the "actuator zero position" and the "zero position" defined in the model) [m ; m; m; deg; deg; deg]
    
    // Control parameters
    Eigen::Vector<double,8> w;            // Weights for the controller (w0t, w0r, w1t, w1r, w2t, w2r) []
    double lambda;                        // Gain for the controller []
    double trajSpeedMax;                  // Maximal velocity for trapezoidal velocity generation [m/s]
    double trajAccMax;                    // Maximal acceleration for trapezoidal velocity generation [m/s^2]

    Eigen::Vector3d force;                // External point force applied at the end-effector [N]
    double noiseStd;                      // Standard deviation of the measurement noise [m]
  };
void printParameter(parameters &p);
int loadParameters (const std::string &path, std::vector<parameters> &vParameters);

#endif //LOAD_PARAMETERS_H