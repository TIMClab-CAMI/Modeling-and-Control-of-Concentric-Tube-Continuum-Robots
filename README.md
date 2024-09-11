## Warning !
This is a temporary version of the code. It is functional yet, the documentation might be incomplete and only simple demonstration codes are available. A code revision will come shortly, providing demonstration codes to reproduce the results of our paper.

# Modeling and Control of Concentric Tube Continuum Robots 
## This code is associated to the following paper, which we request to be explicitly cited in all forms of communication of your work:

> Quentin Boyer, Sandrine Voros, Pierre Roux, François Marionnet, Kanty Rabenorosoa and M. Taha Chikhaoui, "On High Performance Control of Concentric Tube Continuum Robots Through Parsimonious Calibration," in IEEE Robotics and Automation Letters, doi: 10.1109/LRA.2024.3455906

It provides an efficient numerical implementation of the kineto-static torsionally compliant model for Concentric Tube Continuum Robots (CTCRs), including external loads, as described in the paper Rucker et al., “A Geometrically Exact Model for Externally Loaded Concentric-Tube Continuum Robots”, in IEEE Transactions on Robotics, along with the computation of the robot Jacobian matrix and compliance matrix. This model is implemented in C++ and uses Eigen library for Vector and Matrix operations, and OpenMP for parallel computation. 
A closed-loop control scheme based on our model implementation is then implemented, using the Generalized Damped Least Squares (GDLS) method to account for several concurrent tasks.
For more details, please refer to the paper available at: https://hal.science/hal-04685717v1.

![Github_figTrajObstacle_simu](https://github.com/user-attachments/assets/8174e7dd-2640-455a-84cc-a94f5881f410)

## Modeling

The CTCR model is formulated as a multipoint Boundary Value Problem (BVP). The resolution is implemented using a shooting method.
The robot is first divided into segments delimited by the ends of the tubes and by tube precurvature discontinuities. The Initial Value Problem (IVP) is solved using a single step of 4-th order Runge-Kutta integration scheme for each segment. The values of the state variables after the discontinuities are algebraically computed using the values before the discontinuities and the continuity equations. The BVP is solved using the Gauss-Newton algorithm, and the robot Jacobian matrix and compliance matrix can optionally be computed along the model. The use of templated functions allows heavy optimization by determining which parts of the model are required by the user at compile-time. Parallel computing is used to compute the jacobian matrices of the IVP based on the finite differences method. The model implementation does not use dynamic allocation and can easily be adopted in a real-time context. We validated it on Linux Xenomai 3.1 for our experimental setup.

## Control

A closed-loop control scheme based on our model implementation is then implemented. Using the (GDLS) method, several concurrent tasks are included in the control scheme, such as trajectory tracking, damping, actuation constraints, and obstacle avoidance. The control performances are demonstrated in simulation for a 3D trajectory representing a significant portion of the robot workspace, with an obstacle along the trajectory, and in presence of external forces.

## Structure of the code
* The model implementation is provided as a library, which sources are in the "src" folder.
    * The "segmenting" function divides the CTCR into segments delimited by the ends of the tubes and by tube precurvature discontinuities.
    * The "odeCtr" function implements the Ordinary Differential Equations.
    * The "odeIntCtrRK4" function uses the classical Runge-Kutta 4 integration scheme to integrate the ODEs along one segment.
    * The "solveIVP" function solves the IVP by computing the forward integration for each segment and algebraically calculating the state variable values for the next segment based on the values of the previous segment and the continuity equations.
    * The "ComputeIvpJacobianMatrices" computes the IVP jacobian matrices using the finite differences method and parallel computation.
    * The "bcError" function computes the boundary conditions residuals.
    * The "Compute" function in the class "CtrModel" solves the BVP using the Gauss-Newton algorithm to compute the whole shape of the robot, along with (optionally) the robot Jacobian matrix and the compliance matrix.
* Demonstation codes are located in the "demo" folder and show how to use the library in your own application.
    * The demo "001_model_computation" is a simple example of how to compute the kineto-static model.
    * The demo "002_control_line" is a simple control example, performing trajectory tracking of the end-effector in position along a straight line.
    * Additional demonstration codes will be provided soon to reproduce the results shown in our paper.
 
## Prerequisites
This code has been tested on Ubuntu 20.04, Ubuntu 22.04, and Windows 10 using WSL2 with the Ubuntu 22.04 distribution. It will be tested for Windows with Visual Studio in the incoming code revision.
* CMake
* Eigen (version 3.4.0 is included in the repository)
* Gnuplot and Boost are only used to plot the CTCR in the demonstrations. They can be removed if no plotting is required.

(Linux) Install the following packages :  cmake build-essential libboost-all-dev gnuplot
```sh
sudo apt install cmake build-essential libboost-all-dev gnuplot
```

## Configure and build
Create a build folder in the repository, and run CMake in it.
```sh
mkdir build
cd build
cmake ..
```
Compile the library and the demonstation codes
```sh
make -j
```
Execute the demonstation code
```sh
demo/001_model_computation
```

## Licence
This project is licensed under the GPL v3.0 License - see the [LICENSE](https://github.com/TIMClab-CAMI/Modeling-and-Control-of-Concentric-Tube-Continuum-Robots/blob/main/LICENSE) file for details

## Contributing
Feel free to submit pull requests and use the issue tracker to start a discussion about any bugs you encounter. Please provide detailed description of the versions of your operating system, tools, and libraries for any software related bugs.

## Acknowledgements
This work was supported by grants ANR-11-LABX-0004-01, ANR-21-ESRE-0015, ANR-17-CE19-0005, ANR-17-EURE-0002, and ANR-20-CE33-0001.
