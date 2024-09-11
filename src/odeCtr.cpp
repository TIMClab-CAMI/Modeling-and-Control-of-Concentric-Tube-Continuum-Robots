#include <iostream>
#include "Eigen/Dense"
#include "ctrConstants.h"
#include "odeCtr.h"

using namespace Eigen;
using namespace CTR_CONST;

#define UNUSED(expr) (void)(expr)

void odeCtr(
    double s,
    Vector<double,nStateVar>& y,
    const Vector<double,n>& Kxy,
    const Vector<double,n>& Kz,
    const Vector<double,n>& Ux,
    const Vector3d& f,

    Vector<double,nStateVar>& y_s_out){


    UNUSED(s); // to remove the "unused parameter" compiler warning
    // s is not used here with a point force at the end-effector but it is needed in the more general case

    //Unpack state vector
    // 0:2  : r
    // 3:11 : R
    // 12:14 : u1
    // 14:13+n (14:16): uz (14:u1z both in "u1" and in "uz")
    // 14+n (17:18): 12+2*n : theta(i+1) (theta1 = 0)

    Matrix3d R = Map<Matrix<double,3,3,RowMajor> >(&y[3]);
    Vector3d u1 = Map<Vector3d>(&y[12]);
    Vector3d uz = Map<Vector3d>(&y[14]);
    Vector3d theta(0, y[17], y[18]);

    //Pack state vector derivative
    Map<Vector3d> r_s(&y_s_out[0]);
    Map<Matrix<double,3,3,RowMajor>> R_s(&y_s_out[3]);
    Map<Vector2d> u1xy_s(&y_s_out[12]);
    Map<Vector3d> uz_s(&y_s_out[14]);

    r_s = R.col(2); // = R * e3
    R_s = R * hat(u1); 

    Vector3d thetha_s;
    thetha_s << 0, uz(1) - uz(0), uz(2) - uz(0);
    Vector3d e3;
    e3 << 0,0,1;
    
    Matrix3d K_i[n];
    Vector3d Sum_i[n];
    for(int i = 0; i < n; i++){

        K_i[i] = Matrix3d(Vector3d(Kxy(i), Kxy(i),Kz(i)).asDiagonal());

        double ci = cos(theta(i));
        double si = sin(theta(i));

        Matrix3d R_theta_i;
        R_theta_i <<  ci, -si, 0,
                      si, ci,  0,
                      0,  0,   1;

        Matrix3d d_R_thetai_t___dTheta_i;
        d_R_thetai_t___dTheta_i <<  -si,  ci,  0,
                                    -ci,  -si, 0,
                                    0,    0,   0;

        Vector3d u_i = R_theta_i.transpose() * u1 + thetha_s(i) * e3;
        Vector3d U_i;
        U_i << Ux(i), 0, 0;
        Sum_i[i] = R_theta_i * ( K_i[i] * ( thetha_s(i) * d_R_thetai_t___dTheta_i * u1) + hat(u_i)*K_i[i] * (u_i - U_i ));

        Vector<double,n> Kz2 = Kz;
        if (Kz2(i)==0){
            Kz2(i)=1;//  % to avoid dividing by zero when tube doesn't exist
        }
        uz_s(i) = Kxy(i) / Kz2(i) * (-u_i(1)) * Ux(i);
    }
    Matrix3d K = K_i[0] + K_i[1] + K_i[2];
    Vector3d SUM = Sum_i[0] + Sum_i[1] + Sum_i[2];
    Matrix3d i_K = Matrix3d(Vector3d(1/K(0,0),1/K(1,1),1/K(2,2)).asDiagonal());
    
    // For a point force (modeled as a Dirac distribution at the tip), the integral of the force distribution is simply the magnitude of the point force
    Vector3d forceIntegral = f; 
    Vector3d u1xy_s_3d = -i_K * SUM -i_K * hat(e3) * R.transpose() * forceIntegral;
    
    u1xy_s = u1xy_s_3d(seqN(0,2));
    y_s_out[17] = thetha_s(1);
    y_s_out[18] = thetha_s(2);

}