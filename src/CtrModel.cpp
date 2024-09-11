#include <iostream>
#include "Eigen/Dense"
#include "CtrModel.h"
#include "ComputeIvpJacobianMatrices.h"
#include "pinv.h"

#include "segmenting.h"
#include "solveIVP.h"
#include "bcError.h"

#include <chrono>
#include <omp.h>

#define rt_printf(x) printf(x)

using namespace CTR_CONST;
using namespace Eigen;

constexpr int nIterationsMax = 1000;
constexpr double maxNormB = 1e-10;

int nThread = omp_get_max_threads();

CtrModel::CtrModel(parameters p){

    const Vector<double,3> G = p.E.cwiseQuotient(2 * (p.mu + Vector<double,3>::Ones())); // shear modulus
    const Vector<double,n> I = (pi/4) * (p.rOut.array().pow(4) - p.rIn.array().pow(4)); // 2nd moment of inertia
    const Vector<double,n> J = 2.0 * I;   // polar moment of inertia

    Kxy = p.E.cwiseProduct(I); // bending stiffness
    Kz  =   G.cwiseProduct(J); // torsionnal stiffness

    Ux  = p.Ux; // precurvature along x-axis
    l   = p.l; // tube length
    l_k = p.l_k; // precurved length

    offset = p.offset;
    offset(seqN(n,n)) *= 180.0 * pi; // convert rotational offset from degree to radian

    // Set arbitrary initial configuration close to the "fully deployed" configuration with a small margin between each tube base
    constexpr double margin = 5e-3; 
    q = offset + Vector<double, 2 * n>(-3.0 * margin, -2.0 * margin, -margin, 0.0, 0.0, 0.0);

    f = p.force; // external point force applied at the end-effector
    yu0 = Vector<double, NB_YU0>::Zero();

    Compute<LOAD_J>(q); // Compute the model for the initial configuration to initialize yu0, X, and J.
}

CtrModel::~CtrModel(){
  //dtor
}
template <COMPUTATION_OPTION opt> 
int CtrModel::Compute(Vector<double,NB_Q> argQ){
    int nbIteration = 0;
    q = argQ;
    Vector<double,NB_YU0> yu0_tilde = yu0;

    Matrix<double,nStateVar,nSegMax*nIntPoints> yTot_out;
    Vector<double,NB_BC> b_out;
    Matrix<double,6,CTR_CONST::NB_Q> Eq_out;
    Matrix<double,6,CTR_CONST::NB_YU0> Eu_out;
    Matrix<double,CTR_CONST::NB_BC,CTR_CONST::NB_Q> Bq_out;
    Matrix<double,CTR_CONST::NB_BC,CTR_CONST::NB_YU0> Bu_out;

    int nIter = 0;
    
    // first iteration without J and/or C, just to compute the BC residuals and the Bu matrix
    if(ComputeIvpJacobianMatrices<LOAD>(q,yu0_tilde,Kxy,Kz,Ux,l,l_k,f,yTot_out,b_out,Eq_out,Eu_out,Bq_out,Bu_out, nThread) != 0){
      std::cout << "CtrModel::Compute()>> ComputeIvpJacobianMatrices() returned non-zero !" << std::endl;
      return -1;
    }
    yu0_tilde -= pinv(Bu_out) * b_out; // update guess using Gauss-Newton
    nIter++;
    nbIteration++;
    do{ // first pass with yu0 from previous computation
      if(ComputeIvpJacobianMatrices<opt>(q,yu0_tilde,Kxy,Kz,Ux,l,l_k,f,yTot_out,b_out,Eq_out,Eu_out,Bq_out,Bu_out, nThread) != 0){
        std::cout << "CtrModel::Compute()>> ComputeIvpJacobianMatrices() returned non-zero !" << std::endl;
        return -1;
      }
      yu0_tilde -= pinv(Bu_out) * b_out; // update guess using Gauss-Newton
      nIter++;
      nbIteration++;
    }while(b_out.norm() > maxNormB && nIter < nIterationsMax);
    if(b_out.hasNaN() || b_out.norm() > maxNormB){ // if the first pass fails
      std::cout << "CtrModel::Compute()>> Model failed to converge after " << nbIteration << " iterations."
        "Potential bifurcation detected, switching sign of initial guess for torsions and trying again." << std::endl;
      //switch sign of torsion but keep x,y curvature
      yu0_tilde(0) = yu0(0);
      yu0_tilde(1) = yu0(1);
      yu0_tilde(2) = -yu0(2);
      yu0_tilde(3) = -yu0(3);
      yu0_tilde(4) = -yu0(4);
      nIter = 0;
      do{ //second pass, with opposite torsion
        if(ComputeIvpJacobianMatrices<opt>(q,yu0_tilde,Kxy,Kz,Ux,l,l_k,f,yTot_out,b_out,Eq_out,Eu_out,Bq_out,Bu_out, nThread) != 0){
          std::cout << "CtrModel::Compute()>> ComputeIvpJacobianMatrices() returned non-zero !" << std::endl;
          return -1;
        }
        yu0_tilde -= pinv(Bu_out) * b_out; // update guess using Gauss-Newton
        nIter++;
        nbIteration++;
      }while(b_out.norm() > maxNormB && nIter < nIterationsMax);
    }
    if(b_out.hasNaN() || b_out.norm() > maxNormB){ // if the second pass fail
      std::cout << "CtrModel::Compute()>> Model failed to converge even when switching sign of torsions !" << std::endl;
      return -1;
    }
    
    // The model converged sucessfully : Update member variables 

    yu0 = yu0_tilde;

    segmentedData segmented_out;
    if(segmenting(q,Kxy,Kz,Ux,l,l_k,segmented_out)!=0){
      std::cout << "CtrModel::Compute()>> segmenting returned non-zero !" << std::endl;
      return -1;
    }
    segmented = segmented_out;

    yTot = yTot_out;
    X = getXFromYtot(yTot_out,segmented_out);
    R = getRvFromYtot(yTot_out,segmented_out).reshaped(3,3);

    if(opt == LOAD_J){
      Matrix<double,6,6> RR = Matrix<double,6,6>::Zero();
      RR(seq(0,2),seq(0,2)) = R;
      RR(seq(3,5),seq(3,5)) = R;

      J = RR * (Eq_out - Eu_out * pinv(Bu_out) * Bq_out); 
    }

  return nbIteration;
}

template int CtrModel::Compute<LOAD>(Vector<double,NB_Q> argQ);
template int CtrModel::Compute<LOAD_J>(Vector<double,NB_Q> argQ);