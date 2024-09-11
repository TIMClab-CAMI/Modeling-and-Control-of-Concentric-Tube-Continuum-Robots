#include <iostream>
#include "CtrModel.h"

using namespace Eigen;
using namespace CTR_CONST;

int main(int, char**){
  // Example 1 : simply compute the model

  // Load parameters corresponding to CTR
  std::vector<parameters> vParameters;
	if(loadParameters("../parameters/parameters.csv",vParameters) != 0){
    return -1;
  }
  parameters &pNominal =  vParameters[0];

  CtrModel ctr(pNominal);
  // Declare actuation variables q
  Vector<double,NB_Q> q;
  q << -0.3, -0.2, -0.1, 0, 0, 0; // arbitrary initial configuration
  // Compute model
  ctr.Compute<LOAD_J>(q);

  // Get end-effector position
  Vector<double,3> X = ctr.GetX();
  std::cout << "X = [" << X.transpose() << "]" << std::endl;
  return 0;
}
