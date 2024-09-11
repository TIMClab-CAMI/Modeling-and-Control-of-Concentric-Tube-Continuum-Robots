#include <numeric>      // std::iota
#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "ctrConstants.h"
#include "segmenting.h"

using namespace Eigen;
using namespace CTR_CONST;

void sort_indexes(
  const Vector<double,2*n+1>  &v,
  Vector<double,2*n+1>  &v_sorted_out,
  Vector<int,2*n+1>  &index_out){

  int n = v.size();
  std::vector<int> idx(n);
  iota(idx.begin(), idx.end(), 0);

  std::stable_sort(idx.begin(), idx.end(),
    [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  for(int i = 0; i < n; i++){
    index_out(i) = idx[i];
    v_sorted_out(i) = v(idx[i]);
  }
}
void diff(
  const Vector<double,2*n+1>  &v,
  Vector<double,2*n>  &diff_v_out){

  int n = v.size();
  for(int i = 0; i < n-1; i++){
    diff_v_out(i) = v(i+1)-v(i);
  }
}
int segmenting(
  const Vector<double,NB_Q> &q,
  const Vector<double,n> &Kxy,
  const Vector<double,n> &Kz,
  const Vector<double,n> &Ux,
  const Vector<double,n> &l,
  const Vector<double,n> &l_k,

  segmentedData &segmented_out){

  // out : [S,iEnd,EE,II,GG,JJ,UUx]
  Vector<double,n> B;
  B << q(Eigen::seqN(0,n));
  Vector<double,n> d1;
  d1 = l+B;      //% position of tip of the tubes
  Vector<double,n> d2 = d1-l_k;   //% position of the point where tube bending starts
  Vector<double,2*n+1> points; // number of segment delimitations = 1 + d2.size() + d1.size()
  points << 0, d2, d1;

  Vector<double,2*n+1> points_srt;
  Vector<int,2*n+1> index;

  sort_indexes(points,points_srt,index);

  Vector<double,2*n> L;
  diff(points_srt,L);
  for(int i = 0; i < 2*n; i++){
    L(i) = 1e-5 * floor(L(i)*1e5);
  }
  for (int i = 0; i < n-1; i++){
    if (B(i)>B(i+1)){
      std::cout << "proximal end of tube " << i+1 << " is clashing into tube " << i+2 << std::endl;
      std::cout << "B = " << std::endl << B << std::endl;
      std::cout << "q = " << std::endl << q << std::endl;
      return -1;
    }
    if(d1(i)<d1(i+1)){
      std::cout << "distal end of tube " << i+1 << " is clashing into tube " << i+2 << std::endl;
      std::cout << "d1 = " << std::endl << d1 << std::endl;
      std::cout << "q = " << std::endl << q << std::endl;
      return -1;
    }
  }
  int i0;
  Vector<int,n> iCurved;
  segmented_out.iEnd = Vector<int,n>::Zero();
  for (int i = 0; i < n; i++){
    for (int j = 0; j < 2*n+1; j++){
      if(index(j)==0){
        i0 = j;
      }
      else if(index(j) == i+1){
        iCurved(i) = j;
      }
      else if(index(j) == n+i+1){
        segmented_out.iEnd(i) = j;
      }
    }

  }
  /*Optionnal zero-padding : just for easier debugging*/
  segmented_out.S = Vector<double,nSegMax>::Zero();
  segmented_out.KKxy = Matrix<double,n,nSegMax>::Zero();
  segmented_out.KKz = Matrix<double,n,nSegMax>::Zero();
  segmented_out.UUx = Matrix<double,n,nSegMax>::Zero();
  double sumL = 0;
  int iSegment = 0;
  for (int j = i0; j < 2*n; j++){
      sumL += L(j);
      segmented_out.S(iSegment) = sumL;
      for (int i = 0; i < n; i++){
        if(j < segmented_out.iEnd(i)){ // if tube exists
          segmented_out.KKxy(i,iSegment) = Kxy(i);
          segmented_out.KKz(i,iSegment) = Kz(i);
        }
        if(j >= iCurved(i) && j < segmented_out.iEnd(i)){ // if tube is curved
          segmented_out.UUx(i,iSegment) = Ux(i);
        }
      }
      iSegment++;
  }
  segmented_out.iEnd.array() -= (i0+1); // remove segments before 0 so iEnd correspond to the index of the last segment where ith tube exists
  return 0;
}
