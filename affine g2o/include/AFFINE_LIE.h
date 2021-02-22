#ifndef AFFINE_LIE_H_
#define AFFINE_LIE_H_



#include <iostream>
#include <math.h>
#include "g2o/stuff/misc.h"
#include <Eigen/Core>
#include <Eigen/Geometry>


using namespace std;

typedef Eigen::Matrix<double, 4, 1, Eigen::ColMajor> Vector4d;
typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> Vector6d;
//typedef Eigen::Matrix<double, 2, 2, Eigen::ColMajor> Matrix2d;


double computeFactorial(int input)
{
  double res;
  res = 1.;
  if (input == 0)
  {
    return res;
  }
  for (int i=1; i<input+1; i++)
  {
    res *= i;
  }
  
  return res;
}

Eigen::Matrix3d expo(Eigen::Matrix3d matIn, int ord)
{
  Eigen::Matrix3d res;
  res.setIdentity();
  Eigen::Matrix3d input;
  input.setIdentity();
  for (int i=1; i<ord+1; i++)
  {
    //cout << i << endl;
    res += (1./computeFactorial(i))*input*matIn;
    input = input*matIn;
  }

  return res;
}

Eigen::Matrix3d logarit(Eigen::Matrix3d matIn, int ord)
{
  Eigen::Matrix3d res;
  res << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  Eigen::Matrix3d input;
  input.setIdentity();
  Eigen::Matrix3d ind;
  ind.setIdentity();
  for (int i=1; i<ord+1; i++)
  {
    res = res + (pow(-1, i+1)/i)*input*(matIn - ind);
    input = input*(matIn -ind);
  }

  return res;
}


class AFFINE_LIE
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    AFFINE_LIE()
    {
      _sl.setIdentity();
      _trans.fill(0.);
    }
    //constructor from general linear part and translation
    AFFINE_LIE(const Eigen::Matrix2d& sl_, const Eigen::Vector2d& t_): _sl(sl_), _trans(t_){}
    //constructor from elements
    AFFINE_LIE(const double& x11, const double& x12, const double& x21, const double& x22, const double& t_x, const double& t_y)
    {
      _sl << x11, x12,
             x21, x22;
      _trans << t_x, t_y;
    }


    AFFINE_LIE(const Eigen::Affine2d& affIni)
    {
      _sl = affIni.linear();
      _trans = affIni.translation();
    }

    //constructor from Lie algebra exponential map
    AFFINE_LIE(const Vector6d& update)
    {
      Vector4d sl_alge;
      for (int i=0; i<4; i++)
        sl_alge[i] = update[i];
      Eigen::Vector2d trans_alge;
      for (int i=0; i<2; i++)
        trans_alge[i] = update[i+4];
      
      Eigen::Matrix3d mat_alge;
      mat_alge << sl_alge[0]+sl_alge[1], sl_alge[2]+sl_alge[3], trans_alge[0],
                  sl_alge[2]-sl_alge[3], sl_alge[0]-sl_alge[1], trans_alge[1],
                  0.,0.,0.;
      
      Eigen::Matrix3d expoMat;
      expoMat = expo(mat_alge, 10);
      _sl = expoMat.block<2, 2>(0, 0);
      _trans = expoMat.block<2, 1>(0, 2);
    }

    //logarithm map
    Vector6d log() const
    {
      Vector6d res;
      Eigen::Matrix3d affState;
      affState.setIdentity();
      affState.block<2, 2>(0, 0) = _sl;
      affState.block<2, 1>(0, 2) = _trans;
      Eigen::Matrix3d affRes;
      affRes = logarit(affState, 10000);
      
      res[0] = (affRes(0, 0) + affRes(1, 1))/2.;
      res[1] = (affRes(0, 0) - affRes(1, 1))/2.;
      res[2] = (affRes(0, 1) + affRes(1, 0))/2.;
      res[3] = (affRes(0, 1) - affRes(1, 0))/2;
      res[4] = affRes(0, 2);
      res[5] = affRes(1, 2);

      return res;
    }

    Eigen::Vector2d map(const Eigen::Vector2d& xy ) const
    {
      return _sl*xy+_trans;
    }

    AFFINE_LIE inverse() const
    {
      return AFFINE_LIE(_sl.inverse(), -_sl.inverse()*_trans);
    }

    double operator[](int i) const
    {
      assert(i<6);
      if (i<4)
      {
        return _sl.data()[i];
      }
      if (i<6)
      {
        return _trans[i-4];
      }
    }

    double& operator[](int i)
    {
      assert(i<6);
      if (i<4)
      {
        return _sl.data()[i];
      }
      if (i<6)
      {
        return _trans[i-4];
      }
    }

    AFFINE_LIE operator *(const AFFINE_LIE& other) const
    {
      AFFINE_LIE ret;
      ret._sl = _sl*other._sl;
      ret._trans = _sl*other._trans + _trans;
      return ret;
    }

    AFFINE_LIE& operator *=(const AFFINE_LIE& other)
    {
      AFFINE_LIE ret = (*this)*other;
      *this = ret;
      return *this;
    }

    inline const Eigen::Vector2d& translation() const {return _trans;}

    inline Eigen::Vector2d& translation() {return _trans;}

    inline const Eigen::Matrix2d& generalLinear() const {return _sl;}

    inline Eigen::Matrix2d& generalLinear() {return _sl;}

protected:
  Eigen::Matrix2d _sl;
  Eigen::Vector2d _trans;

};


#endif