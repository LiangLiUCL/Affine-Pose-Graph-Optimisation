#ifndef AFFINE_VERTEX_H_
#define AFFINE_VERTEX_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h> 
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "AFFINE_LIE.h"

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/eigen_types.h>
#include "g2o/solvers/csparse/linear_solver_csparse.h"

using namespace std;

typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> Vector6d;
typedef Eigen::Matrix<double,6,6> Matrix6d;
typedef Eigen::Transform<double, 2, Eigen::Affine, Eigen::ColMajor> AffineE;



class VertexAffine2d: public g2o::BaseVertex<6, AFFINE_LIE>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexAffine2d(){}
    bool read (istream& is)
    {
        double affMat[6];
        for ( int i=0; i<6; i++ )
            is>>affMat[i];
        setEstimate(AFFINE_LIE(affMat[0], affMat[1], affMat[2], affMat[3], affMat[4], affMat[5]));
        return true;
    }

    bool write (ostream& os) const
    {
        //os<<id()<< " ";
        Eigen::Vector2d trans = _estimate.translation();
        Eigen::Matrix2d genLinear = _estimate.generalLinear();

        os<<genLinear(0,0)<<" "<<genLinear(0,1)<<" "<<genLinear(1, 0)<<" "<<genLinear(1,1)<<" "<<trans[0]<<" "<<trans[1]<<" "<<endl;
        return true;
    }

    virtual void setToOriginImpl()
    {
        _estimate = AFFINE_LIE();
    }

    virtual void oplusImpl ( const double* update_ )
    {
        Eigen::Map<Vector6d> update(const_cast<double*>(update_));
        AFFINE_LIE s(update);

        _estimate = s*_estimate;
    }

    /*virtual bool setEstimateDataImpl(const double* est)
    {
        _estimate = AFFINE_LIE(est[0], est[1], est[2], est[3], est[4], est[5]);
        return true;
    }

    virtual bool getEstimateData(double* est) const {
        Vector6d v(est);
        v = _estimate.toVector();
        return true;
      }
      
      virtual int estimateDimension() const { return 6; }

      virtual bool setMinimalEstimateDataImpl(const double* est){
        return setEstimateData(est);
      }

      virtual bool getMinimalEstimateData(double* est) const {
        return getEstimateData(est);
      }*/
};

#endif