#ifndef AFFINE_EDGE_H_
#define AFFINE_EDGE_H_

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


class EdgeAffine2d: public g2o::BaseBinaryEdge<6, AFFINE_LIE, VertexAffine2d, VertexAffine2d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //EdgeAffine2d(){}
    bool read (istream& is)
    {
        double affMat[6];
        for ( int i=0; i<6; i++ )
            is>>affMat[i];
        setMeasurement(AFFINE_LIE(affMat[0], affMat[1], affMat[2], affMat[3], affMat[4], affMat[5]));
        for ( int i=0; i<information().rows() && is.good(); i++ )
            for ( int j=i; j<information().cols() && is.good(); j++ )
            {
                is >> information() ( i,j );
                if ( i!=j )
                    information() ( j,i ) = information() ( i,j );
            }
        return true;
    }

    bool write ( ostream& os ) const
    {
        VertexAffine2d* v1 = static_cast<VertexAffine2d*> (_vertices[0]);
        VertexAffine2d* v2 = static_cast<VertexAffine2d*> (_vertices[1]);
        os<<v1->id()<<" "<<v2->id()<<" ";
        AFFINE_LIE m = _measurement;
        Eigen::Vector2d trans = m.translation();
        Eigen::Matrix2d genLinear = m.generalLinear();
        os<<genLinear(0,0)<<" "<<genLinear(0,1)<<" "<<genLinear(1,0)<<" "<<genLinear(1,1)<<" "<<trans[0]<<" "<<trans[1]<<" ";
        // information matrix 
        for ( int i=0; i<information().rows(); i++ )
            for ( int j=i; j<information().cols(); j++ )
            {
                os << information() ( i,j ) << " ";
            }
        os<<endl;
        return true;
    }

    virtual void computeError()
    {
        AFFINE_LIE v1 = (static_cast<VertexAffine2d*> (_vertices[0]))->estimate();
        AFFINE_LIE v2 = (static_cast<VertexAffine2d*> (_vertices[1]))->estimate();
        _error = (_measurement.inverse()*v1.inverse()*v2).log();
    }

    /*virtual void linearizeOplus()
    {
        AFFINE_LIE v1 = (static_cast<VertexAffine2d*> (_vertices[0]))->estimate();
        AFFINE_LIE v2 = (static_cast<VertexAffine2d*> (_vertices[1]))->estimate();

    }*/

    void initialEstimate(const g2o::OptimizableGraph::VertexSet& from_, g2o::OptimizableGraph::Vertex* /*to_*/) {
    VertexAffine2d *from = static_cast<VertexAffine2d*>(_vertices[0]);
    VertexAffine2d *to   = static_cast<VertexAffine2d*>(_vertices[1]);

    if (from_.count(from) > 0) {
      to->setEstimate(from->estimate() * _measurement);
    } else
      from->setEstimate(to->estimate() * _measurement.inverse());
  }

    /*virtual void setMeasurement(AFFINE2& m){
        _measurement = m;
        _inverseMeasurement = m.inverse();
      }

    virtual bool setMeasurementData(const double* d){
        _measurement = AFFINE_LIE(d[0], d[1], d[2], d[3], d[4], d[5]);
        _inverseMeasurement = _measurement.inverse();
        return true;
      }

    virtual bool getMeasurementData(double* d) const {
        Vector6d v=_measurement.toVector();
        d[0] = v[0];
        d[1] = v[1];
        d[2] = v[2];
        d[3] = v[3];
        d[4] = v[4];
        d[5] = v[5];
        return true;
      }

    virtual int measurementDimension() const {return 6;}

    virtual bool setMeasurementFromState() {
        const VertexAffine2d* v1 = static_cast<const VertexAffine2d*>(_vertices[0]);
        const VertexAffine2d* v2 = static_cast<const VertexAffine2d*>(_vertices[1]);
        _measurement = v1->estimate().inverse()*v2->estimate();
        _inverseMeasurement = _measurement.inverse();
        return true;
      }
protected:
      AFFINE2 _inverseMeasurement;*/
};

#endif