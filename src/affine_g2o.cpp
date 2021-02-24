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


Vector6d affineVec(Eigen::Affine2d errMatrix)
{
    double x = errMatrix.matrix()(0, 2);
    double y = errMatrix.matrix()(1, 2);
    double theta = atan2(-errMatrix.matrix()(0, 1) + errMatrix.matrix()(1, 0), errMatrix.matrix()(0, 0) + errMatrix.matrix()(1, 1));
    Eigen::Rotation2D<double> rotTheta(theta);
    Eigen::Affine2d affDec = rotTheta.inverse()*errMatrix;
    double phi = atan2(affDec.matrix()(0, 1) + affDec.matrix()(1, 0), affDec.matrix()(1, 1) - affDec.matrix()(0, 0));
    phi = phi/2;
    double s_x = affDec.matrix()(0, 0) + affDec.matrix()(1, 1) - (affDec.matrix()(1, 1) - affDec.matrix()(0, 0))/cos(2*phi);
    s_x = s_x/2;
    double s_y = affDec.matrix()(0, 0) + affDec.matrix()(1, 1) - s_x;

    Vector6d result;
    result << x, y, theta, phi, s_x, s_y;
    return result;
}


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
        Eigen::Vector2d trans = _estimate.translation();
        Eigen::Matrix2d genLinear = _estimate.generalLinear();

        os<<genLinear(0,0)<<" "<<genLinear(0,1)<<" "<<trans[0]<<" "<<genLinear(1,0)<<" "<<genLinear(1,1)<<" "<<trans[1]<<" "<<endl;
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

};


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
        os<<genLinear(0,0)<<" "<<genLinear(0,1)<<" "<<trans[0]<<" "<<genLinear(1,0)<<" "<<genLinear(1,1)<<" "<<trans[1]<<" ";
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

};



int main(int argc, char** argv)
{

    if ( argc != 2 )
    {
        cout<<"Usage: affine_g2o odometry.g2o"<<endl;
        return 1;
    }
    ifstream fin ( argv[1] );
    if ( !fin )
    {
        cout<<"file "<<argv[1]<<" does not exist."<<endl;
        return 1;
    }

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6,6>> BlockSolverType;
    typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;
    auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));  

    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm ( solver );

    int vertexCnt = 0, edgeCnt = 0;

    vector<VertexAffine2d*> vectices;
    vector<EdgeAffine2d*> edges;
    while ( !fin.eof() )
    {
        string name;
        fin>>name;
        if ( name == "Vertex" )
        {
            VertexAffine2d* v = new VertexAffine2d();
            int index = 0;
            fin>>index;
            v->setId( index );
            v->read(fin);
            optimizer.addVertex(v);
            vertexCnt++;
            vectices.push_back(v);
            if ( index==0 )
                v->setFixed(true);
        }
        else if ( name=="Edge" )
        {
            EdgeAffine2d* e = new EdgeAffine2d();
            int idx1, idx2;
            fin>>idx1>>idx2;
            e->setId( edgeCnt++ );
            e->setVertex( 0, optimizer.vertices()[idx1] );
            e->setVertex( 1, optimizer.vertices()[idx2] );
            e->read(fin);
            optimizer.addEdge(e);
            edges.push_back(e);
        }
        if ( !fin.good() ) break;
    }

    cout<<"read total "<<vertexCnt<<" vertices, "<<edgeCnt<<" edges."<<endl;

    cout<<"prepare optimizing ..."<<endl;
    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    cout<<"calling optimizing ..."<<endl;
    optimizer.optimize(100);

    cout<<"saving optimization results ..."<<endl;
    ofstream fout("result_pgo.g2o");
    for ( VertexAffine2d* v:vectices )
    {
        fout<<"Vertex ";
        v->write(fout);
    }
    for ( EdgeAffine2d* e:edges )
    {
        fout<<"Edge ";
        e->write(fout);
    }
    fout.close();
    return 0;
}