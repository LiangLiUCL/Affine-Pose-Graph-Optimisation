#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h> 
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include "AFFINE_LIE.h"
#include "AFFINE_VERTEX.h"
#include "AFFINE_EDGE.h"

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
#include "g2o/stuff/sampler.h"
#include "g2o/stuff/command_args.h"
#include "g2o/core/factory.h"

using namespace std;
using namespace g2o;

typedef Eigen::Matrix<double, 4, 1, Eigen::ColMajor> Vector4d;
typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> Vector6d;


int main(int argc, char** argv)
{
  // command line parsing
  int nodesPerLevel;
  int numLaps;
  double radius;
  std::vector<double> noiseTranslation;
  std::vector<double> noiseRotation;
  std::vector<double> noiseScaling;
  string outFilename;
  CommandArgs arg;
  arg.param("o", outFilename, "-", "output filename");
  arg.param("nodesPerLevel", nodesPerLevel, 50, "how many nodes per lap on the sphere");
  arg.param("laps", numLaps, 5, "how many times the robot travels around the sphere");
  arg.param("radius", radius, 100., "radius of the sphere");
  arg.param("noiseTranslation", noiseTranslation, std::vector<double>(), "set the noise level for the translation, separated by semicolons without spaces e.g: \"0.1;0.1;0.1\"");
  arg.param("noiseRotation", noiseRotation, std::vector<double>(), "set the noise level for the rotation, separated by semicolons without spaces e.g: \"0.001;0.001;0.001\"");
  arg.parseArgs(argc, argv);

  if (noiseTranslation.size() == 0) {
    cerr << "using default noise for the translation" << endl;
    noiseTranslation.push_back(0.03);//(0.02);//(0.01);
    noiseTranslation.push_back(0.03);//(0.02);//(0.01);
  }
  cerr << "Noise for the translation:";
  for (size_t i = 0; i < noiseTranslation.size(); ++i)
    cerr << " " << noiseTranslation[i];
  cerr << endl;
  if (noiseRotation.size() == 0) {
    cerr << "using default noise for the rotation" << endl;
    noiseRotation.push_back(0.010);//(0.008);//(0.005);
    noiseRotation.push_back(0.010);//(0.008);//(0.005);
    noiseRotation.push_back(0.010);//(0.008);//(0.005);
    noiseRotation.push_back(0.010);//(0.008);//(0.005);
  }
  cerr << "Noise for the rotation:";
  for (size_t i = 0; i < noiseRotation.size(); ++i)
    cerr << " " << noiseRotation[i];
  cerr << endl;

  Eigen::Matrix2d transNoise = Eigen::Matrix2d::Zero();
  for (int i = 0; i < 2; ++i)
    transNoise(i, i) = std::pow(noiseTranslation[i], 2);

  Eigen::Matrix4d rotNoise = Eigen::Matrix4d::Zero();
  for (int i = 0; i < 4; ++i)
    rotNoise(i, i) = std::pow(noiseRotation[i], 2);

  Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();
  information.block<2,2>(4,4) = transNoise.inverse();
  information.block<4,4>(0,0) = rotNoise.inverse();

  vector<VertexAffine2d*> vectices;
  vector<EdgeAffine2d*> odometryEdges;
  vector<EdgeAffine2d*> edges;
  vector<VertexAffine2d*> vecticesNoise;
  int id = 0;
  double scale = 1.;
  for (int f=0; f<numLaps; f++)
  {
    for (int n=0; n<nodesPerLevel; n++)
    {
      VertexAffine2d* v = new VertexAffine2d;
      v->setId(id++);
      Eigen::Rotation2D<double> rot(-M_PI + 2*n*M_PI / nodesPerLevel);
      scale = 1 - 0.7*id/(numLaps * nodesPerLevel);

      Eigen::Affine2d affIni;
      affIni = Eigen::Scaling(scale)*rot;
      affIni.translation() = affIni.linear() * Eigen::Vector2d(radius, 0);
      AFFINE_LIE affInp(affIni);
      v->setEstimate(affInp);
      vectices.push_back(v);
    }
  }


  //generate odometry edges
  for (size_t i=1; i<vectices.size(); ++i)
  {
    VertexAffine2d* prev = vectices[i-1];
    VertexAffine2d* cur  = vectices[i];
    AFFINE_LIE t = prev->estimate().inverse() * cur->estimate();
    EdgeAffine2d* e = new EdgeAffine2d;
    e->setVertex(0, prev);
    e->setVertex(1, cur);
    e->setMeasurement(t);
    e->setInformation(information);
    odometryEdges.push_back(e);
    edges.push_back(e);
    Vector6d toAlg = e->measurement().log();
    AFFINE_LIE noisyMeasurement(toAlg);
    //cout << noisyMeasurement.generalLinear() << endl;
  }


  //generate loop closure edges
  for (int f = 1; f < numLaps; ++f)
  {
    for (int nn = 0; nn < nodesPerLevel; ++nn)
    {
      VertexAffine2d* from = vectices[(f-1)*nodesPerLevel + nn];
      for (int n = -1; n <= 1; ++n)
      {
        if (f == numLaps-1 && n == 1)
          continue;
        VertexAffine2d* to   = vectices[f*nodesPerLevel + nn + n];
        AFFINE_LIE t = from->estimate().inverse() * to->estimate();
        EdgeAffine2d* e = new EdgeAffine2d;
        e->setVertex(0, from);
        e->setVertex(1, to);
        e->setMeasurement(t);
        e->setInformation(information);
        edges.push_back(e);
      }
    }
  }

  GaussianSampler<Eigen::Vector2d, Eigen::Matrix2d> transSampler;
  transSampler.setDistribution(transNoise);
  GaussianSampler<Vector4d, Eigen::Matrix4d> rotSampler;
  rotSampler.setDistribution(rotNoise);

  //noise for all edges
  for (size_t i = 0; i < edges.size(); ++i)
  {
    EdgeAffine2d* e = edges[i];
    Eigen::Matrix2d gtLin = e->measurement().generalLinear();
    Eigen::Vector2d gtTrans = e->measurement().translation();

    Vector4d linN = rotSampler.generateSample();
    Eigen::Vector2d transN = transSampler.generateSample();
    Vector6d toAlg = e->measurement().log();
    Vector6d noiseAlg;
    noiseAlg.block<4, 1>(0, 0) = linN;
    noiseAlg.block<2, 1>(4, 0) = transN;
    toAlg = toAlg + noiseAlg;
    //cout << toAlg.transpose() << endl;
    AFFINE_LIE noisyMeasurement(toAlg);
    e->setMeasurement(noisyMeasurement);
  }



  AFFINE_LIE accVertex(-1, 0, 0, -1, -100, 0);
  VertexAffine2d* v1 = new VertexAffine2d;
  v1->setId(0);
  v1->setEstimate(accVertex);
  vecticesNoise.push_back(v1);

  //update the vectices based on the odometry edge
  for (size_t i = 0; i < vectices.size()-1; i++)
  {
    EdgeAffine2d* edgeI = edges[i];
    AFFINE_LIE edgeIAff = edgeI->measurement();

    accVertex = accVertex * edgeIAff;

    VertexAffine2d* v = new VertexAffine2d;
    v->setId(i+1);
    v->setEstimate(accVertex);
    vecticesNoise.push_back(v);
  }

  // write output
  ofstream fileOutputStream;
  if (outFilename != "-") {
    cerr << "Writing into " << outFilename << endl;
    fileOutputStream.open(outFilename.c_str());
  } else {
    cerr << "writing to stdout" << endl;
  }


  ostream& fout = outFilename != "-" ? fileOutputStream : cout;
  for (size_t i = 0; i < vecticesNoise.size(); ++i) {
    VertexAffine2d* v = vecticesNoise[i];
    fout << "Vertex" << " " << v->id() << " ";
    v->write(fout);
  }

  for (size_t i = 0; i < edges.size(); ++i) {
    EdgeAffine2d* e = edges[i];
    VertexAffine2d* from = static_cast<VertexAffine2d*>(e->vertex(0));
    VertexAffine2d* to = static_cast<VertexAffine2d*>(e->vertex(1));
    fout << "Edge" << " ";// << from->id() << " " << to->id() << " ";
    e->write(fout);
  }

  return 0;

}