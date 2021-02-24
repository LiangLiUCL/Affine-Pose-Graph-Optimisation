# Affine Pose Graph Optimisation for Fetoscopic Mosaicking

This is an extension of [G2O](https://github.com/RainerKuemmerle/g2o) that can optimise 2D affine pose graph using Lie group theory.


## Requirements for G2O

-   C++14 compiler (CI pipeline runs with gcc, clang and MSVC)
-   cmake             <http://www.cmake.org>
-   Eigen3            <http://eigen.tuxfamily.org>

On Ubuntu / Debian these dependencies are resolved by installing the
following packages.

-   cmake
-   libeigen3-dev

### Optional requirements for G2O

-   suitesparse       <http://faculty.cse.tamu.edu/davis/suitesparse.html>
-   Qt5               <http://qt-project.org>
-   libQGLViewer      <http://www.libqglviewer.com>

On Ubuntu / Debian these dependencies are resolved by installing the
following packages.

-   libsuitesparse-dev
-   qtdeclarative5-dev
-   qt5-qmake
-   libqglviewer-dev-qt5

## Install G2O

On Ubuntu, install G2O with:

```bash
$ git clone https://github.com/RainerKuemmerle/g2o
$ cd g2o
$ mkdir build
$ cd build
$ cmake ../
$ make
```

## Compilation

```bash
$ mkdir build
$ cd build
$ cmake ../
$ make
```


## Usage

###Optimisation:
```
$ affine_g2o ../data/*.g2o (* is the filename)
```

### Simulation data creation
```
sim_generate -o *.g2o (* is the filename)
```

##TODO

Add comments to the code.
Add more explanation to this README file.




