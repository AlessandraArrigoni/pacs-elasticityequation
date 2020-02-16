# pacs-elasticityequation
## A finite element method for scalar and vectorial interface problems.
### Authors: Alessandra Arrigoni and Sara Francesca Pichierri.

This repository contains the code to solve scalar elliptic problems and the 2D linear elasticity equation on a domain with a single interface conforming to the mesh using a Petrov-Galerkin finite element method. Two different algebraic formulations are implemented, namely the "Basic" and the "Symmetric" methods.

The code was developed and run on Ubuntu 18.04.2 and Ubuntu 18.04.4 with the `g++` compiler, version 7.4.0.
It is based on the library GetFEM++, version 5.3 which can be downloaded from the [official website](http://getfem.org/index.html "getFEM website") and installed following the instructions provided by the developers.

Each subfolder associated to combination problem/method is organised as follows:`inputData` contains the input files to perform the tests, `output_vtk` will store the mesh and the solutions in the `.vtk` format, `outputData` will store the matrices with the `.mm` extension and the files with the errors in the L^2 and H^1 norms. It also contains a `main_test.cc` and a `Makefile` to compile, according to the `C++14` standard, the corresponding problem/method with the following options:
* `all` : builds the program with few output information at runtime and generating only the `.vtk` files;
* `clean` : as usual, deleting also all the output files in `output_vtk` and `outputData`;
* `debug` : builds the program with extensive output information at runtime and generating also the `.mm` and errors files;
* `optimised`: as for `all` with the `-O3` flag active;
* `test` : as for `optimised` generating also the `.mm` and errors files.

#### Compile
1. Choose the desired combination problem/method and go to the associated directory (es. `linearElasticityBasic`)
2. In `Makefile` change `GETFEM_PATH` with the location of `GetFEM++` code
3. In `Makefile` change `GETFEM_LIB` and `LDLIBS` if you chose to install the `GetFEM++`, LAPACK, BLAS, Qhull libraries in locations different than the default one
4. Save and run `make` with the desired option.

#### Run tests
1. Choose the desired combination problem/method and go to the associated directory (es. `linearElasticityBasic`)
2. Compile with the chosen option for `make`
3. Select the test in `inputData` folder and modify the parameters `spatialDiscretizationX` and `spatialDiscretizationY` to set different grid refinement levels
4. Save and run the executable with `-f inputData/<testname>`. If no file is provided, the code runs with the default test `quadratic`.

To execute the tests presented in the report, type:
* in `laplacianBasic` folder : `./LaplacianBasic -f inputData/periodicDiff` with `spatialDiscretizationX = 80` and `spatialDiscretizationY = 160`
* in `laplacianSymmetric` folder : `./LaplacianSym -f inputData/periodicq1` with `spatialDiscretizationX = 80` and `spatialDiscretizationY = 160`
* in `linearElasticityBasic` folder : `./ElasticityBasic -f inputData/cubicq0` with `spatialDiscretizationX = 80` and `spatialDiscretizationY = 160`
* in `linearElasticitySymmetric` folder : `./ElasticitySym -f inputData/quadraticq1` with `spatialDiscretizationX = 80` and `spatialDiscretizationY = 160`

#### Write input files
Here are some notes and tips on how to write an input file that can be read correctly by the `LifeV::Parser` employed in the project:
* for the "basic" method, the parameter `u_BC` for the Dirichlet boundary conditions on Omega1 describe also the interface jump q0; use a **boolean condition** with a suitable tolerance to identify the x-coordinates where to apply this jump. The product of several boolean conditions is not allowed, and the parser stops processing the file as soon as the condition is false; this is the reason why it is always at the end of the string.
* for the "symmetric" method, the jump q0 is described by an *ad hoc* parameter (`qzero`) that is defined in the same way for both subdomains, according to the chosen sign convention ("Omega1 - Omega2", in this case).
* the parameter `du_BC` refers to the Neumann boundary conditions. In our tests we decided to assign Dirichlet boundary conditions on all sides, so it is always set to 0 for Omega1 and contains the jump q1 (with the "Omega1 - Omega2" sign convention) for Omega2.
* for the **vectorial** problem the strings for the first and second component are separated by a "**;**".
* do **not** put additional blank lines at the beginning of the file or blank spaces in the definition of the strings, e.g. *correct*: `4*y*(y-1)-y*y*(x>(0.5-1e-5))` or `sin(2*pi*x);sin(2*pi*y)`, *wrong*: `4*y*(y-1) - y*y* (x > (0.5-1e-5))` and `sin(2*pi*x) ; sin(2*pi*y)`.
* do **not** write apostrophes (**'**) or other special characters in **comments** since the parser stops processing the file after encountering them.

#### Build documentation
Reference manuals for the code can be generated with Doxygen, version 1.8.17 by running `doxygen Doxyfile` in the current directory. Two different formats are available:
* `doc/html` contains the manual in the html format; to open it, point a HTML browser to the `index.html` file;
* `doc/latex` contains the LaTeX documentation; to generate the `.pdf` file, go to that folder and run `make pdf`.
