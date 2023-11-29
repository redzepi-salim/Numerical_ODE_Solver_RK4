# Numerical_ODE_Solver_RK4

This is a numeric solver for higher order ordinary differential equations based on Runge Kutta 4. Created for Uni Class

## Dependencies

### Pugixml

Installation should be done through VCPKG, following the standard procedure.

### Matplotlibcpp

The code uses apart from standard C++ libraries and pugixml, matplotlib-cpp for graphic representations, which can be found on the following link.

[https://github.com/lava/matplotlib-cpp](https://github.com/lava/matplotlib-cpp)

### Installation

Installation instructions for Ubuntu are given in the README.md in the matplotlib-cpp repository. Because the library is essentially a wapper of a python library, python needs to be integrated to run it. Here are instructions for setting up environment in Visual Studio 2022, Python 3.10. (Windows)

1. Install the latest version of Python making sure to add check the box that adds Python to path (refers to environment variables). Make sure that python can be accessed via command prompt:

python --version

2. Install numpy and matplotlib using commands:

python -m pip install matplotlib

python -m pip install numpy

3. After creating a new project in VS (make sure to choose release and not debug, also in x64) and downloading matplotlibcpp.h from the repository add it to the project (the source files).

4. In order to add additional include directories, go to project>properties>C/C++>Additional include dirrectories and add the path to python and numpy include directories. In my case that was

> C:\Users\Korisnik\AppData\Local\Programs\Python\Python310\include;
> C:\Users\Korisnik\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\core\include

5. After that go to project>properties>linker>General>Additional library directories and add python libs. In my case that was

> C:\Users\Korisnik\AppData\Local\Programs\Python\Python310\libs

6. Then go to project>properties>linker>input>additional dependencies and add python310.lib making sure not to delete %. In my case that was

> C:\Users\Korisnik\AppData\Local\Programs\Python\Python310\libs\python310.lib;$(CoreLibraryDependencies);%(AdditionalDependencies)
