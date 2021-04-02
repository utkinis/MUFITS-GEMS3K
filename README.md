# USEREOS
[![DOI](https://zenodo.org/badge/334063009.svg)](https://zenodo.org/badge/latestdoi/334063009)

Code for the paper "Using user-supplied modules for fluid properties prediction with the MUFITS reservoir simulator"

## Installation and usage
To run the examples presented in the paper, you need the MUFITS simulator executable file, `H64.EXE`, and a user-supplied module in the form of a shared library with name `USEREOS.DLL`. There are three ways to get the `USEREOS` modules:

1. Get the pre-built binaries and source code from the [USEREOS section of MUFITS website](http://www.mufits.imec.msu.ru/example-usereos.html);
2. Download pre-built binaries and source code from [Github releases](https://github.com/utkinis/USEREOS/releases);
3. Build the shared libraries directly from the source code.

MUFITS executables are available for 64-bit Windows, OS X and Debian Linux. Pre-built USEREOS shared libraries are available only for Windows. For other platforms, it is necessary to build USEROS modules from source.

## Build instructions
Building the modules from source requires Fortran and C++ compilers. Examples 1 and 2 don't have any dependencies. Example 3 depends on GEMS3K, numerical code for calculating chemical equilibria. GEMS3K can be downloaded from its [homepage](http://gems.web.psi.ch/GEMS3K/).

To access the latest version of the code, you need to either download the source code as an archive or clone the repository. Cloning the repository requires working [Git](https://git-scm.com/) installation.

To clone the repository, open the command-line terminal and enter the following command:

```
> git clone https://github.com/utkinis/USEREOS.git
```
