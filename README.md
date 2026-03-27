ATOMPAW
================

The software `ATOMPAW` generates projector and basis functions which are needed for
performing electronic structure calculations (Density-Functional Theory)
based on the Projector Augmented-Wave (PAW) method.

### What does AtomPAW do?

The program is applicable to materials throughout the periodic table.
For each element, the user inputs the atomic number, the electronic configuration,
a choice of basis functions, and an augmentation radius.
The program produces output files containing the projector and basis functions
and the corresponding matrix elements in a format which can be read by several
DFT codes ([abinit](https://www.abinit.org),
[quantum expresso](https://www.quantum-espresso.org),
[gpaw](https://wiki.fysik.dtu.dk/gpaw),
[pwpaw](http://users.wfu.edu/natalie/papers/pwpaw),
[onetep](http://www.onetep.org),
and all codes that can read atomic PAW setups in the
[PAW-XML](https://esl.cecam.org/Paw-xml) format).


Most of the relevant information can be found on the
ATOMPAW [official website](http://users.wfu.edu/natalie/papers/pwpaw).

Many documentation files can be found in the doc directory.
See especially the 
[~/doc/atompaw-usersguide.pdf](https://github.com/atompaw/atompaw/blob/master/doc/atompaw-usersguide.pdf)
 file.

### License

See `COPYING file`

### Installation

#### Installation via CMake

Installation with `CMake` is the most straightforward method.
Once inside the ATOMPAW source directory, simply run:  

```
mkdir build && cd build
cmake ..
make
[make install]
```  

The linear algebra library (`BLAS`/`LAPACK`) should be automatically detected. If not, you can add the following on the `cmake ..` line:  

```
-DBLAS_ROOT=/PATH/TO/BLAS -DLAPACK_ROOT=PATH/TO/LAPACK
```

The [libxc](https://libxc.gitlab.io) library (collection of exchange-correlation functionals) should be automatically detected, although it is optional. If not, add the following on the `cmake ..` line:  

```
-DLIBXC_ROOT=/PATH/TO/LIBXC
```

#### Installation via Autotools

If you obtained the sources directly from the git repository,
you will first need to generate the `configure` script by running:  

```
./bootstrap.sh
```
(This step is not necessary if you downloaded the sources from the ATOMPAW
[website](http://users.wfu.edu/natalie/papers/pwpaw))  

Then run:  

```
mkdir build && cd build
../configure --prefix=PATH/TO/ATOMPAW [options]
make
[make install]
```  

Most commonly used options (for a complete list, run `./configure --help`):  

- A `BLAS`/`LAPACK` library is required. If not installed in a standard location, specify it with:  
  `--with-linalg-prefix=PATH/TO/LINEAR/ALGEBRA`  
- To link against the [libxc](https://www.tddft.org/programs/libxc/) library, which provides a collection of 
  exchange-correlation functionals, use:  
  `--enable-libxc --with-libxc-prefix=PATH/TO/LIBXC`.

#### Installation via Homebrew (MacOS)

If you are using macOS, you can use the [Homebrew package manager](https://brew.sh) to easily install the latest version of ATOMPAW.
Once Homebrew is installed, simply run:  

```
brew tap atompaw/repo
brew install atompaw
```

or directly:

```
brew install atompaw/repo/atompaw
```

> Notes:  
> 
> - Always use the latest Homebrew version (use brew upgrade).  
> - Bottles (compiled versions) are not always provided. If so, ATOMPAW will be built on the fly during installation process.
