ATOMPAW
================

The software `ATOMPAW` generates projector and basis functions which are needed for
performing electronic structure calculations (Density-Functional Theory)
based on the Projector Augmented-Wave (PAW) method.

### What does AtomPAW?

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


If you obtained the sources directly from the git repository,
you will first need to generate the configure script by running
```
./bootstrap.sh
```
(Not needed if you downloaded the sources from ATOMPAW
[website](http://users.wfu.edu/natalie/papers/pwpaw))  

Then run:  
```
./configure --prefix=PATH/TO/ATOMPAW [options]
make
make install
```
Most common options (complete list: `./configure --help`):
- A `blas/lapack` library is required. If not present in a standard directory, use:  
  `--with-linalg-prefix=PATH/TO/LINEAR/ALGEBRA`
- To link with [libxc](https://www.tddft.org/programs/libxc/) collection of 
  exchange-correlation functionals, use:  
  `--enable-libxc --with-libxc-prefix=PATH/TO/LIBXC`.
