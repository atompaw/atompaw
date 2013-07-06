PROG =	atompaw

SRCS =	abinitinterface.f90 aeatom.f90 anderson_driver.f90 atomdata.f90 \
	atompaw_lib.f90 atompaw_prog.f90 atompaw_report.f90 blockdavidson_mod.f90 \
	excor.f90 exx_mod.f90 exx_pseudo.f90 exxdata.f90 fock.f90  \
	general_mod.f90 globalmath.f90 gridmod.f90 hf_mod.f90 hf_pseudo.f90 \
	interpolation_mod.f90 ldagga_mod.f90 numerov_mod.f90 paw_sub.f90 \
	pkginfo.f90 pseudo.f90 pseudo_sub.f90 pseudodata.f90 pwscfinterface.f90 \
	radialsr.f90 report_mod.f90 search_sort.f90 libxc_mod.F90 xmlinterface.f90


OBJS =	abinitinterface.o aeatom.o anderson_driver.o atomdata.o \
	atompaw_lib.o atompaw_prog.o atompaw_report.o blockdavidson_mod.o \
	excor.o exx_mod.o exx_pseudo.o exxdata.o fock.o \
	general_mod.o globalmath.o gridmod.o hf_mod.o hf_pseudo.o \
	interpolation_mod.o ldagga_mod.o numerov_mod.o paw_sub.o \
	pkginfo.o pseudo.o pseudo_sub.o pseudodata.o pwscfinterface.o \
	radialsr.o report_mod.o search_sort.o libxc_mod.o xmlinterface.o

include make.inc

all:	$(PROG)

$(PROG):	$(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

include graphatom.inc

compile:
	$(F90) $(LDFLAGS) -o $(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod
	rm -f graphatom graphatom.o

.SUFFIXES: $(SUFFIXES) .f90

%.o :	%.f90
	$(F90) $(F90FLAGS) -c $<

%.o :	%.F90
	$(F90) $(F90FLAGS) -c $<

abinitinterface.o: atomdata.o globalmath.o gridmod.o excor.o \
	pseudo.o interpolation_mod.o pkginfo.o libxc_mod.o
aeatom.o: atomdata.o excor.o exx_mod.o general_mod.o globalmath.o gridmod.o \
	hf_mod.o ldagga_mod.o
anderson_driver.o: globalmath.o
atomdata.o:
atompaw_lib.o:
atompaw_prog.o: globalmath.o aeatom.o atomdata.o atompaw_report.o \
	pseudo.o abinitinterface.o pwscfinterface.o xmlinterface.o pkginfo.o \
	libxc_mod.o
atompaw_report.o: atomdata.o fock.o globalmath.o gridmod.o pseudo.o pseudodata.o \
	libxc_mod.o
blockdavidson_mod.o: globalmath.o search_sort.o
excor.o: atomdata.o globalmath.o gridmod.o libxc_mod.o
exx_mod.o: anderson_driver.o atomdata.o exxdata.o fock.o general_mod.o \
	globalmath.o gridmod.o hf_mod.o report_mod.o
exx_pseudo.o: anderson_driver.o atomdata.o exx_mod.o fock.o globalmath.o \
	gridmod.o pseudodata.o pseudo_sub.o
exxdata.o:
fock.o: atomdata.o gridmod.o globalmath.o
general_mod.o: atomdata.o numerov_mod.o globalmath.o gridmod.o radialsr.o
globalmath.o:
graphatom.o: globalmath.o aeatom.o atomdata.o graphatom_report.o
graphatom_report.o: atomdata.o gridmod.o report_mod.o
gridmod.o: globalmath.o
hf_mod.o: anderson_driver.o atomdata.o fock.o general_mod.o \
	globalmath.o gridmod.o report_mod.o
hf_pseudo.o: atomdata.o hf_mod.o fock.o globalmath.o \
	gridmod.o paw_sub.o pseudodata.o pseudo_sub.o
interpolation_mod.o:
ldagga_mod.o: atomdata.o anderson_driver.o excor.o general_mod.o globalmath.o \
	gridmod.o report_mod.o
libxc_mod.o: globalmath.o
numerov_mod.o: gridmod.o blockdavidson_mod.o
paw_sub.o: globalmath.o gridmod.o
pkginfo.o:
pseudo.o: globalmath.o atomdata.o aeatom.o excor.o exx_pseudo.o hf_pseudo.o \
	numerov_mod.o paw_sub.o pseudodata.o pseudo_sub.o radialsr.o
pseudo_sub.o: globalmath.o gridmod.o pseudodata.o excor.o fock.o
pseudodata.o: gridmod.o atomdata.o
pwscfinterface.o: globalmath.o atomdata.o excor.o gridmod.o \
	interpolation_mod.o pseudo.o
radialsr.o: globalmath.o gridmod.o atomdata.o
report_mod.o: atomdata.o gridmod.o excor.o
search_sort.o:
xmlinterface.o: atomdata.o globalmath.o gridmod.o excor.o pseudo.o pkginfo.o
