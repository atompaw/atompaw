GOBJS =	aeatom.o anderson_driver.o atomdata.o blockdavidson_mod.o \
	numerov_mod.o excor.o exxdata.o exx_mod.o fock.o general_mod.o \
	globalmath.o hf_mod.o\
	graphatom.o graphatom_report.o gridmod.o ldagga_mod.o radialsr.o \
	report_mod.o search_sort.o

graphatom:	$(GOBJS)
	$(F90) $(LDFLAGS) -o $@ $(GOBJS) $(LIBS)
