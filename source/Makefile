# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

sim_dsam:		sim_dsam.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

atomicWeights:		atomicWeights.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

decay:			decay.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

depthProfile:		depthProfile.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

detectionEfficiency:	detectionEfficiency.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

globals:		globals.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

massTable:		massTable.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

multipleScattering:	multipleScattering.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

particle:		particle.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

sequence:		sequence.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

stoppingDataTables:	stoppingDataTables.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"

rootDataContainer:	rootDataContainer.o
			$(LD) $(LDFLAGS) $^ $(LIBS) -lMathMore $(OutPutOpt)$@
			@echo "$@ done"


clean:
			@rm -f $(OBJS) core sim_dsam.o atomicWeights.o decay.o depthProfile.o detectionEfficiency.o globals.o massTable.o multipleScattering.o particle.o sequence.o stoppingDataTables.o rootDataContainer.o
 
distclean:      	clean
			-@mv -f linearIO.root linearIO.roott
			@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp 
			   *.root *.ps *.so *.lib *.dll *.d .def so_locations
			@rm -rf cxx_repository
			-@mv -f linearIO.roott linearIO.root
			-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
