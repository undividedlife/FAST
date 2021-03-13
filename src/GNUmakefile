FC := gfortran
FOPTS := -g -O2 -I./
LIBST := libfast.a
#
RM := rm
AR := ar
CP := cp
#
.PHONY: all clean
#
SRCS := Base.F90 PhysConst.F90 Bisect.F90 Tricub.F90 Interp.F90
LIBS := $(SRCS:.F90=.o)
MODS := $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)
#
all: $(LIBS) $(SRCS) $(MODS)
	$(AR) rcs $(LIBST) $(LIBS)
clean:
	$(RM) -f $(LIBS) $(MODS) $(LIBST)
#
# Dependencies
#
Base.o: Base.F90
	$(FC) -c $(FOPTS) $<
PhysConst.o: PhysConst.F90
	$(FC) -c $(FOPTS) $<
Bisect.o: Bisect.F90 Base.o
	$(FC) -c $(FOPTS) $<
Tricub.o: Tricub.F90 Base.o
	$(FC) -c $(FOPTS) $<
Interp.o: Interp.F90 Base.o Bisect.o Tricub.o
	$(FC) -c $(FOPTS) $<
