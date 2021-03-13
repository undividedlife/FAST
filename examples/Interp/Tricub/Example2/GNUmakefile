FC := gfortran
CXX := g++
FCOPTS := -g -O2 -I./ -I../../../../src
RM := /bin/rm
#
EXE := Tricub_Example2.x
SRC := $(EXE:.x=.F90)
OBJ := $(EXE:.x=.o)
#
CXXSRCS := randomval.cpp
CXXOBJS := $(shell echo $(CXXSRCS:.cpp=.o))
#
all: $(OBJ) $(CXXOBJS)
	$(FC) $(FCOPTS) -o $(EXE) $(OBJ) $(CXXOBJS) -L../../../../src -lfast
clean:
	$(RM) -f $(EXE) $(OBJ) $(CXXOBJS) testmonoc.dat
#
$(OBJ): $(SRC)
	$(FC) $(FCOPTS) -c $<
randomval.o: randomval.cpp
	$(CXX) $(FCOPTS) -c $<
