CXX := g++
RM := /bin/rm
#
MAINEXE := example1.x
MAINSRC := $(shell echo $(MAINEXE:.x=.cpp))
MAINOBJ := $(shell echo $(MAINEXE:.x=.o))
#
CXXSRCS := libtricubic.cpp tricubic_utils.cpp
CXXOBJS := $(shell echo $(CXXSRCS:.cpp=.o))
#
all: $(MAINOBJ) $(CXXOBJS)
	$(CXX) -o $(MAINEXE) $(MAINOBJ) $(CXXOBJS)
clean:
	$(RM) -f $(MAINEXE) $(MAINOBJ) $(CXXOBJS)
#
$(MAINOBJ): $(MAINSRC) $(CXXOBJS)
	$(CXX) -I./ -c $<
libtricubic.o: libtricubic.cpp tricubic_utils.o
	$(CXX) -I./ -c $<
tricubic_utils.o: tricubic_utils.cpp
	$(CXX) -I./ -c $<
