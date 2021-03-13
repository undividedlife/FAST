FC := gfortran
FOPTS := -g -O2 -I./ -I../../../src
RM := /bin/rm
#
#.PHONY: clean
#
EXE := Bisect_Example2.x
SRC := $(EXE:.x=.F90)
OBJ := $(EXE:.x=.o)
#
#
all: $(SRC) $(OBJ)
	$(FC) $(FOPTS) -o $(EXE) $(OBJ) -L../../../src -lfast
clean:
	$(RM) -f $(EXE) $(OBJ) fort.*
#
$(OBJ): $(SRC)
	$(FC) $(FOPTS) -c $<
