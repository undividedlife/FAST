FC := gfortran
FOPTS := -g -O2 -I./ -I../../../src
RM := /bin/rm
#
.PHONY: all clean
#
EXE := Interp_Example1.x
SRC := $(EXE:.x=.F90)
OBJ := $(EXE:.x=.o)
#
all:: $(SRC) $(OBJ)
	$(FC) $(FOTPS) -o $(EXE) $(OBJ) -L../../../src -lfast
clean:
	$(RM) -f $(EXE) $(OBJ)
#
# Dependencies
#
$(OBJ): $(SRC)
	$(FC) -c $(FOPTS) $<
