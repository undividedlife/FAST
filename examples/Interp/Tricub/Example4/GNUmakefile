FC := gfortran
FCOPTS := -g -O2 -I./ -I../../../../src
RM := /bin/rm
#
EXE := Tricub_Example4.x
SRC := $(EXE:.x=.F90)
OBJ := $(EXE:.x=.o)
#
all: $(OBJ)
	$(FC) $(FCOPTS) -o $(EXE) $(OBJ) -L../../../../src -lfast
clean:
	$(RM) -f $(EXE) $(OBJ)
#
$(OBJ): $(SRC)
	$(FC) $(FCOPTS) -c $<
