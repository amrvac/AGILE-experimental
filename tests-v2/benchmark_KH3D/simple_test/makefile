# Simple Makefile for OpenACC Fortran test

# Choose your compiler:
FC       = nvfortran
#FC      = gfortran   # Uncomment if you want to try GCC OpenACC support

# Compiler flags:
FFLAGS   = -acc -Minfo=accel
# For GCC instead, use:
#FFLAGS  = -fopenacc -O2

# Target program
TARGET   = test
SRC      = test.f90

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o *.mod

