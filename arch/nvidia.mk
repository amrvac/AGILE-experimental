arch := nvidia

compile = mpif90
#compile = nvfortran
f90_flags += -cpp -Mfree

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -fopenmp
endif

ifdef OPENACC
$(info Enabling OpenACC)
#f90_flags += -Wall -acc=gpu -Minfo=all -Mvect=levels:5 -Minline
f90_flags += -Wall -acc=wait -gpu=debug -Minfo=accel -Mvect=levels:5 -Minline -Mbounds -g -O0
#f90_flags += -Wall -Minline -Mbounds -g -O0
enabled += OPENACC
endif

#ifdef DEBUG
#$(info Enable debugging symbols)
#enabled += DEBUG
#f90_flags += -g -O0
#else
#f90_flags += -O3 -fast
#endif

f90_flags += -g -O0

link_flags += $(f90_flags) -lnvToolsExt

