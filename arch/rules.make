VACPP := $(AMRVAC_DIR)/src/vacpp.pl

# Disable built-in make rules
.SUFFIXES:


# How to compile modules into object files
mod_%.o: mod_%.F90 amrvac.h
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to get .mod files from modules. Modules are automatically updated, and
# only need to be explicitly generated when they were manually removed.
mod_%.mod: mod_%.F90 mod_%.o
	@test -f $@ || $(F90) $(F90FLAGS) -c $(@:.mod=.F90) -o $(@:.mod=.o) $(addprefix -I,$(INC_DIRS))
m_%.mod: m_%.F90 m_%.o
	@test -f $@ || $(F90) $(F90FLAGS) -c $(@:.mod=.F90) -o $(@:.mod=.o) $(addprefix -I,$(INC_DIRS))

# How to generate object files
%.o: %.F90 amrvac.h
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# How to translate .t source files to normal Fortran
%.F90: %.t
	$(VACPP) $(PPFLAGS) -d=$(NDIM) $< > $(@)

# How to generate executables
%: %.o
	$(LINK) $(F90FLAGS) $^ -o $@ $(addprefix -L,$(LIB_DIRS)) $(addprefix -l,$(LIBS))


