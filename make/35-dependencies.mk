# dep_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/dep/%.dep)

# TODO: is there alternative to  fortdepend without dependencies?
ifeq (, $(shell which fortdepend))
$(error "fortdepend not found. Check the readme, or pip install fortdepend.")
endif

ifdef USE_MPIWRAPPERS
  fortdepend_flags += -DUSE_MPIWRAPPERS
endif
ifdef OPENACC
  fortdepend_flags += -D_OPENACC
endif

# make sure that config is read first
ifdef CONFIG_READ
$(build_dir)/dependencies.mk: $(f90_files) $(build_dir)/f90/amrvac.h | $(build_dir)
	@echo "Regenerating depencies"
	@fortdepend $(fortdepend_flags) -s -f $(f90_files) -i mpi openacc -b $(build_dir)/obj -w -o $@

$(local_build_dir)/dependencies.mk: $(local_f90_files) | $(local_build_dir)
	@echo "Regenerating local depencies"
	@fortdepend $(fortdepend_flags) -s -f $(local_f90_files) -i mpi openacc -b $(local_build_dir)/obj -w -o $@

# Precompiling and dependency tracking is not needed if we're cleaning.
ifndef disable_precompile
include $(build_dir)/dependencies.mk
endif
endif
