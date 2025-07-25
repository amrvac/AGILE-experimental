# dep_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/dep/%.dep)

# TODO: is there alternative to  fortdepend without dependencies?
ifeq (, $(shell which fortdepend))
$(error "fortdepend not found. Check the readme, or pip install fortdepend.")
endif

# make sure that config is read first
ifdef CONFIG_READ
$(build_dir)/dependencies.mk: $(f90_files) $(build_dir)/f90/amrvac.h | $(build_dir)
	@echo "Regenerating depencies"
	@fortdepend -f $(f90_files) -i mpi openacc -b $(build_dir)/obj -w -o $@

# Precompiling and dependency tracking is not needed if we're cleaning.
ifndef disable_precompile
include $(build_dir)/dependencies.mk
endif
endif

