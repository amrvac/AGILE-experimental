# obj_files = $(patsubst %.f90, $(build_dir)/obj/%, $(notdir $(f90_files)))
# TODO: get rid of amrvac.h? YES!

obj_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/obj/%.o)
local_obj_files := $(local_f90_files:$(local_build_dir)/f90/%.f90=$(local_build_dir)/obj/%.o)

compile_flags += -I$(build_dir)/f90

$(build_dir)/f90/amrvac.h:
	@touch $@

# The GNU Fortran compiler has a -J flag to place the .mod files,
# but this is not supported by other compilers. Running the compile
# command in the f90 directory should fix this.

$(build_dir)/obj/%.o: $(build_dir)/f90/%.f90 $(build_dir)/f90/amrvac.h
	@mkdir -p $(@D)
	@echo -e "Compiling $(_magenta)$(notdir $<)$(_reset)"
	@cd $(build_dir)/f90; $(compile) $(compile_flags) $< -o $@

$(local_build_dir)/obj/%.o: $(local_build_dir)/f90/%.f90 $(build_dir)/f90/amrvac.h $(obj_files)
	@mkdir -p $(@D)
	@echo -e "Compiling $(_magenta)$<$(_reset)"
	@cd $(local_build_dir)/f90; $(compile) $(compile_flags) $< -o $@
