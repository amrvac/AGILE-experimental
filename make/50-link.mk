obj_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/obj/%.o)
local_obj_files := $(local_f90_files:$(local_build_dir)/f90/%.f90=$(local_build_dir)/obj/%.o)

$(build_dir)/obj/amrvac.a: $(obj_files)
	@echo -e "Linking $(_green)$(notdir $@)$(_reset)"
	@$(link) $(link_flags) -static $^ -o $@

amrvac: $(local_build_dir)/obj/amrvac
	@rm -f amrvac
	@ln -s $< $@

$(local_build_dir)/obj/amrvac: $(local_obj_files) $(build_dir)/obj/amrvac.a
	@echo -e "Linking $(_green)$(notdir $@)$(_reset)"
	@$(link) $(link_flags) $^ -o $@
