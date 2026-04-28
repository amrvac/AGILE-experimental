$(build_dir)/obj/amrvac.a: $(obj_files)
	@echo -e "Linking $(_green)$(notdir $@)$(_reset)"
	@$(link) $(link_flags) -static $^ -o $@

amrvac: $(local_build_dir)/amrvac
	@rm -f amrvac
	@ln -s $< $@

$(local_build_dir)/amrvac: $(local_obj_files) $(build_dir)/obj/amrvac.a
	@echo -e "Linking $(_green)$(notdir $@)$(_reset)"
	@$(link) $(link_flags) $^ -o $@
