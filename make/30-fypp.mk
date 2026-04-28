ifdef DEBUG
fypp_flags += -DDEBUG
endif

ifeq (, $(shell which fypp))
$(error "fypp not found. Check the readme, or pip install fypp.")
endif

fypp_flags += -n

fypp_flags += -M $(amrvac)/make

$(info Fypp flags: $(fypp_flags))

local_source_files := mod_usr.fpp
local_f90_files := $(patsubst %.fpp, $(local_build_dir)/f90/%.f90, \
		$(notdir $(local_source_files)))

source_files := $(shell find $(amrvac)/src -name '*.fpp')
f90_files := $(patsubst %.fpp, $(build_dir)/f90/%.f90, \
	     $(notdir $(source_files)))

.PRECIOUS: $(f90_files)

define fypp_rule
-include $(patsubst %.fpp, $(build_dir)/f90/%.d, $(notdir $(1)))

$(patsubst %.fpp, $(build_dir)/f90/%.f90, $(notdir $(1))): $(1)
	@mkdir -p $$(@D)
	@echo -e "Precompiling $(_magenta)$$(<:$(amrvac)/src/%=%)$(_reset)"
	@python $(amrvac)/make/fypp-deps.py $$< $$@
	@fypp $(fypp_flags) $$<  $$@
endef

define local_fypp_rule
-include $(patsubst %.fpp, $(local_build_dir)/f90/%.d, $(notdir $(1)))

$(patsubst %.fpp, $(local_build_dir)/f90/%.f90, $(notdir $(1))): $(1)
	@mkdir -p $$(@D)
	@echo -e "Precompiling $(_magenta)./$<$(_reset)"
	@python $(amrvac)/make/fypp-deps.py $$< $$@
	@fypp $(fypp_flags) $$<  $$@
endef

$(foreach f, $(source_files), $(eval $(call fypp_rule, $(f))))
$(foreach f, $(local_source_files), $(eval $(call local_fypp_rule, $(f))))
