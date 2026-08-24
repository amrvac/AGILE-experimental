# Optional external libraries enabled through config flags.

ifneq ($(filter -DCOMPRESS_ZFP%,$(fypp_flags)),)

ifndef ZFP_DIR
$(error ZFP_DIR is not set; export ZFP_DIR=?/zfp/install)
endif

$(info Enabling ZFP compression, dir: $(ZFP_DIR))
f90_flags  += -I$(ZFP_DIR)/include
link_flags += -L$(ZFP_DIR)/lib64 -lzFORp -lzfp -Wl,-rpath,$(ZFP_DIR)/lib64

endif
