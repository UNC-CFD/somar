# -*- Mode: Makefile -*-
export makefiles += examples_somar

# This makefile produces an executable for each
# name in the `ebase' variable
ebase := somar

# Don't mess with this block
CHOMBO_HOME = ../Chombo/lib
MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)
TRAPFPE=FALSE
include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config
LibNames :=  AMRElliptic AMRTools AMRTimeDependent BoxTools

# If you want to ignore a directory, add it to the list.
base_dir =  .
src_dirs = $(shell find ../* -maxdepth 256 -type d \
             -not \( -name .svn   -prune \) \
             -not \( -name Chombo -prune \) \
             -not \( -name exec   -prune \) \
             -not \( -name tests  -prune \) \
             -not \( -name post   -prune \) \
             -not \( -name docs   -prune \) \
             | tr '\n' ' ')
inc_dirs = $(CHOMBO_HOME)/include $(src_dirs)
cxxcppflags += $(addprefix -I,$(inc_dirs)) -I/usr/include

# INPUT = inputs.debug.iceman

include $(CHOMBO_HOME)/mk/Make.example

# .phony: idb
# idb:
# 	$(ECHO)err=0;$(foreach i,$(ebase),echo "Running Intel debugger on $i for configuration $(config) ...";idb $(foreach i,$(inc_dirs),-directory $i) -directory $(CHOMBO_HOME)/src -args ./$(_app_configs) $(RUNFLAGS) $(VERBOSE) $(INPUT);stat=$$?;echo "... $i finished with status $$stat";echo "";if [ $$stat -ne 0 ]; then err=$$stat ; fi ; ) exit $$err

