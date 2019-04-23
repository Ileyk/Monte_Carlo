# TO DO
# 1. What abour if not all my files are modules? And what is not a module?
# 2. How can I use files stored in a src folder to produce an executable
#    in the folder where the makefile is run?
# 3. How to make clean automatically in the end?
# 4. Make all? Make parallel?

# Beware the order of the objects! The most fundamental first and the main last
# It sounds like the order should be the same as when we compile directly.
objects = glbl_prmtrs.o miscellaneous.o mod_kracken.o IO.o mc.o

ifndef MC_DIR
$(error MC_DIR is not set)
endif

# Here, I read the compiling options
# The user should specify the compiling mode here, and the corresponding
# set of compiling options will be read from the file $(MC_DIR)/arch/
ARCH ?= default
export ARCH
# Include compilation rules
include $(MC_DIR)/arch/$(ARCH).defs

mc: $(objects)
	$(F90) $(F90FLAGS) -o $@ $^
%.mod:%.90 %.o
	$(F90) $(F90FLAGS) -c -o $@ $^
%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# In case there is a file named "clean" in the folder, it will not be used
# as a prerequisite file
# https://stackoverflow.com/questions/2145590/what-is-the-purpose-of-phony-in-a-makefile
.PHONY: clean

clean:
	rm *.o *.mod
